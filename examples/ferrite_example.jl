using Hyperelastics
using DifferentiationInterface
import ForwardDiff
using Ferrite, Tensors, TimerOutputs, ProgressMeter, IterativeSolvers

struct NeoHookeanParameters
    μ::Float64
    λ::Float64
end

function ψ(C, mp)
    I = [tr(C), 0.5 * (tr(C)^2 - tr(C^2)), det(C)]
    J = sqrt(I[3])
    W = StrainEnergyDensity(NeoHookean(InvariantForm()), I, mp) - mp.μ*log(J) + mp.λ/2*log(J)^2
    return W
end

function constitutive_driver(C, mp::NeoHookeanParameters)
    # Compute all derivatives in one function call
    ∂²Ψ∂C², ∂Ψ∂C = Tensors.hessian(y -> ψ(y, mp), C, :all)
    S = 2.0 * ∂Ψ∂C
    ∂S∂C = 2.0 * ∂²Ψ∂C²
    return S, ∂S∂C
end;

function assemble_element!(ke, ge, cell, cv, fv, mp, ue, ΓN)
    # Reinitialize cell values, and reset output arrays
    reinit!(cv, cell)
    fill!(ke, 0.0)
    fill!(ge, 0.0)

    b = Vec{3}((0.0, -0.5, 0.0)) # Body force
    tn = 0.1 # Traction (to be scaled with surface normal)
    ndofs = getnbasefunctions(cv)

    for qp in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, qp)
        # Compute deformation gradient F and right Cauchy-Green tensor C
        ∇u = function_gradient(cv, qp, ue)
        F = one(∇u) + ∇u
        C = tdot(F) # F' ⋅ F
        # Compute stress and tangent
        S, ∂S∂C = constitutive_driver(C, mp)
        P = F ⋅ S
        I = one(S)
        ∂P∂F = otimesu(I, S) + 2 * otimesu(F, I) ⊡ ∂S∂C ⊡ otimesu(F', I)

        # Loop over test functions
        for i in 1:ndofs
            # Test function and gradient
            δui = shape_value(cv, qp, i)
            ∇δui = shape_gradient(cv, qp, i)
            # Add contribution to the residual from this test function
            ge[i] += (∇δui ⊡ P - δui ⋅ b) * dΩ

            ∇δui∂P∂F = ∇δui ⊡ ∂P∂F # Hoisted computation
            for j in 1:ndofs
                ∇δuj = shape_gradient(cv, qp, j)
                # Add contribution to the tangent
                ke[i, j] += (∇δui∂P∂F ⊡ ∇δuj) * dΩ
            end
        end
    end

    # Surface integral for the traction
    for face in 1:nfaces(cell)
        if (cellid(cell), face) in ΓN
            reinit!(fv, cell, face)
            for q_point in 1:getnquadpoints(fv)
                t = tn * getnormal(fv, q_point)
                dΓ = getdetJdV(fv, q_point)
                for i in 1:ndofs
                    δui = shape_value(fv, q_point, i)
                    ge[i] -= (δui ⋅ t) * dΓ
                end
            end
        end
    end
end;

function assemble_global!(K, g, dh, cv, fv, mp, u, ΓN)
    n = ndofs_per_cell(dh)
    ke = zeros(n, n)
    ge = zeros(n)

    # start_assemble resets K and g
    assembler = start_assemble(K, g)

    # Loop over all cells in the grid
    @timeit "assemble" for cell in CellIterator(dh)
        global_dofs = celldofs(cell)
        ue = u[global_dofs] # element dofs
        @timeit "element assemble" assemble_element!(ke, ge, cell, cv, fv, mp, ue, ΓN)
        assemble!(assembler, global_dofs, ge, ke)
    end
end;

function solve()
    reset_timer!()

    # Generate a grid
    N = 10
    L = 1.0
    left = zero(Vec{3})
    right = L * ones(Vec{3})
    grid = generate_grid(Tetrahedron, (N, N, N), left, right)

    # Material parameters
    E = 10.0
    ν = 0.3
    μ = E / (2(1 + ν))
    λ = (E * ν) / ((1 + ν) * (1 - 2ν))
    mp = NeoHookeanParameters(μ, λ)

    # Finite element base
    ip = Lagrange{3,RefTetrahedron,1}()
    qr = QuadratureRule{3,RefTetrahedron}(1)
    qr_face = QuadratureRule{2,RefTetrahedron}(1)
    cv = CellVectorValues(qr, ip)
    fv = FaceVectorValues(qr_face, ip)

    # DofHandler
    dh = DofHandler(grid)
    add!(dh, :u, 3) # Add a displacement field
    close!(dh)

    function rotation(X, t)
        θ = pi / 3 # 60°
        x, y, z = X
        return t * Vec{3}((
            0.0,
            L / 2 - y + (y - L / 2) * cos(-θ) - (z - L / 2) * sin(-θ),
            L / 2 - z + (y - L / 2) * sin(-θ) + (z - L / 2) * cos(-θ)
        ))
    end

    dbcs = ConstraintHandler(dh)
    # Add a homogeneous boundary condition on the "clamped" edge
    dbc = Dirichlet(:u, getfaceset(grid, "right"), (x, t) -> [0.0, 0.0, 0.0], [1, 2, 3])
    add!(dbcs, dbc)
    dbc = Dirichlet(:u, getfaceset(grid, "left"), (x, t) -> rotation(x, t), [1, 2, 3])
    add!(dbcs, dbc)
    close!(dbcs)
    t = 0.5
    Ferrite.update!(dbcs, t)

    # Neumann part of the boundary
    ΓN = union(
        getfaceset(grid, "top"),
        getfaceset(grid, "bottom"),
        getfaceset(grid, "front"),
        getfaceset(grid, "back"),
    )

    # Pre-allocation of vectors for the solution and Newton increments
    _ndofs = ndofs(dh)
    un = zeros(_ndofs) # previous solution vector
    u = zeros(_ndofs)
    Δu = zeros(_ndofs)
    ΔΔu = zeros(_ndofs)
    apply!(un, dbcs)

    # Create sparse matrix and residual vector
    K = create_sparsity_pattern(dh)
    g = zeros(_ndofs)

    # Perform Newton iterations
    newton_itr = -1
    NEWTON_TOL = 1e-8
    NEWTON_MAXITER = 30
    prog = ProgressMeter.ProgressThresh(NEWTON_TOL, "Solving:")

    while true
        newton_itr += 1
        # Construct the current guess
        u .= un .+ Δu
        # Compute residual and tangent for current guess
        assemble_global!(K, g, dh, cv, fv, mp, u, ΓN)
        # Apply boundary conditions
        apply_zero!(K, g, dbcs)
        # Compute the residual norm and compare with tolerance
        normg = norm(g)
        ProgressMeter.update!(prog, normg; showvalues=[(:iter, newton_itr)])
        if normg < NEWTON_TOL
            break
        elseif newton_itr > NEWTON_MAXITER
            error("Reached maximum Newton iterations, aborting")
        end

        # Compute increment using conjugate gradients
        @timeit "linear solve" IterativeSolvers.cg!(ΔΔu, K, g; maxiter=1000)

        apply_zero!(ΔΔu, dbcs)
        Δu .-= ΔΔu
    end

    # Save the solution
    @timeit "export" begin
        vtk_grid("hyperelasticity", dh) do vtkfile
            vtk_point_data(vtkfile, dh, u)
        end
    end

    print_timer(title="Analysis with $(getncells(grid)) elements", linechars=:ascii)
    return u
end

u = solve();
