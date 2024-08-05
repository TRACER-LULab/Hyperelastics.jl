using Gridap
using Hyperelastics
using LineSearches: BackTracking
E = 10.0
ν = 0.3
const μ = E / 2 / (1 + ν)
const λ = 2 * μ * ν / (1 - 2 * ν)
const ps = (μ = μ, λ = λ)

F(∇u) = one(∇u) + ∇u
J(F) = det(F)
B(x) = VectorValue(0.0, -0.5, 0.0)
T(x) = VectorValue(-0.1, 0.0, 0.0)

function ψ(F)
    C = F ⊙ F
    I = [tr(C), 0.5 * (tr(C)^2 - tr(C^2)), det(C)]
    j = J(F)
    W = StrainEnergyDensity(NeoHookean(InvariantForm()), I, ps) - ps.μ * log(j) + ps.λ / 2 * log(j)^2
    return W
end

S(F) = ∇(ψ)(F)

domain = (0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
partition = (10, 10, 10)
model = CartesianDiscreteModel(domain, partition)

labels = get_face_labeling(model)
add_tag_from_tags!(labels, "rightFace", [2, 4, 6, 8, 14, 16, 18, 20, 26])
add_tag_from_tags!(labels, "leftFace", [1, 3, 5, 7, 13, 15, 17, 19, 25])

reffe = ReferenceFE(lagrangian, VectorValue{3,Float64}, 1)
V = TestFESpace(model, reffe, conformity=:H1, dirichlet_tags=["leftFace", "rightFace"])

Ω = Triangulation(model)
Γ = BoundaryTriangulation(model)

degree = 2
dΩ = Measure(Ω, degree)
dΓ = Measure(Γ, degree)

energy(u) = ∫(ψ ∘ F ∘ ∇(u)-(B∘u)⋅u)*dΩ - ∫((T∘u)⋅u)*dΓ # You also will need to add the external loads here
res(u, v) = gradient(energy, u)
jac(u, du, v) = hessian(energy, u) # I have realized that you have to explicitly build the jacobian like this


nls = NLSolver(
    show_trace=true,
    method=:newton,
    linesearch=BackTracking()
)

g0(x) = VectorValue(0.0, 0.0, 0.0)
g1(x) = VectorValue(
    0.0,
    (0.5 + (x[2] - 0.5) * cos(π / 3) - (x[3] - 0.5) * sin(π / 3) - x[2]) / 2,
    (0.5 + (x[2] - 0.5) * sin(π / 3) + (x[3] - 0.5) * cos(π / 3) - x[3]) / 2
)

U = TrialFESpace(V, [g0, g1])

#FE problem
op = FEOperator(res, jac, U, V)
solver = FESolver(nls)

x0 = zeros(Float64, num_free_dofs(V))
uh = FEFunction(U, x0)
uh, = solve!(uh, solver, op)

writevtk(Ω, "results_$(lpad(step,3,'0'))", cellfields=["uh" => uh])
