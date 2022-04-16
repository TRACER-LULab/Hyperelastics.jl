using Hyperelastics
using GalacticOptim
using Optim
using Gridap
using LineSearches: BackTracking

W = NeoHookean((μ=15e3,))

# Deformation
F(∇u) = one(∇u) + ∇u' # 1+δu

# Right Cauchy Deformation Tensor
C(F) = (F') * F

# Cauchy Stress Tensor
σ(∇u) = (1.0 / J(C(F(∇u)))) * F(∇u) ⋅ S(∇u) ⋅ (F(∇u))'

# Invariants
I₁(C) = tr(C)
I₂(C) = 1 / 2 * (tr(C)^2 - tr(C^2))
I₃(C) = det(C)
J(C) = sqrt(I₃(C))

# Consituitive Laws
function S(∇u)
    Cinv = inv(C(F(∇u)))
    1.0 * (one(∇u) - Cinv) + 100.0 * log(J(F(∇u))) * Cinv
end

function dS(∇du, ∇u)
    Cinv = inv(C(F(∇u)))
    _dE = dE(∇du, ∇u)
    100.0 * (Cinv ⊙ _dE) * Cinv + 2 * (1.0 - 100.0 * log(J(F(∇u)))) * Cinv ⋅ _dE ⋅ (Cinv')
end

# Model
domain = (0, 1, 0, 1)
partition = (20, 20)
model = CartesianDiscreteModel(domain, partition)

# Define Boundaries
labels = get_face_labeling(model)
add_tag_from_tags!(labels, "diri_0", [1,3,7])
add_tag_from_tags!(labels, "diri_1", [2,4,8])

# Setup Integration
degree = 2
Ω = Triangulation(model)
dΩ = Measure(Ω, degree)

# Weak Form
res(u, v) = ∫( (dE∘))