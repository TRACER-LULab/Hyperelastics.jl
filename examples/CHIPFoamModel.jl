using Hyperelastics
using DifferentiationInterface
import ForwardDiff
# ---
function J̄(J, (;pg, C10, φ₀, K))
    p̃ = pg + C10 * (φ₀^(1 / 3) * (4J - 4 + 5φ₀) / (J - 1 + φ₀) - (4J - 1) * (4J + 1) / (3J^(4 / 3)))
    Jm = exp(-p̃ / K)
    return J / Jm
end

function dJ̄dJ1((;pg, C10, φ₀, K))
    cbrt_φ₀ = cbrt(φ₀)
    cbrt_φ₀2 = cbrt_φ₀^(2)
    (-C10 - 4C10 * (cbrt_φ₀2) + K * (cbrt_φ₀2)) / (K * exp((5.0C10 - pg - 5.0C10 * (cbrt_φ₀)) / K) * (cbrt_φ₀2))
end
#
using Symbolics
@variables J pg C10 φ₀ K
p̃ = pg + C10 * (φ₀^(1 // 3) * (4J - 4 + 5φ₀) / (J - 1 + φ₀) - (4J - 1) * (4J + 1) / (3J^(4 // 3)))
Jm = exp(-p̃ / K)
J̄ = J / Jm
dJ̄dJ = Symbolics.derivative(J̄, J)
dJ̄dJ_1 = substitute(dJ̄dJ, Dict(J => 1, ))
