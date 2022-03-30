function MooneyRivlin(p)
    (;C10, C01) = p
    return (λ⃗) -> C10*(I₁(λ⃗)-3)+C01*(I₂(λ⃗) - 3)
end

function NeoHookean(p)
    (; μ) = p
    return (λ⃗) -> μ * (I₁(λ⃗) - 3)
end

function Gent(p)
    (; μ, Jₘ) = p
    return (λ⃗) -> -μ*Jₘ/2*log(1 - (I₁(λ⃗)-3)/Jₘ)
end