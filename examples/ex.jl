using Tullio
function eq_3(λ⃗, p)
    (; μ, N, κ) = p
    λch  = √(I₁(λ⃗)/3)
    ℒinv = BergstromApproximation
    σJ = κ*log(prod(λ⃗))
    @show σJ
    @tullio σ[i] :=  μ * √(N) * (λ⃗[i]^2 - λch^2)/(λch)*ℒinv(λch/sqrt(N))
    @show σ
    return σ.+σJ
end
