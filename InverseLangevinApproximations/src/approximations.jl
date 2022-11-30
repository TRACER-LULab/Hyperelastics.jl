CohenRounded3_2(y) = y * (3 - y^2) / (1 - y^2)

ArrudaApproximation(y) = 3y+9/5*y+297/175*y^3+297/175*y^5+1539/875*y^7+126117/67375*y^9

PadeApproximation_3_2(y) = y * (3 - 36 / 35 * y^2) / (1 - 33 / 35 * y^2)

PusoApproximation(y) = 3 * y / (1 - y^3)

TreloarApproximation(y) = 3 * y / (1 - (3 / 5 * y^2 + 36 / 175 * y^4 + 108 / 875 * y^6))

WarnerApproximation(y) = 3 * y / (1 - y^2)

KuhnGrunApproximation(y) = 3y + 9y^3 / 5 + 297y^5 / 175 + 1539y^7 / 875 + 126117y^9 / 67375 + 43733439y^11 / 21896875 + 231321177y^13 / 109484375 + 20495009043y^15 / 9306171875 + 1073585186448381y^17 / 476522530859375 + 4387445039583y^19 / 1944989921875

function BergstromApproximation(y)
    if abs(y) < 0.84136
        return 1.31446tan(1.58986y) + 0.91209y
    elseif 0.84136 <= abs(y) < 1.0
        return 1 / (sign(y) - y)
    end
end

PadeApproximation_1_4(y) = 3y / (1 - 3y^2 / 5 - 36y^4 / 175)

PadeApproximation_1_2(y) = 3y / (1 - 3y^2 / 5)

PadeApproximation_5_0(y) = 3y + 9y^3 / 5 + 297y^5 / 175

PadeApproximation_3_0(y) = 3y + 9y^3 / 5

Jedynak2017(y) = y*(3-773/768*y^2-1300/1351*y^4+501/340*y^6-678/1385*y^8)/(1-y)/(1+866/853*y)
