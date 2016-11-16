#=
Nuclear burning.
=#

# provides c0, G, m_H, k_B, sigma_sb, a_rad
include("constants.jl")

function epsilon_pp(rho, T, X)
    #=
    Proton-proton chain luminosity density.
    Kippenhahn chapter 18 page 194, from Agnulo+1999
    =#
    # approximate the pp2, pp3 correction as a logistic curve from 1 to 1.5
    # compare with Fig 18.7 in Kippenhahn for small 4He fractions
    f(x) = 1 + 0.5 ./ (1 + e .^ (-3 .* (x - 2)))    
    phi = f(T / 1e7)     # pp2, pp3 correction
    f11 = 1       # weak screening approximation
    T9 = T / 1e9
    g11 = 1 + 3.82 * T9 + 1.51 * T9 .^ 2 + 0.144 * T9 .^ 3 - 0.0114 * T9 .^ 4
    factor = 2.57e4 .* phi .* f11 .* g11 .* rho .* X ^ 2
    return factor .* T9 .^ (-2 / 3) .* e .^ (-3.381 ./ (T9 .^ (1/3)))
end # epsilon_pp

function epsilon_cno(rho, T, X, Y)
    #=
    CNO luminosity density
    Kippanhahn chapter 18 page 196
    =#
    T9 = T / 1e9
    X_CNO = (1 - X - Y) * 0.7
    g141 = 1 - 2 .* T9 + 3.41 .* T9 .^ 2 - 2.43 .* T9 .^ 3
    factor = 8.24e25 .* g141 .* X_CNO .*  X .* rho
    return factor .* T9 .^ (-2 / 3) .* e .^ (-15.231 .* T9 .^ (-1 / 3) - (T9 / 0.8) .^ 2)
end # epsilon_cno
