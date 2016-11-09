#=
Nuclear burning.
=#

# provides c0, G, m_H, k_B, sigma_sb, a_rad
include("constants.jl")

function epsilon_pp(P, T, mu, X)
    #=
    Proton-proton chain luminosity density.
    Kippenhahn chapter 18 page 194, from Agnulo+1999
    =#
    phi = 1.5     # pp2, pp3 correction
    f11 = 1       # weak screening approximation
    T9 = T / 1e9
    g11 = 1 + 3.82 * T9 + 1.51 * T9 .^ 2 + 0.144 * T9 .^ 3 - 0.0114 * T9 .^ 4
    rho = rho_ideal(P, T, mu)
    factor = 2.57e4 * phi * f11 * g11 .* rho * X ^ 2
    return factor .* T9 .^ (-2 / 3) .* e .^ (-3.381 ./ (T9 .^ (1/3)))
end # epsilon_pp
