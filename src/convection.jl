#=
Functions for dealing with convective energy transport.
=#

using ValidatedNumerics.RootFinding.find_roots

include("constants.jl")
include("physics.jl")

function nabla(P, T, dP, dT)
    #=
    Gradient (d logP / d logT)
    =#
    return T * dP / (P * dT)
end # nabla

function nabla_ad()
    return 0.4
end # nabla_ad

function nabla_rad(m, l, P, T, kappa)
    #=
    Gradient (d logP / d logT) in fully radiative case.
    =#
    return 3 * kappa .* l .* P ./ (16 * pi * a_rad * c0 * G * m .* T .^ 4)
end # nabla_rad

function nabla_mlt(m, r, l, P, T, mu, kappa)
    #=
    Gradient from mixing length theory.
    ml is the mixing length
    =#
    ml = 1.5 * 5.5e9 # 1.5 times the solar ML at R/2
    c_P = 5 / 2 * k_B / (m_H * mu)
    rho = rho_ideal(P, T, mu)
    g = acceleration(m, r)
    Hp = scale_height(m, r, P, T, mu)
    delta = 1 # for an ideal gas
    grad_rad = nabla_rad(m, l, P, T, kappa)
    U = 3 * a_rad * c0 * T .^ 3 / (c_P * rho .^2 .* kappa * ml ^ 2) .*
        sqrt(8 * Hp ./ (g * delta))
    W = grad_rad - nabla_ad()
    
    f(xi) = (xi - U) ^ 3 + 8 * U / 9 * (xi ^ 2 - U ^ 2 - W)
    res = find_roots(f, 0, grad_rad)
    # println(res)
    return res[1].interval.lo
end # nabla_mlt
