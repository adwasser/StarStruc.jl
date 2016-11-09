#=
Derivatives of (r, l, P, T) with respect to the mass coordinate.
=#

include("constants.jl")
include("physics.jl")
include("nuclear.jl")
include("convection.jl")

function drdm(m, r, l, P, T, mu)
    #=
    Derivative of r with respect to m, from mass conservation.
    =#
    rho = rho_ideal(P, T, mu)
    return 1 / (4 * pi * r .^ 2 .* rho)
end # drdm

function dPdm(m, r, l, P, T)
    #=
    Derivative of P with respect to m, from hydrostatic equilibrium.
    =#
    return -G * m / (4 * pi * r .^ 4)
end # dPdm

function dldm(m, r, l, P, T, X, mu)
    #=
    Derivative of l with respect to m, from energy conservation.
    =#
    return epsilon_pp(P, T, mu, X)
end # dldm

function dTdm(m, r, l, P, T, kappa)
    #=
    Derivative of T with respect to m, from... thermodynamics.
    =#
    factor = dPdm(m, r, l, P, T) .* P ./ T
    grad_rad = nabla_rad(m, l, P, T, kappa)
    stable =  grad_rad < 0.4
    # println("grad_rad ", grad_rad)
    if stable
        # println("stable")
        grad = grad_rad
    else
        # println("unstable")
        # mixing length theory convection transport
        grad = nabla_mlt(m, r, l, P, T, mu, kappa)
    end # if
    # println("grad ", grad)
    return factor .* grad
end # dTdm
