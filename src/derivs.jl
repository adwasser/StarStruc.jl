#=
Derivatives of (r, l, P, T) with respect to the mass coordinate.
=#

include("constants.jl")
include("physics.jl")
include("nuclear.jl")
include("convection.jl")
include("opacity.jl")

# TODO: implement an opacity that changes with composition
const logkappa_spl = opal_test_spl()

function drdm(m, r, l, P, T, rho, mu)
    #=
    Derivative of r with respect to m, from mass conservation.
    =#
    return 1 / (4 * pi * r .^ 2 .* rho)
end # drdm

function dldm(m, r, l, P, T, rho, X, Y)
    #=
    Derivative of l with respect to m, from energy conservation.
    =#
    return epsilon_pp(rho, T, X) + epsilon_cno(rho, T, X, Y)
end # dldm

function dPdm(m, r, l, P, T)
    #=
    Derivative of P with respect to m, from hydrostatic equilibrium.
    =#
    return -G * m / (4 * pi * r .^ 4)
end # dPdm

function dTdm(m, r, l, P, T, mu, kappa)
    #=
    Derivative of T with respect to m, from... thermodynamics.
    =#
    factor = dPdm(m, r, l, P, T) .* P ./ T
    grad_rad = nabla_rad(m, l, P, T, kappa)
    stable =  grad_rad < 0.4
     if stable
        @debug("dTdm: stable to convection")
        grad = grad_rad
    else
        @debug("dTdm: unstable to convection")
        # mixing length theory convection transport
        grad = nabla_mlt(m, r, l, P, T, mu, kappa)
    end # if
     return factor .* grad
end # dTdm

function deriv(m, r, l, P, T, X, Y)
    #=
    Subroutine for shootf.
    Calulates the derivatives of (r, l, P, T) with respect to mass.
    See Kippenhahn chapter 10, page 89.
    =#
    @debug("deriv: y = ", [m, r, l, P, T])
    mu = mu_from_composition(X, Y)
    rho = rho_ideal(P, T, mu)
    kappa = 10 ^ (logkappa_spl(log10(rho), log10(T)))
    return [drdm(m, r, l, P, T, rho, mu),
            dldm(m, r, l, P, T, rho, X, Y),
            dPdm(m, r, l, P, T),
            dTdm(m, r, l, P, T, mu, kappa)]
end # deriv
