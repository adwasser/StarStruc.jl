#=
Derivatives of (r, l, P, T) with respect to the mass coordinate.
=#

include("constants.jl")
include("physics.jl")
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
    grad_rad = del_rad(m, l, P, T, kappa)
    grad_ad = del_ad()
    stable =  grad_rad < grad_ad
    if stable
        # @debug("dTdm: stable to convection")
        grad = grad_rad
    else
        # @debug("dTdm: unstable to convection")
        grad = grad_ad
    end # if
    factor = - G .* m .* T ./ (4 * pi .* r .^ 4 .* P)
     return factor .* grad
end # dTdm

function deriv(m, r, l, P, T, X, Y)
    #=
    Subroutine for shootf.
    Calulates the derivatives of (r, l, P, T) with respect to mass.
    See Kippenhahn chapter 10, page 89.
    =#
    # @debug("deriv: m, y = ", m, ", ", [r, l, P, T])
    mu = mu_from_composition(X, Y)
    rho = rho_ideal(P, T, mu)
    kappa = try
        opacity(logkappa_spl, rho, T)
    catch
        @debug("deriv: Opacity interpolation failed!")
        @debug("deriv: y = ", [r, l, P, T])
        @debug("deriv: logT = ", log10(T))
        @debug("deriv: logR = ", log10(rho / (T / 1e6)^3))
        throw(DomainError)
    end
    dy = [drdm(m, r, l, P, T, rho, mu),
          dldm(m, r, l, P, T, rho, X, Y),
          dPdm(m, r, l, P, T),
          dTdm(m, r, l, P, T, mu, kappa)]
    # @debug("deriv: dy = ", dy)
    return dy
end # deriv
