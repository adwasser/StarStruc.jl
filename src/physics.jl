#=
Miscellaneous physics, may be reorganized soon.
=#

# provides c0, G, m_H, k_B, sigma_sb, a_rad
include("constants.jl")
include("nuclear.jl")

function mu_from_composition(X, Y)
    #=
    Calculates the dimensionless mean molecular weight for a fully ionized gas.
    
    Parameters
    ----------
    X : hydrogen mass fraction
    Y : helium mass fraction
    
    Returns
    -------
    mu : mean molecular weight
    =#
    return 2 / (1 + 3 * X + 0.5 * Y)
end # mu_from_composition

function beta(Pgas, T)
    #=
    Calculates the ratio between gas pressure and total pressure.

    Parameters
    ----------
    Pgas : gas pressure in dyne cm^-2
    T : temperature in K

    Returns
    -------
    beta : Pgas / P
    =#
    Prad = a_rad * T^4 / 3
    return Pgas / (Pgas + Prad)
end # beta

function rho_ideal(Pgas, T, mu)
    #=
    Calculates the density of an ideal fully ionized gas.

    Parameters
    ----------
    Pgas : gas pressure in dyne cm^-2
    T : temperature in K
    mu : dimensionless mean molecular weight

    Returns
    -------
    rho : density in g cm^-3
    =#
    Prad = a_rad * T^4 / 3
    P = Pgas + Prad
    return P * m_H * mu / (k_B * T)
end # rho_ideal

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
    return 3 * kappa * l * P / (16 * pi * a_rad * c0 * G * m * T ^ 4)
end # nabla_rad

function T_surface(R, L_surface)
    #=
    Calculates the effective temperature, given the surface luminosity and the 
    stellar radius.
    =#
    return (L_surface / (4 * pi * R ^ 2 * sigma_sb)) ^ (1 / 4)
end # T_surface
    
function P_surface(M, R, kappa)
    #=
    Calculates the surface pressure, given the total mass, radius, and surface
    opacity.
    =#
    g = G * M / (R ^ 2)
    return 2 * g / (3 * kappa)
end # P_surface

function drdm(m, r, l, P, T, mu)
    #=
    Derivative of r with respect to m, from mass conservation.
    =#
    rho = rho_ideal(P, T, mu)
    return 1 / (4 * pi * r ^ 2 * rho)
end # drdm

function dPdm(m, r, l, P, T)
    #=
    Derivative of P with respect to m, from hydrostatic equilibrium.
    =#
    return -G * m / (4 * pi * r ^ 4)
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
    factor = dPdm(m, r, l, P, T) * P / T
    grad_rad = lambda_rad(m, l, P, T, kappa)
    stable =  grad_rad < 0.4
    if stable
        grad = grad_rad
    else
        # mixing length convection transport
    end # if
    return factor * grad
end # dTdm
   
