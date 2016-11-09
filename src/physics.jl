#=
Miscellaneous physics
=#

# provides c0, G, m_H, k_B, sigma_sb, a_rad, Msun, Rsun, Lsun
include("constants.jl")

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
    Prad = a_rad * T .^ 4 / 3
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
    Prad = a_rad * T .^ 4 / 3
    P = Pgas .+ Prad
    rho = P .* m_H .* mu ./ (k_B * T)
    return rho
end # rho_ideal

function acceleration(m, r)
    #=
    Calculates the local acceleration due to gravity.
    =#
    return G * M ./ (R .^ 2)
end # acceleration

function scale_height(m, r, P, T, mu)
    #=
    Calculates the pressure scale height.
    =#
    g = acceleration(m, r)
    rho = rho_ideal(P, T, mu)
    return P ./ (rho .* g)
end # scale_height
    
function T_surface(R, L_surface)
    #=
    Calculates the effective temperature, given the surface luminosity and the 
    stellar radius.
    =#
    return (L_surface ./ (4 * pi * R .^ 2 * sigma_sb)) .^ (1 / 4)
end # T_surface
    
function P_surface(M, R, kappa)
    #=
    Calculates the surface pressure, given the total mass, radius, and surface
    opacity.
    =#
    g = acceleration(M, R)
    return 2 * g ./ (3 * kappa)
end # P_surface

