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
    return G * m ./ (r .^ 2)
    
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

function del_ad()
    return 0.4
end # nabla_ad

function del_rad(m, l, P, T, kappa)
    #=
    Gradient (d logP / d logT) in fully radiative case.
    =#
    return 3 * kappa .* l .* P ./ (16 * pi * a_rad * c0 * G * m .* T .^ 4)
end # nabla_rad

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
