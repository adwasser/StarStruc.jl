#=
All units are cgs.
=#

module StarStruc

# provides c0, m_H, k_B, sigma_sb, a_rad
# include("constants.jl")
include("constants.jl")
include("opacity.jl")

function mu(X, Y)
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
end

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
end
    
function rho(Pgas, T, X, Y)
    #=
    Calculates the density of an ideal fully ionized gas.

    Parameters
    ----------
    Pgas : gas pressure in dyne cm^-2
    T : temperature in K
    X : hydrogen mass fraction
    Y : helium mass fraction

    Returns
    -------
    rho : density in g cm^-3
    =#
    Prad = a_rad * T^4 / 3
    P = Pgas + Prad
    return P * m_H * mu(X, Y) / (k_B * T)
end


end # module

