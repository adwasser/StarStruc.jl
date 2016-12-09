#=
All units are cgs.
=#

module StarStruc

using Plots
using LaTeXStrings
using Logging

@Logging.configure(level=INFO)

# provides c0, G, m_H, k_B, sigma_sb, a_rad
include("constants.jl")

type Star
    M::AbstractFloat
    m::AbstractFloat
    mf::AbstractFloat
    X::AbstractFloat
    Y::AbstractFloat
    grid::Array{AbstractFloat}
    r::Array{AbstractFloat}
    l::Array{AbstractFloat}
    P::Array{AbstractFloat}
    T::Array{AbstractFloat}
end

Star(M) = Star(M, 0.01 * Msun, 0.8 * M, 0.7, 0.28,
               zeros(1), zeros(1), zeros(1), zeros(1), zeros(1))
Star(M, X, Y) = Star(M, 0.01 * Msun, 0.8 * M, X, Y,
                     zeros(1), zeros(1), zeros(1), zeros(1), zeros(1))                     

include("opacity.jl")
include("shootf.jl")
include("nr.jl")
include("tests.jl")

function starplot(star::Star)
    p1 = plot(star.grid / Msun, star.r / Rsun, label="r",
              xlabel=L"Mass [$M_\odot$]", ylabel=L"Radius [$R_\odot$]")
    p2 = plot(star.grid / Msun, star.l / Lsun, label="l",
              xlabel=L"Mass [$M_\odot$]", ylabel=L"Luminosity [$L_\odot$]")
    p3 = plot(star.grid / Msun, star.P, label="P", yscale=:log10,
              xlabel=L"Mass [$M_\odot$]", ylabel=L"Pressure [dyne cm$^{-2}$]")
    p4 = plot(star.grid / Msun, star.T, label="T", yscale=:log10,
              xlabel=L"Mass [$M_\odot$]", ylabel="Temperature [K]")
    plot(p1, p2, p3, p4, layout=(2, 2))
end

export shootf
export shootf!
export profiles
export Msun
export Rsun
export Star
export starplot

end # module

