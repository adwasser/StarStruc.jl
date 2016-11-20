#=
All units are cgs.
=#

module StarStruc

# provides c0, G, m_H, k_B, sigma_sb, a_rad
include("constants.jl")
include("opacity.jl")
include("shootf.jl")
include("odesolve.jl")
include("tests.jl")

export shootf

end # module

