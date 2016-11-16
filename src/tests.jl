using Plots
using LaTeXStrings
using Base.Test

# provides c0, G, m_H, k_B, sigma_sb, a_rad, Msun, Rsun, Lsun
include("constants.jl")
include("physics.jl")
include("opacity.jl")
include("nuclear.jl")
include("convection.jl")
include("derivs.jl")
include("shootf.jl")

logkappa_spl = opal_test_spl()

# solar values
mc = 0.001 * Msun
Pc = 2.4e17 # dyne / cm2
Tc = 1.6e7 # K
X = 0.73
Y = 0.24
Z = 1 - X - Y
mu = mu_from_composition(X, Y)
    
function test_nuclear(savefig=nothing)
    temps = logspace(6.5, 8)
    X = 0.7
    Y = 0.28
    Z = 0.02
    mu = mu_from_composition(X, Y)
    rho = rho_ideal(Pc, temps, mu)
    pp = epsilon_pp(rho, temps, X)
    cno = epsilon_cno(rho, temps, X, Y)
    
    print("pp/CNO at T6=18, rho=80: ", epsilon_pp(80, 18e6, X) / epsilon_cno(80, 18e6, X, Y))
    
    plot(temps, pp, label="pp",
         xlabel="Temperature [K]", xscale=:log,
         ylabel=L"Luminosity density [g cm$^{-3}$]", yscale=:log)
    plot!(temps, cno, label="CNO")
    plot!(temps, pp + cno, label="pp + CNO")
    if ~is(savefig, nothing)
        png(savefig)
    end # if
end

function test_dTdm()
    mu = 0.6
    rhoc = rho_ideal(Pc, Tc, mu)
    kappac = 10 ^ logkappa_spl(log10(rhoc), log10(Tc))
    rc, lc, _, __ = load1(Pc, Tc, X, Y)
    mc = 0.001 * Msun
    dTdm(mc, rc, lc, Pc, Tc, kappac)
end

function test_load()
    # load1
    center = load1(mc, Pc, Tc, X, Y)
    println("center ", center)
    # load2
    surface = load2(Msun, Rsun, Lsun, X, Y)
    println("surface ", surface)
    Rs, Ls, Ps, Ts = surface
    rhos = rho_ideal(Ps, Ts, mu)
    println(rhos)
    kappa = 10 ^ logkappa_spl(log10(rhos), log10(Ts))
    # test that the opacity gives a pressure consistant to within 0.001 percent
    @test_approx_eq_eps Ps P_surface(Msun, Rsun, kappa) (Ps * 1e-5)
end    

function test_deriv()
    rc, lc, _, _ = load1(mc, Pc, Tc, X, Y)
    println("deriv at center: ", deriv(mc, rc, lc, Pc, Tc, X, Y))
    rs, ls, Ps, Ts = load2(Msun, Rsun, Lsun, X, Y)
    println("deriv at surface: ", deriv(Msun, rs, ls, Ps, Ts, X, Y))    
end
