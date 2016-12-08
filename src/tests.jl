using Plots
using LaTeXStrings
using Base.Test

# provides c0, G, m_H, k_B, sigma_sb, a_rad, Msun, Rsun, Lsun
include("constants.jl")
include("physics.jl")
include("opacity.jl")
include("derivs.jl")
include("nr.jl")
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

function test_profile(profile; mass=5)
    X = 0.7
    Y = 0.28
    m = 0.1 * Msun
    M = mass * Msun
    mf = 0.8 * M
    Rs, Ls, Pc, Tc = init_guess(M, X, Y)
    Rs, Ls, Ps, Ts = load2(M, Rs, Ls, X, Y)
    println("surface: ", [Rs, Ls, Ps, Ts])
    Rc, Lc, Pc, Tc = load1(m, Pc, Tc, X, Y)
    println("center: ", [Rc, Lc, Pc, Tc])

    rhos = rho_ideal(Ps, Ts, mu)
    rhoc = rho_ideal(Pc, Tc, mu)
    kappac = opacity(logkappa_spl, rhoc, Tc)
    kappas = opacity(logkappa_spl, rhos, Ts)
    println("kappa at center: ", kappac)
    println("log R at center: ", log10(rhoc / (Tc * 1e-6) ^ 3))
    println("kappa at surface: ", kappas)
    println("log R at surface: ", log10(rhos / (Ts * 1e-6) ^ 3))

    x1, y1, x2, y2 = profiles(m, M, Rs, Ls, Pc, Tc, X, Y, mf)

    idx_dict = Dict("r" => 1, "l" => 2, "P" => 3, "T" => 4)
    idx = idx_dict[profile]
    if idx == 1
        norm = Rsun
    elseif idx == 2
        norm = Lsun
    else
        norm = 1
    end
    if all(y1[1:end, idx] .== y1[1, idx])
        println("outward flatlined")
    end
    if all(y2[1:end, idx] .== y2[1, idx])
        println("inward flatlined")
    end

    println(profile, "_out(mf) = ", y1[end, idx]) 
    println(profile, "_in(mf) = ", y2[end, idx])   
    plot(x1 / Msun, y1[1:end, idx] / norm, label=profile * "out", yscale=:log10)
    plot!(x2 / Msun, y2[1:end, idx] / norm, label=profile * "in", yscale=:log10)
    # return x1, x2, y1[1:end, idx], y2[1:end, idx]
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

function test_opacity(M=5 * Msun)
    X = 0.7
    Y = 0.28
    mu = mu_from_composition(X, Y)
    
    m = 0.01 * Msun
    # M = 5 * Msun
    mf = 0.8 * M
    Rs, Ls, Pc, Tc = init_guess(M, X, Y)
    Rs, Ls, Ps, Ts = load2(M, Rs, Ls, X, Y)
    println("surface: ", [Rs, Ls, Ps, Ts])
    Rc, Lc, Pc, Tc = load1(m, Pc, Tc, X, Y)
    println("center: ", [Rc, Lc, Pc, Tc])
    rhos = rho_ideal(Ps, Ts, mu)
    rhoc = rho_ideal(Pc, Tc, mu)
    temps = logspace(Ts, Tc)
    dens = logspace(rhos, rhoc)
    kappac = opacity(logkappa_spl, rhoc, Tc)
    kappas = opacity(logkappa_spl, rhos, Ts)
    println("kappa at center: ", kappac)
    println("log R at center: ", log10(rhoc / (Tc * 1e-6) ^ 3))
    
    println("kappa at surface: ", kappas)
    println("log R at surface: ", log10(rhos / (Ts * 1e-6) ^ 3))
end
