#=
Implementation of shootf, the Numerical Recipes
algorithm for solving a boundary value problem
by shooting to a fixed point.
=#

# provides c0, G, m_H, k_B, sigma_sb, a_rad, Msun, Rsun, Lsun
include("constants.jl")
include("physics.jl")
include("derivs.jl")
include("nr.jl")

function load1(m, Pc, Tc, X, Y)
    #=
    Subroutine for shootf.
    Loads the (r, l, P, T) center (m = 0) values for starting the outward 
    integration.
    =#
    mu = mu_from_composition(X, Y)
    # m = 0.01 * Msun # needs tuning
    rhoc = rho_ideal(Pc, Tc, mu)
    @debug("load1: rhoc = ", rhoc)
    rc = (3 * m / (4 * pi * rhoc)) ^ (1 / 3)
    lc = dldm(m, 0, 0, 0, Tc, rhoc, X, Y) * m
    @debug("load1: center values = ", [rc, lc, Pc, Tc])
    return [rc, lc, Pc, Tc]
end # load1

function load2(M, Rs, Ls, X, Y; max_iterations=1000)
    #=
    Subroutine for shootf.
    Loads the (r, l, P, T) surface (m = M) values for starting the inward
    integration.
    =#
    mu = mu_from_composition(X, Y)
    converged = false
    rho_guess = 1e-8 # need a less arbitrary guess?
    Ts = T_surface(Rs, Ls)
    Ps = 0
    count = 1
    while ~converged & (count < 1000)
        kappa = 10 ^ logkappa_spl(log10(rho_guess), log10(Ts))
        Ps = P_surface(M, Rs, kappa)
        rho = rho_ideal(Ps, Ts, mu)
        relerr = (rho - rho_guess) / rho
        if abs(relerr) < 1e-6 # need less arbitrary convergence condition
            converged = true
        else
            rho_guess = rho_guess + relerr * rho / 6
            count += 1
        end # if
    end # while
    if count >= max_iterations
        @warn("load2: Exceeded ", max_iterations, " iterations!")
    end # if
    @debug("load2: surface values = ", [Rs, Ls, Ps, Ts])
    return [Rs, Ls, Ps, Ts]
end # load2

function init_guess(M, X, Y)
    #=
    Gives a reasonable guess for Pc, Tc, R, and L for a star of mass M
    =#
    mu = mu_from_composition(X, Y)

    # mass-radius relation
    if M / Msun > 1.3 # then CNO burning dominates
        z1 = 0.8
        z2 = 0.64
    else # then pp chain dominates
        z1 = 0.5
        z2 = 0.1
    end # if
    R = Rsun * (M / Msun) ^ z1 * (mu / 0.61) ^ z2
    # mass-luminosity relation
    L = Lsun * (M / Msun) ^ 3 * (mu / 0.61) ^ 4
    # Pc, Tc guesses, essentially dimensional analysis with fudges
    Pc = 2 * G * M ^ 2 / (pi * R ^ 4)
    Tc = 8 * 0.02 * mu * m_H * G * M / (k_B * R)
    return [R, L, Pc, Tc]
end # init_guess

function score(m, M, Rs, Ls, Pc, Tc, X, Y, mf)
    #=
    Function to zero.  Integrates out and in.
    =#
    D(x, y) = deriv(x, y[1], y[2], y[3], y[4], X, Y)
    # load the center BC guess
    yc = load1(m, Pc, Tc, X, Y)
    # load the surface BC guess
    ys = load2(M, Rs, Ls, X, Y)
    # integrate out from center to fixed point
    @debug("score: Integrating from m = ", m, " to mf = ", mf)
    x1, y1 = ode45(D, yc, [m, mf])
    # x1, y1 = rk(D, yc, [m, mf], 1e-3 * M)
    # integrate in from surface to fixed point
    @debug("score: Integrating from M = ", M, " to mf = ", mf)
    x2, y2 = ode45(D, ys, [M, mf])
    # x2, y2 = rk(D, ys, [M, mf], -1e-3 * M)
    return y1[end] - y2[end]
end # score

function shootf(m, M, X, Y; fixed_point=0.8)
    #= 
    Shooting to a fixed point
    Go from initial guesses to a solution to the two-point boundary value
    problem.
    =#
    y0 = init_guess(M, X, Y)
    @info("shootf: Initial guess for M = ", M, ":")
    @info("shootf:     y0 = ", y0)
    mf = fixed_point * M
    f(y) = score(m, M, y[1], y[2], y[3], y[4], X, Y, mf)
    Rs, Ls, Pc, Tc = newton(f, y0)
    return [Rs, Ls, Pc, Tc]
end # shootf
