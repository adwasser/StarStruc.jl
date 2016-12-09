#=
Implementation of shootf, the Numerical Recipes
algorithm for solving a boundary value problem
by shooting to a fixed point.
=#

# provides c0, G, m_H, k_B, sigma_sb, a_rad, Msun, Rsun, Lsun
include("constants.jl")
include("physics.jl")
include("opacity.jl")
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
    rc = (3 * m / (4 * pi * rhoc)) ^ (1 / 3)
    lc = dldm(m, 0, 0, 0, Tc, rhoc, X, Y) * m
    # @debug("load1: center values = ", [rc, lc, Pc, Tc])
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
    Ts = try
        T_surface(Rs, Ls)
    catch
        @warn("load2: Rs = ", Rs)
        @warn("load2: Ls = ", Ls)        
    end
    Ps = 0
    count = 1
    while ~converged & (count < 100)
        kappa = opacity(logkappa_spl, rho_guess, Ts)
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
    # @debug("load2: surface values = ", [Rs, Ls, Ps, Ts])
    return [Rs, Ls, Ps, Ts]
end # load2

function init_guess(M, X, Y)
    #=
    Gives a reasonable guess for Pc, Tc, R, and L for a star of mass M
    =#
    yhalf = [Rsun * 10 ^ 0.062544715341505616, Lsun * 10 ^ -0.36630780815983022, 10 ^ 19.700356657315954, 10 ^ 7.2863105891013289]
    y2 = [Rsun * 10 ^ 0.36072469919806871, Lsun * 10 ^ 1.3021335874028637, 10 ^ 17.186879800538229, 10 ^ 7.3557279873577759]
    y5 = [Rsun * 10 ^ 0.59322648100945885, Lsun * 10 ^ 2.8743140385067893, 10 ^ 16.752833342638592, 10 ^ 7.4547058342916337]
    y4 = [Rsun * 10 ^ 0.52962615387538248, Lsun * 10 ^ 2.50275963393719, 10 ^ 16.855407680647915, 10 ^ 7.429616124439673]
    y10 = [Rsun * 10 ^ 0.68191099104456998, Lsun * 10 ^ 3.8447263619458569, 10 ^ 16.517451455410988, 10 ^ 7.5012272142698366] 
    if M == Msun
        return [Rsun, Lsun, 2.4e17, 1.6e7]
    elseif M == 0.5 * Msun
        @info("init_guess: from MESA-Web")
        return yhalf
    elseif M == 2 * Msun
        @info("init_guess: from MESA-Web")
        return y2
    elseif (M >= 4 * Msun) & (M <= 5 * Msun)
        @info("init_guess: from MESA-Web")
        dm = M - 4 * Msun
        dydm = (y5 - y4) / Msun
        return y4 + dydm * dm
    elseif M == 10 * Msun
        @info("init_guess: from MESA-Web")
        return y10
    end
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

function profiles(m, M, Rs, Ls, Pc, Tc, X, Y, mf; n=1000)
    #=
    Integrates outward from m to mf and inward from M to mf to construct
    r, l, P, and T profiles from the given initial conditions.
    =#
    D(x, y) = deriv(x, y[1], y[2], y[3], y[4], X, Y)
    # load the center BC guess
    yc = load1(m, Pc, Tc, X, Y)
    # load the surface BC guess
    ys = load2(M, Rs, Ls, X, Y)
    # create the mass grid
    mass_out = collect(linspace(m, mf, n))
    mass_in = collect(linspace(M, mf, 100 * n))
    # integrate out from center to fixed point
    x1, y1 = rk(D, yc, mass_out)
    # integrate in from surface to fixed point
    x2, y2 = rk(D, ys, mass_in)
    return x1, y1, x2, y2
end

function profiles(star::Star, Rs, Ls, Pc, Tc; n=1000)
    return profiles(star.m, star.M, Rs, Ls, Pc, Tc, star.X, star.Y, star.mf)
end

function score(m, M, Rs, Ls, Pc, Tc, X, Y, mf; n=Int(1e3))
    #=
    Function to zero.  Integrates out and in.
    =#
    x1, y1, x2, y2 = profiles(m, M, Rs, Ls, Pc, Tc, X, Y, mf)
    yf1 = y1[end, 1:end]
    yf2 = y2[end, 1:end]
    return (yf1 - yf2) ./ yf1
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

function shootf(star::Star)
    return shootf(star.m, star.M, star.X, star.Y)
end

function shootf!(star::Star)
    bc = shootf(star)
    x1, y1, x2, y2 = profiles(star, bc...)
    star.grid = vcat(x1, reverse(x2))
    star.r = vcat(y1[1:end, 1], reverse(y2[1:end, 1]))
    star.l = vcat(y1[1:end, 2], reverse(y2[1:end, 2]))
    star.P = vcat(y1[1:end, 3], reverse(y2[1:end, 3]))
    star.T = vcat(y1[1:end, 4], reverse(y2[1:end, 4]))
    return x1, y1, x2, y2
end
