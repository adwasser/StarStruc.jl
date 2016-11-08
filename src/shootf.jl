#=
Implementation of shootf, the Numerical Recipes
algorithm for solving a boundary value problem
by shooting to a fixed point.
=#

# provides c0, G, m_H, k_B, sigma_sb, a_rad
include("constants.jl")
include("physics.jl")
include("opacity.jl")

logkappa_spl = opal_test_spl()

function load1(Pc, Tc, X, Y)
    #=
    Subroutine for shootf.
    Loads the (r, l, P, T) center (m = 0) values for starting the outward 
    integration.
    =#
    mu = mu_from_composition(X, Y)
    m = 0.01 # needs tuning
    rhoc = rho_ideal(Pc, Tc, mu)
    rc = (3 * m / (4 * pi * rhoc)) ^ (1 / 3) 
    lc = epsilon_pp(Pc, Tc, mu, X) * m
    return [rc, lc, Pc, Tc]
end # load1

function load2(M, Rs, Ls, X, Y)
    #=
    Subroutine for shootf.
    Loads the (r, l, P, T) surface (m = M) values for starting the inward
    integration.
    =#
    mu = mu_from_composition(X, Y)
    converged = false
    rho_guess = 1e-6 # need a less arbitrary guess
    while ~converged
        Ts = T_surface(Rs, Ls)
        kappa = 10 ^ logkappa_spl(log10(rho_guess), log10(Ts))
        Ps = P_surface(M, R, kappa)
        rho = rho_ideal(Ps, Ts, mu)
        relerr = (rho - rho_guess) / rho
        if abs(relerr) < 0.1 # need less arbitrary convergence condition
            converged = true
        else
            rho_guess = rho_guess + relerr * rho / 2
        end # if
    end # while
    return [Rs, ls, Ps, Ts]
end # load2

function deriv(m, r, l, P, T, X, Y)
    #=
    Subroutine for shootf.
    Calulates the derivatives of (r, l, P, T) with respect to mass.
    See Kippenhahn chapter 10, page 89.
    =#
    mu = mu_from_composition(X, Y)
    rho = rho_ideal(P, T, mu)
    kappa = 10 ^ (logkappa_spl(log10(rho), log10(T)))
    return [drdm(m, r, l, P, T, mu),
            dldm(m, r, l, P, T, X, mu),
            dPdm(m, r, l, P, T),
            dTdm(m, r, l, P, T, kappa)]
end # deriv

function shootf(M, mu, )
    #= 
    Shooting to a fixed point.
    Go from initial guesses to a solution to the two-point boundary value
    problem.
    =#
end # shootf
