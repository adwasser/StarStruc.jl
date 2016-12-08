#=
Root finding and ODE solver.
See Numerical Recipes 16.1 - 16.2
=#

using Calculus
using ODE

function fdjac{T<:AbstractFloat}(f::Function, x::Array{T, 1}; epsilon=1e-6)
    #=
    Forward difference Jacobian
    =#
    n = length(x)
    y = f(x)
    jacobian = Array{Float64}(n, n)
    for i = 1:n
        temp = x[i]
        h = abs(x[i]) * epsilon
        if h == 0.0
            h = epsilon
        end # if
        x[i] = temp + h
        h = x[i] - temp
        ystep = f(x)
        x[i] = temp
        dy = (ystep - y) ./ h
        jacobian[1:end, i] = dy
    end # for
    return jacobian
end # fdjac

function newton(f, x0; epsilon=1e-2, max_count=1000, init_damping=1e-3)
    #=
    Multidimensional Newton-Raphson root finding.
    TODO: figure out scaled stopping critereon
    =#
    x = x0
    count = 1
    @info("newton: Starting Newton method with x0 = ", x0)
    damping = init_damping
    while true
        @info("newton: count = ", count)
        if count > max_count
            @warn("newton: Number of iterations has exceeded ", max_count)
            break
        end # if
        # approximate Jacobian
        # J = fdjac(f, x)
        J = Calculus.finite_difference_jacobian(f, x)
        for i=1:size(J)[1]
            @debug("newton: J", i, " = ", J[i, 1:end])
        end
        # invert for dx
        invJ = inv(J)
        for i=1:size(J)[1]
            @debug("newton: invJ", i, " = ", invJ[i, 1:end])
        end
        dx = -invJ * f(x)
        frac = dx ./ x
        @info("newton: dx / x = ", frac)
        # adaptive damping
        if all(abs(frac) .< 0.1)
            damping = 1.0
        elseif all(abs(frac) .< 1)
            damping = 0.1
        elseif all(abs(frac) .< 10)
            damping = 0.01
        else
            damping = 0.001
        end
        @info("newton: damping = ", damping)
        close_enough = all(abs(frac) .< epsilon)
        x = x + dx .* damping
        @info("newton: x = ", x)        
        if close_enough
            @info("newton: Zero-point found at x = ", x)
            return x
        end # if
        count += 1
    end # while
end # newton

function rkstep{T<:AbstractFloat}(f::Function,
                                  y::Array{T, 1},
                                  x::T,
                                  h::T)
    k1 = f(x, y)
    k2 = f(x + h / 2, y + h / 2 .* k1)
    k3 = f(x + h / 2, y + h / 2 .* k2)
    k4 = f(x + h, y + h .* k3)
    dy = h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return y + dy
end # rkstep

function rk{T<:AbstractFloat}(f::Function,
                              y0::Array{T, 1},
                              xrange::Array{T, 1})
    #=
    Classic fourth-order Runge-Kutta method
    f(x, y) = dy / dx
    Integrate from x0 to xf, y0 = y(x0)
    =#
    # @debug("rk: starting Runge-Kutta, y0 = ", y0, ", h = ", h)
    x0 = xrange[1]
    xf = xrange[end]
    n = length(xrange)
    m = length(y0)
    x = Array{Float64}(n)
    y = Array{Float64}(n, m)
    x[1] = x0
    y[1, 1:end] = y0
    for (i, xstep) in enumerate(xrange[2:end])
        j = i + 1
        h = xstep - x[j - 1]
        yy = rkstep(f, y[j - 1, 1:end], x[j - 1], h)
        x[j] = xstep
        y[j, 1:end] = yy
    end # while
    return x, y
end # rk

function rkck{T<:AbstractFloat}(f::Array{Function, 1},
                                x0::AbstractFloat,
                                y0::Array{T, 1},
                                h0::AbstractFloat,
                                xf::AbstractFloat;
                                epsilon=1e-3)
    #=
    Cash-Karp embedded Runge-Kutta method with adaptive step size.
    f(x, y) = dy / dx
    Integrate from x0 to xf, y0 = y(x0)
    =#
    # Butcher tableau for Cash-Karp
    a = [0, 1 / 5, 3 / 10, 3 / 5, 1, 7 / 8]
    b = zeros(6, 5)
    b[2, 1] = 1/5
    b[3, 1:2] = [3 / 40, 9 / 40]
    b[4, 1:3] = [3 / 10, -9 / 10, 6 / 5]
    b[5, 1:4] = [-11 / 54, 5 / 2, -70 / 27, 35 / 27]
    b[6, :] = [1631 / 55296, 175 / 512, 575 / 13824, 44275 / 110592, 253 / 4096]
    c = [37 / 378, 0, 250 / 621, 125 / 594, 0, 512 / 1771]
    c_embedded = [2825 / 27648, 0, 18575 / 48384, 13525 / 55296, 277 / 14336, 1 / 4]

    # number of equations
    m = length(y0)
    
    x_list = []
    y_list = []
    h_list = []
    push!(x_list, x0)
    push!(y_list, y0)
    push!(h_list, h0)
    while x_list[end] < xf
        xx = x_list[end]
        yy = y_list[end]
        hh = h_list[end]
        k = zeros(6, m) # todo, speed up by putting column loop inside
        k[1, :] = map(g -> hh * g(xx, yy), f)
        for i = 2:6
            for j = 1:m
                k[i, j] = hh * f[j](xx + a[i] * hh, yy + dot(vec(b[i, 1:i - 1]), k[1:i - 1, j]))
            end # for
        end # for
        delta = abs(map(i -> dot(c - c_embedded, k[:, i]), 1:m))
        delta_target = epsilon * (abs(yy) + abs(vec(k[1, :])))
        if any(delta .> delta_target)
            # reject step, try a smaller step size
            h_list[end] = 0.9 * hh * minimum(delta_target ./ delta) ^ 0.25
            continue
        else
            # accept step
            push!(h_list, 0.9 * hh * minimum(delta_target ./ delta) ^ 0.2)
        end # if
        # println("accepted")
        x_next = xx + hh
        y_next = map(i -> yy[i] + dot(c, k[:, i]), 1:m)         
        push!(x_list, x_next)
        push!(y_list, y_next)
    end # while
    x = Array{Float64}(x_list)
    h = Array{Float64}(h_list)
    y = Array{Float64}(length(y_list), length(y_list[1]))
    for i in 1:length(y_list)
        y[i, 1:end] = y_list[i]
    end #for
    return x, y, h
end # rkck

