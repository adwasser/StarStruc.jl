#=
ODE solver.
See Numerical Recipes 16.1 - 16.2
=#

function rk{T<:AbstractFloat}(f::Array{Function, 1},
                              x0::AbstractFloat,
                              y0::Array{T, 1},
                              h::AbstractFloat,
                              xf::AbstractFloat)
    #=
    Classic fourth-order Runge-Kutta method
    f(x, y) = dy / dx
    Integrate from x0 to xf, y0 = y(x0)
    =#
    n = round(Int, (xf - x0) / h)
    m = length(y0)
    x = Array{Float64}(n)
    y = Array{Float64}(n, m)
    x[1] = x0
    y[1, 1:end] = y0
    for i = 2:n
        xx = x[i - 1]
        yy = transpose(y[i - 1, 1:end])
        k1 = map(g -> g(xx, yy), f)
        k2 = map(g -> g(xx + h / 2, yy + h / 2 * k1), f)
        k3 = map(g -> g(xx + h / 2, yy + h / 2 * k2), f)
        k4 = map(g -> g(xx + h, yy + h * k3), f)
        x[i] = xx + h
        y[i, 1:end] = yy + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
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
