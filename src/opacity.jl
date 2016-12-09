#=
Reading in and interpolating over opacity tables.
=#

using Dierckx

# regex for matching whitespace
whitespace = r" +"

function split_tables(filename)
    #=
    Split the convoluted opacity tables into array of table strings.

    Parameters
    ----------
    filename : string, file to load

    Returns
    -------
    tables : array of table strings
    =#
    return 0
end

function load_table(table_string)
    #=
    Load the convoluted opacity tables into a sample of opacity at (density, temperature) samples.

    Parameters
    ----------
    table_string : string with table data

    Returns
    -------
    spl : Spline2D
    =#

    lines = split(table_string, '\n')[1: end - 1]
    
    # get mass fractions
    X, Y, Z = [float(split(s, '=')[2]) for s in split(lines[1], ' ')[end - 4: end - 2]]

    # get logR coords
    logR = float(split(lines[5], whitespace)[2: end - 1])

    # get logT coords and kappa
    log_rho = []
    log_T = []
    log_kappa = []
    log_R = []
    for (i, line) in enumerate(lines[7: end])
        entries = float(split(line, whitespace))
        logT = entries[1]
        temp = 10 ^ logT
        for (j, entry) in enumerate(entries[2: end])
            if entry == 9.999
                continue
            end
            push!(log_R, logR[j])
            dens = 10 ^ logR[j] * (temp * 1e-6) ^ 3
            push!(log_T, logT)
            push!(log_rho, log10(dens))
            push!(log_kappa, entry)
        end
    end
    log_rho = map(Float64, log_rho)
    log_T = map(Float64, log_T)
    log_R = map(Float64, log_R)
    log_kappa = map(Float64, log_kappa)
    return Spline2D(log_R, log_T, log_kappa, s=1)
end

function opal_test_spl()
    table_string = readstring(open(joinpath(Pkg.dir("StarStruc"), "data", "test.txt")))
    spl = load_table(table_string)
    return spl
end

function opacity(spl, rho, T)
    R = rho ./ (T / 1e6) .^ 3
    return 10 .^ spl(log10(R), log10(T))
end
