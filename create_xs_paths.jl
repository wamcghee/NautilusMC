function read_ace_header(filename)
    open(filename) do f
        idx = 1
        line = readline(f)
        if length(split(line)) == 3
            readline(f)
            line = readline(f)
        end
        HZ, AW, TZ, HD = split(line, " ", keepempty = false)
        suffix = HZ[end-2 : end]
        if suffix == "00c"
            T = 293.6
        elseif suffix == "01c"
            T = 600.0
        elseif suffix == "02c"
            T = 900.0
        elseif suffix == "03c"
            T = 1200.0
        elseif suffix == "04c"
            T = 2500.0
        elseif suffix == "05c"
            T = 0.1
        elseif suffix == "06c"
            T = 250.0
        end
        return HZ, T, parse(Float64, AW)
    end
end

xs_dir = "C:\\Users\\warner\\Desktop\\My Stuff\\Personal Learning\\Monte Carlo\\Lib80x\\Lib80x"
#=
xs_paths = []
for (root, dirs, files) in walkdir(xs_dir)
    for filename in joinpath.(root, files)
        HZ, T, AW = read_ace_header(filename)
        push!(xs_paths, [HZ, AW, T, filename])
        #println(HZ)
    end
end
=#
ZAs = Int64[]
Ts = Float64[]
xs_paths = String[]
for (root, dirs, files) in walkdir(xs_dir)
    for filename in joinpath.(root, files)
        HZ, T, AW = read_ace_header(filename)
        ZA = parse(Int, HZ[1:end-4])
        push!(ZAs, ZA)
        push!(Ts, T)
        push!(xs_paths, filename)
        #println(HZ)
    end
end

order = sortperm(ZAs)
xs_paths = [ZAs[order] Ts[order] xs_paths[order]]

xs_dir_filename = "C:\\Users\\warner\\Desktop\\My Stuff\\Personal Learning\\Monte Carlo\\Lib80x\\xsdir.txt"
using DelimitedFiles
writedlm(xs_dir_filename, xs_paths, ',')
