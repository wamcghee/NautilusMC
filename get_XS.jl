
using DelimitedFiles
include("C:\\Users\\warner\\Desktop\\My Stuff\\Personal Learning\\Monte Carlo\\NautilusMC\\NautilusMC\\Materials.jl")
include("C:\\Users\\warner\\Desktop\\My Stuff\\Personal Learning\\Monte Carlo\\NautilusMC\\NautilusMC\\read_ace.jl")


# Set up grid index mapping
E_min = 1.0E-10
logE_min = log10(E_min)
E_max = 1.0E+02
logE_max = log10(E_max)
log_span = logE_max - logE_min

n_grid = 10001 # Number of grid mapping spaces to use

function grid_index(E::Number)
    #=
    Returns index in grid index map for the given energy.

    :param E: Energy of interest
    =#
    f = (log10(E) - logE_min) / log_span # log fraction through the map
    return convert(Int64, floor(f * n_grid)) + 1
end

function grid_energy(idx::Int)
    #=
    Returns energy at grid index of interest.

    :param idx: Index of interest
    =#
    f = (idx - 1) / n_grid # log fraction through the map
    return exp10(logE_min + f*log_span)
end

# Read in table of cross section directories.
xsdir = readdlm("C:\\Users\\warner\\Desktop\\My Stuff\\Personal Learning\\Monte Carlo\\Lib80x\\xsdir.txt", ',')

# Define materials
#=
water = define_Material(
        density = 1.0,
        nuclides = ["O-nat", "H-nat"],
        atomic_fracs = [1.0, 2.0])

UO2 = define_Material(
        density = 10.0,
        nuclides = ["O-nat", "U-234"],
        atomic_fracs = [2.0, 1.0])

=#

B4C = define_Material(
        density = 2.52,
        nuclides = ["B-nat", "C-nat"],
        atomic_fracs = [4.0, 1.0])

Gd2O3 = define_Material(
        density = 7.40,
        nuclides = ["Gd-nat", "O-nat"],
        atomic_fracs = [2.0, 3.0])

#materials = Material[water, UO2]
materials = Material[B4C, Gd2O3]

# Create a list of incident neutron data needed
xs_nuclides = Tuple{Int64, Float64}[] # Pairs of ZA and T.
for i = 1:length(materials) # loop over materials
    for j = 1:length(materials[i].nuclides) # loop over nulcides in material
        ZA_T = (materials[i].nuclides[j], materials[i].temperature)
        # Add the nulcide identifier if it's not in the list yet
        if !(ZA_T in xs_nuclides)
            push!(xs_nuclides, ZA_T)
        end
    end
end

# Extract the paths to the data files needed. Keep nuclide ordering consistent
xs_nuclides_new = Tuple{Int64, Float64}[] # new ordering of xs_nuclides
xs_paths = String[] # paths to data files
for i = 1:length(xsdir[:,1])
    ZA = xsdir[i, 1]
    T = xsdir[i, 2]
    ZA_T = (ZA, T)
    if ZA_T in xs_nuclides
        push!(xs_nuclides_new, ZA_T)
        push!(xs_paths, xsdir[i, 3])
    end
end
xs_nuclides = xs_nuclides_new

# Load the nuclear data and populate the energy grid map
xs_data = AceData[] # list containing nuclear data
grid_map = zeros(Int64, (length(xs_nuclides), n_grid)) # energy grid map
for i = 1:length(xs_nuclides)
    acedata = read_ace(xs_paths[i])
    push!(xs_data, acedata)
    println(xs_nuclides[i])
    nuclide_E_grid = E_grid(acedata)
    for j = 1:n_grid
        grid_E = grid_energy(j)
        grid_map[i, j] = findlast(E -> E < grid_E, nuclide_E_grid)
    end

end


function ??_T_array(material::Material, E::Number)
    nuclide_idxs = findall(ZA_T -> (ZA_T[1] in material.nuclides) & (ZA_T[2] == material.temperature), xs_nuclides)
    map_lo = grid_index(E)
    idxs_lo = grid_map[nuclide_idxs, map_lo]
    idxs_hi = grid_map[nuclide_idxs, map_lo + 1] .+ 1
    idxs = find_idx.(xs_data[nuclide_idxs], E, idxs_lo, idxs_hi)
    interp_fracs = interp_frac_log.(xs_data[nuclide_idxs], E, idxs)
    ??_T0 = ??_T_idx.(xs_data[nuclide_idxs], idxs)
    ??_T1 = ??_T_idx.(xs_data[nuclide_idxs], idxs.+1)
    ??_Ts = material.atomic_dens.*exp.(log.(??_T0) + log.(??_T1./??_T0).*interp_fracs) * 1.0E-24
    return ??_Ts
end

function ??_a_array(material::Material, E::Number)
    nuclide_idxs = findall(ZA_T -> (ZA_T[1] in material.nuclides) & (ZA_T[2] == material.temperature), xs_nuclides)
    map_lo = grid_index(E)
    idxs_lo = grid_map[nuclide_idxs, map_lo]
    idxs_hi = grid_map[nuclide_idxs, map_lo + 1] .+ 1
    idxs = find_idx.(xs_data[nuclide_idxs], E, idxs_lo, idxs_hi)
    interp_fracs = interp_frac_log.(xs_data[nuclide_idxs], E, idxs)
    ??_a0 = ??_a_idx.(xs_data[nuclide_idxs], idxs)
    ??_a1 = ??_a_idx.(xs_data[nuclide_idxs], idxs.+1)
    ??_as = material.atomic_dens.*exp.(log.(??_a0) + log.(??_a1./??_a0).*interp_fracs) * 1.0E-24
    return ??_as
end

#=
Energy = 10 .^ collect(-10:0.001:1)
??_T = sum.(??_T_array(water, Energy[i]) for i = 1:length(Energy))
using PyPlot
loglog(Energy, ??_T)
=#

using PyPlot
Energy = 10 .^ collect(-10:0.001:1)
??_a = sum.(??_a_array(B4C, Energy[i]) for i = 1:length(Energy))
loglog(Energy, ??_a)
??_a = sum.(??_a_array(Gd2O3, Energy[i]) for i = 1:length(Energy))
loglog(Energy, ??_a)
