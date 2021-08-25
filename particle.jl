include("C:\\Users\\warner\\Desktop\\My Stuff\\Personal Learning\\Monte Carlo\\NautilusMC\\NautilusMC\\Materials.jl")
include("C:\\Users\\warner\\Desktop\\My Stuff\\Personal Learning\\Monte Carlo\\NautilusMC\\NautilusMC\\get_XS.jl")


mutable struct Particle
    type::String
    x::Vector{Float64}
    v::Vector{Float64}
    E::Float64
    cell::Cell
    weight::Float64
    alive::Bool
end

function transport(particle)
    material = particle.cell.material
    max_distance = cell.max_distance(x, v)
    # Get total cross sections
    nuclide_idxs = findall(ZA_T -> (ZA_T[1] in material.nuclides) & (ZA_T[2] == material.temperature), xs_nuclides)
    map_lo = grid_index(E)
    idxs_lo = grid_map[nuclide_idxs, map_lo]
    idxs_hi = grid_map[nuclide_idxs, map_lo + 1] .+ 1
    idxs = find_idx.(xs_data[nuclide_idxs], E, idxs_lo, idxs_hi)
    interp_fracs = interp_frac_log.(xs_data[nuclide_idxs], E, idxs)
    σ_T0 = σ_T_idx.(xs_data[nuclide_idxs], idxs)
    σ_T1 = σ_T_idx.(xs_data[nuclide_idxs], idxs.+1)
    Σ_Ts = material.atomic_dens.*exp.(log.(σ_T0) + log.(σ_T1./σ_T0).*interp_fracs) * 1.0E-24
    Σ_T = sum(Σ_T)
    # Move particle
    distance = -log(rand()) / Σ_T
    particle.x += particle.v * distance
    # Determine collision nuclide
    nuclide_probs = cumsum(Σ_Ts) / Σ_T
    collision_nuclide = material.nuclides[findlast(p -> p < nuclide_probs)]
    # Determine collision reaction
    println(collision_nuclide)
    # Do the collision
end
