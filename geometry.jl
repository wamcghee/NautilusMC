
include("C:\\Users\\warner\\Desktop\\My Stuff\\Personal Learning\\Monte Carlo\\NautilusMC\\NautilusMC\\Materials.jl")

mutable struct Trans
    type::String
    params::Vector{Float64}
end


mutable struct Surface
    type::String
    params::Vector{Float64}
    transformations::Vector{Trans}
end


mutable struct Cell
    surfs_within::Vector{Surface}
    surfs_outside::Vector{Surface}
    fill::Any
    transformations::Vector{Trans}
end


mutable struct Universe
    cells::Vector{Cell}
    transformations::Vector{Trans}
end


function trans_backwards(transformations::Vector{Trans},
                         x::Float64,
                         y::Float64,
                         z::Float64,
                         u::Float64 = 10,
                         v::Float64 = 10,
                         w::Float64 = 10)
    if u == 10
        return x, y, z
    else
        return x, y, z, u, v, w
    end
end
