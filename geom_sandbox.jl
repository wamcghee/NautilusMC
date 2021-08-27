
include("C:\\Users\\warner\\Desktop\\My Stuff\\Personal Learning\\Monte Carlo\\NautilusMC\\NautilusMC\\Materials.jl")
include("C:\\Users\\warner\\Desktop\\My Stuff\\Personal Learning\\Monte Carlo\\NautilusMC\\NautilusMC\\Geometry.jl")
include("C:\\Users\\warner\\Desktop\\My Stuff\\Personal Learning\\Monte Carlo\\NautilusMC\\NautilusMC\\Geometry_surface_bounds.jl")
include("C:\\Users\\warner\\Desktop\\My Stuff\\Personal Learning\\Monte Carlo\\NautilusMC\\NautilusMC\\Geometry_surface_inside.jl")
include("C:\\Users\\warner\\Desktop\\My Stuff\\Personal Learning\\Monte Carlo\\NautilusMC\\NautilusMC\\Geometry_surface_distance.jl")

water = define_Material(
    density = 1.0,
    nuclides = ["O-nat", "H-nat"],
    atomic_fracs = [1.0, 2.0],
    color = [0, 200, 0])

B4C = define_Material(
    density = 2.52,
    nuclides = ["B-nat", "C-nat"],
    atomic_fracs = [4.0, 1.0],
    color = [100, 0, 100])

s1 = Surface("box", [-1.0, 1.0, -1.0, 1.0, -1.0, 1.0], [])
s2 = Surface("box", [-2.0, 2.0, -2.0, 2.0, -2.0, 2.0], [])
s3 = Surface("box", [-3.0, 3.0, -3.0, 3.0, -3.0, 3.0], [])

C1 = Cell([s1], [], water, [])
C2 = Cell([s2], [s1], B4C, [])
outside = Cell([], [s2], nothing, [])

Root = Universe([C1, C2, outside], [])

function get_cell(uni::Universe, x::Float64, y::Float64, z::Float64)
    for cell in uni.cells
        inside = true
        for surf in cell.surfs_within
            if ! surface_inside(surf, x, y, z)
                inside = false
                break
            end
        end
        for surf in cell.surfs_outside
            if surface_inside(surf, x, y, z)
                inside = false
                break
            end
        end
        if inside
            return cell
        end
    end
    throw(Error("No cell at position %f, %f, %f", x, y, z))
end


function cell_neighbors(surf::Surface)
