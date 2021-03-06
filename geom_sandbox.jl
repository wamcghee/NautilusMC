
include("C:\\Users\\warner\\Desktop\\My Stuff\\Personal Learning\\Monte Carlo\\NautilusMC\\NautilusMC\\Materials.jl")
include("C:\\Users\\warner\\Desktop\\My Stuff\\Personal Learning\\Monte Carlo\\NautilusMC\\NautilusMC\\geometry.jl")

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

s1 = box([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0])
s2 = box([-2.0, -2.0, -2.0], [2.0, 2.0, 2.0])
s3 = box([-3.0, -3.0, -3.0], [3.0, 3.0, 3.0])

C1 = Cell([s1], [], water, [])
C2 = Cell([s2], [s1], B4C, [])
outside = Cell([], [s2], nothing, [])

Root = Universe([C1, C2, outside], [])
