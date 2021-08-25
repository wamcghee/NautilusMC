
include("C:\\Users\\warner\\Desktop\\My Stuff\\Personal Learning\\Monte Carlo\\Materials.jl")

mutable struct Surface
    type::String
    bounds_func::Function
    inside_func::Function
    dist_to_edge_func::Function
end


struct Cell
    universe::Any
    surfs_within::Vector{Surface}
    surfs_outside::Vector{Surface}
    fill::Any
end


function inside_surf(
    surface::Surface,
    x::Vector{Float64})
    return surface.inside_func(x)
end

function dist_to_surf(
    surface::Surface,
    x::Vector{Float64},
    v::Vector{Float64})
    return surface.dist_to_edge_func(x, v)
end



function inf()
    function bounds_func()
        return [-Inf, -Inf, -Inf], [Inf, Inf, Inf]
    end

    function inside_func(x::Vector{Float64})
        return true
    end

    function dist_to_edge_func(x::Vector{Float64}, v::Vector{Float64})
            return Inf
    end

    return Surface("infinite", bounds_func, inside_func, dist_to_edge_func)

end


function px(X::Float64)
    function bounds_func()
        return [X, -Inf, -Inf], [X, Inf, Inf]
    end

    function inside_func(x::Vector{Float64})
        return (x[1] < X)
    end

    function dist_to_edge_func(x::Vector{Float64}, v::Vector{Float64})
        d = (X - x[1]) / v[1]
        if d > 0
            return d
        else
            return Inf
        end
    end

    return Surface("plane", bounds_func, inside_func, dist_to_edge_func)

end


function py(Y::Float64)
    function bounds_func()
        return [-Inf, Y, -Inf], [Inf, Y, Inf]
    end

    function inside_func(x::Vector{Float64})
        return (x[2] < Y)
    end

    function dist_to_edge_func(x::Vector{Float64}, v::Vector{Float64})
        d = (Y - x[2]) / v[2]
        if d > 0
            return d
        else
            return Inf
        end
    end

    return Surface("plane", bounds_func, inside_func, dist_to_edge_func)

end


function pz(Z::Float64)
    function bounds_func()
        return [-Inf, -Inf, Z], [Inf, Inf, Z]
    end

    function inside_func(x::Vector{Float64})
        return (x[3] < Z)
    end

    function dist_to_edge_func(x::Vector{Float64}, v::Vector{Float64})
        d = (Z - x[3]) / v[3]
        if d > 0
            return d
        else
            return Inf
        end
    end

    return Surface("plane", bounds_func, inside_func, dist_to_edge_func)

end


function planeEQ(A::Float64, B::Float64, C::Float64, D::Float64)
    #=
    Plane with equation Ax + By + Cz + D = 0
    =#
    function bounds_func()
        return [-Inf, -Inf, -Inf], [Inf, Inf, Inf]
    end

    function inside_func(x::Vector{Float64})
        val = A*x[1] + B*x[2] + C*x[3] + D
        return (val < 0)
    end

    function dist_to_edge_func(x::Vector{Float64}, v::Vector{Float64})
        d = -(D + A*x[1] + B*x[2] + C*x[3]) / (A*v[1] + B*v[2] + C*v[3])
        if d > 0
            return d
        else
            return Inf
        end
    end

    return Surface("plane", bounds_func, inside_func, dist_to_edge_func)

end


function planePts(x1::Vector{Float64}, x2::Vector{Float64}, x3::Vector{Float64})
    #=
    Plane through points x1, x2, and x3.
    =#

    u1 = x1 - x2
    u2 = x3 = x2

    A = u1[2]*u2[3] - u1[3]*u2[2]
    B = u1[1]*u2[3] + u1[3]*u2[1]
    C = u1[1]*u2[2] - u1[2]*u2[1]
    D = -A*x[1] - B*x[2] - C*x[3]

    function bounds_func()
        return [-Inf, -Inf, -Inf], [Inf, Inf, Inf]
    end

    function inside_func(x::Vector{Float64})
        vec = x - x2
        val = A*vec[1] + B*vec[2] + C*vec[3]
        return (val < 0)
    end

    function dist_to_edge_func(x::Vector{Float64}, v::Vector{Float64})
        d = -(D + A*x[1] + B*x[2] + C*x[3]) / (A*v[1] + B*v[2] + C*v[3])
        if d > 0
            return d
        else
            return Inf
        end
    end

    return Surface("plane", bounds_func, inside_func, dist_to_edge_func)

end


function box(x1::Vector{Float64}, x2::Vector{Float64})
    #=
    Rectangular prism drawn between corners x and x2
    =#

    function bounds_func()
        return x1, x2
    end

    function inside_func(x::Vector{Float64})
        return sum(x.>x1) + sum(x.<x2) == 6
    end

    function dist_to_edge_func(x::Vector{Float64}, v::Vector{Float64})
        dx1 = (x1[1] - x[1]) / v[1]
        dx2 = (x2[1] - x[1]) / v[1]
        dy1 = (x1[2] - x[2]) / v[2]
        dy2 = (x2[2] - x[2]) / v[2]
        dz1 = (x1[3] - x[3]) / v[3]
        dz2 = (x2[3] - x[3]) / v[3]
        distances = [dx1, dx2, dy1, dy2, dz1, dz2]
        distances = ifelse.(distances .< 0, Inf, distances)
        return minimum(distances)
    end

    return Surface("box", bounds_func, inside_func, dist_to_edge_func)

end
