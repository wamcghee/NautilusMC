include("C:\\Users\\warner\\Desktop\\My Stuff\\Personal Learning\\Monte Carlo\\NautilusMC\\NautilusMC\\Geometry.jl")

function surface_distance(surf::Surface, x, y, z, u, v, w)

    x, y, z, u, v, w = trans_backwards(surf.transformations, x, y, z, u, v, w)

    if surf.type == "inf"
        return Inf


    elseif surf.type == "px"
        X = params[1]
        d = (X - x) / u
        if d < 0
            return Inf
        else
            return d
        end


    elseif surf.type == "py"
        Y = params[1]
        d = (Y - y) / v
        if d < 0
            return Inf
        else
            return d
        end


    elseif surf.type == "pz"
        Z = params[1]
        d = (Z - z) / w
        if d < 0
            return Inf
        else
            return d
        end


    elseif surf.type == "planeEQ"
        A, B, C, D = params
        d = -(D + A*x + B*y + C*z) / (A*u + B*v + C*w)
        if d < 0
            return Inf
        else
            return d
        end


    elseif surf.type == "planePts"
        x1, y1, z1, x2, y2, z2, x3, y3, z3 = params
        A = y2*z3 - y3*z2 - y1*(z3 - z2) + z1*(y3 - y2)
        B = z2*x3 - z3*x2 - z1*(x3 - x2) + x1*(z3 - z2)
        C = x2*y3 - x3*y2 - x1*(y3 - y2) + y1*(x3 - x2)
        D = -A*x1 - B*x2 - C*x3
        d = -(D + A*x + B*y + C*z) / (A*u + B*v + C*w)
        if d < 0
            return Inf
        else
            return d
        end


    elseif surf.type == "box"
        x1, x2, y1, y2, z1, z2 = surf.params
        dx1 = (x1 - x) / u
        dx2 = (x2 - x) / u
        dy1 = (y1 - y) / v
        dy2 = (y2 - y) / v
        dz1 = (z1 - z) / w
        dz2 = (z2 - z) / w
        distances = [dx1, dx2, dy1, dy2, dz1, dz2]
        distances = ifelse.(distances .< 0, Inf, distances)
        return minimum(distances)
    end

end
