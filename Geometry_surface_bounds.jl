include("C:\\Users\\warner\\Desktop\\My Stuff\\Personal Learning\\Monte Carlo\\NautilusMC\\NautilusMC\\Geometry.jl")

function surface_bounds(surf::Surface)


    if surf.type == "inf"
        return [-Inf, -Inf, -Inf], [Inf, Inf, Inf]


    elseif surf.type == "px"
        X = params[1]
        return [X, -Inf, -Inf], [X, Inf, Inf]


    elseif surf.type == "py"
        Y = params[1]
        return [-Inf, Y, -Inf], [Inf, Y, Inf]


    elseif surf.type == "pz"
        Z = params[1]
        return [-Inf, -Inf, Z], [Inf, Inf, Z]


    elseif surf.type == "planeEQ"
        return [-Inf, -Inf, -Inf], [Inf, Inf, Inf]


    elseif surf.type == "planePts"
        return [-Inf, -Inf, -Inf], [Inf, Inf, Inf]



    elseif surf.type == "box"
        x1, x2, y1, y2, z1, z2 = surf.params
        # Sort x
        if x1 < x2
            lo_x = x1
            hi_x = x2
        else
            lo_x = x2
            hi_x = x1
        end
        # Sort y
        if y1 < y2
            lo_y = y1
            hi_y = y2
        else
            lo_y = y2
            hi_y = y1
        end
        # Sort z
        if z1 < z2
            lo_z = z1
            hi_z = z2
        else
            lo_z = z2
            hi_z = z1
        end
        return [lo_x, lo_y, lo_z], [hi_x, hi_y, hi_z]


    end

end
