include("C:\\Users\\warner\\Desktop\\My Stuff\\Personal Learning\\Monte Carlo\\NautilusMC\\NautilusMC\\Geometry.jl")

function surface_inside(surf::Surface, x, y, z)

    x, y, z = trans_backwards(surf.transformations, x, y, z)

    if surf.type == "inf"
        return true


    elseif surf.type == "px"
        X = params[1]
        return x <= X


    elseif surf.type == "py"
        Y = params[1]
        return y <= Y


    elseif surf.type == "pz"
        Z = params[1]
        return z <= Z


    elseif surf.type == "planeEQ"
        A, B, C, D = params
        val = A*x + B*y + C*z + D
        return (val <= 0)


    elseif surf.type == "planePts"
        x1, y1, z1, x2, y2, z2, x3, y3, z3 = params
        A = y2*z3 - y3*z2 - y1*(z3 - z2) + z1*(y3 - y2)
        B = z2*x3 - z3*x2 - z1*(x3 - x2) + x1*(z3 - z2)
        C = x2*y3 - x3*y2 - x1*(y3 - y2) + y1*(x3 - x2)
        D = -A*x1 - B*x2 - C*x3
        val = A*x + B*y + C*z + D
        return (val <= 0)


    elseif surf.type == "box"
        x1, x2, y1, y2, z1, z2 = surf.params
        # X bounds
        if (x1 > x2)
            if x > x1
                return false
            elseif x < x2
                return false
            end
        else
            if x < x1
                return false
            elseif x > x2
                return false
            end
        end
        # Y bounds
        if (y1 > y2)
            if y > y1
                return false
            elseif y < y2
                return false
            end
        else
            if y < y1
                return false
            elseif y > y2
                return false
            end
        end
        # Z bounds
        if (z1 > z2)
            if z > z1
                return false
            elseif z < z2
                return false
            end
        else
            if z < z1
                return false
            elseif z > z2
                return false
            end
        end
        return true

    end

end
