
function ENDF_interp(nbt, int, x, y, q)
    # https://t2.lanl.gov/nis/endf/intro09.html
    if q >= x[end]
        return y[end]
    elseif q < x[1]
        return y[1]
    end
    idx_1 = findfirst(X -> X > q, x)
    idx_0 = idx_1 - 1
    x_1 = x[idx_1]
    x_0 = x[idx_0]
    y_1 = y[idx_1]
    y_0 = y[idx_0]
    interp_idx = findfirst(X -> X > idx_0, nbt)
    interp = int[interp_idx]
    #println(x_1, y_1, x_0, y_0, interp)
    if interp == 1
        return y_0
    elseif interp == 2
        return y_0 + (y_1 - y_0) * (q - x_0) / (x_1 - x_0)
    elseif interp == 3
        return y_0 + (y_1 - y_0) * log(q/x_0) / log(x_1/x_0)
    elseif interp == 4
        return exp(log(y_0) + log(y_1/y_0) * (q - x_0) / (x_1 - x_0))
    elseif interp == 5
        return exp(log(y_0) + log(y_1/y_0) * log(q/x_0) / log(x_1/x_0))
    elseif interp == 6
        B = log((y_1*x_1)/(y_0*x_0)) / (1/sqrt(x_0-0) - 1/sqrt(x_1-0))
        A = exp(B/sqrt(x_0-0))*y_0*x_0
        return (A/B) * exp(-B/sqrt(q-0))
    end
end


function ENDF_interp_func(NBT, INT, ES, P)
    return E -> ENDF_interp(NBT, INT, ES, P, E)
end


#filename ="C:/Users/warner/Downloads/Lib80x/Lib80x/H/1001.800nc"
filename = "C:/Users/warner/Desktop/My Stuff/Personal Learning/Monte Carlo/Lib80x/Lib80x/C/6012.805nc"
file_lines = readlines(filename)
pos = 1

if length(split(file_lines[1])) == 3
    pos = pos + 2
end
HZ, AW, TZ, HD = split(file_lines[pos], " ", keepempty = false)
HK = file_lines[pos+1][1:70]
HM = file_lines[pos+1][71:80]
pos = pos + 2

IZAW = split(join(file_lines[pos:pos+3]))
IZ = parse.(Int, IZAW[1:2:end])
AW = parse.(Float64, IZAW[2:2:end])
pos = pos + 4

NXS = parse.(Int, split(join(file_lines[pos:pos+1])))
LEN_XSS = NXS[ 1] # Length of second block of data (XSS array)
ZA      = NXS[ 2] # 1000 * Z + A
NES     = NXS[ 3] # Number of energies
NTR     = NXS[ 4] # Number of reactions excluding elastic scattering
NR      = NXS[ 5] # Number of reactions having secondary neutrons excluding elastic scattering
NTRP    = NXS[ 6] # Number of photon production reactions
NTYPE   = NXS[ 7] # Number of particle types for which production data is given
NPCR    = NXS[ 8] # Number of delayed neutron precurser families
S       = NXS[ 9] # Excited state
Z       = NXS[10] # Atomic number
A       = NXS[11] # Atomic mass number
pos = pos + 2

JXS = parse.(Int, split(join(file_lines[pos:pos+3])))
pos = pos + 4
XSS = parse.(Float64, split(join(file_lines[pos:end])))

ESZ     = JXS[ 1] # Energy table
NU      = JXS[ 2] # Fission yield data
MTR     = JXS[ 3] # MT array
LQR     = JXS[ 4] # Q-value array
TYR     = JXS[ 5] # Reaction type array
LSIG    = JXS[ 6] # Table of cross section locators
SIG     = JXS[ 7] # Cross sections
LAND    = JXS[ 8] # Table of angular distribution locators
AND     = JXS[ 9] # Angular distributions
LDLW    = JXS[10] # Table of energy distribution locators
DLW     = JXS[11] # Energy distributions
GPD     = JXS[12] # Photon production data
MTRP    = JXS[13] # Photon production MT array
LSIGP   = JXS[14] # Table of photon production cross section locators
SIGP    = JXS[15] # Photon production cross sections
LANDP   = JXS[16] # Table of photon production angular distribution locators
ANDP    = JXS[17] # Photon production angular distributions
LDLWP   = JXS[18] # Table of photon production energy distribution locators
DLWP    = JXS[19] # Photon production energy distributions
YP      = JXS[20] # Table of yield multipliers
FIS     = JXS[21] # Total fission cross section
END     = JXS[22] # Last word of the conventional table (last word of photon production data)
LUNR    = JXS[23] # Probability tables
DNU     = JXS[24] # Delayed neutron yield data
BDD     = JXS[25] # Basic delayed neutron precursor data (decay constant’s, probabilities)
DNEDL   = JXS[26] # Table of delayed neutron energy distribution locators
DNED    = JXS[27] # Delayed neutron energy distributions
PTYPE   = JXS[30] # Particle type array
NTRO    = JXS[31] # Array containing the number of particle production reactions
NEXT    = JXS[32] # Table of particle production locators (IXS array)

# Read ESZ section
Energy      = XSS[ESZ : ESZ + NES - 1]
sigma_t     = XSS[ESZ + NES : ESZ + 2*NES - 1]
sigma_a     = XSS[ESZ + 2*NES : ESZ + 3*NES - 1]
sigma_el    = XSS[ESZ + 3*NES : ESZ + 4*NES - 1]
H_ave       = XSS[ESZ + 4*NES : ESZ + 5*NES - 1]

# Read NU section
function nu_function(KNU)
    if XSS[KNU] == 1
        NC = convert(Int, XSS[KNU + 1])
        C = XSS[KNU + 2 : KNU + 1 + NC]
        l = collect(1:1:NC)
        return E -> sum(C.*E^(l-1))
    elseif XSS[KNU] == 2
        NR = convert(Int, XSS[KNU + 1])
        NBT = convert.(Int, XSS[KNU + 2 : KNU + 1 + NR])
        INT = convert.(Int, XSS[KNU + 2 + NR : KNU + 1 + 2*NR])
        NE = convert(Int, XSS[KNU + 2 + 2*NR])
        ES = XSS[KNU + 3 + 2*NR : KNU + 2 + 2*NR + NE]
        nubar = XSS[KNU + 3 + 2*NR + NE : KNU + 2 + 2*NR + 2*NE]
        if NR == 0
            NBT = [NE + 1]
            INT = [2]
        end
        return E -> ENDF_interp(NBT, INT, ES, nubar, E)
    end
end

# If XSS[NU] > 1, then there is only one array given that starts at NU[2]
if NU > 0
    if XSS[NU] > 0
        if DNU == 0
            nu_total = nu_function(NU)
        else
            nu_prompt = nu_function(NU)
        end
    # Otherwise, prompt and total arrays are given
    else
        nu_prompt = nu_function(NU + 1)
        nu_total = nu_function(NU + 1 + convert(Int, abs(XSS[NU])))
    end
end
if DNU > 0
    nu_delayed = nu_function(DNU)
end

#=
function dnu_function(start)

end


if BDD > 0
    DNP_decay_constants = []
    DNP_probabilities = Dict()
    idx = BDD
    for group = 1:NPCR
        DEC = XSS[idx]
        NR = convert(Int, XSS[idx + 1])
        NBT = convert.(Int, XSS[idx + 2 : idx + 1 + NR])
        INT = convert.(Int, XSS[idx + 2 + NR : idx + 1 + 2*NR])
        NE = convert(Int, XSS[idx + 2 + 2*NR])
        ES = XSS[idx + 3 + 2*NR : idx + 2 + 2*NR + NE]
        P = XSS[idx + 3 + 2*NR + NE : idx + 2 + 2*NR + 2*NE]
        if NR == 0
            NBT = [NE + 1]
            INT = [2]
        end
        append!(DNP_decay_constants, DEC)
        DNP_probabilities[group] = ENDF_interp_func(NBT, INT, ES, P)
        idx = idx + 3 + 2*NR + 2*NE
    end
end
=#


#=
using PyPlot
loglog(Energy, sigma_el)
=#
