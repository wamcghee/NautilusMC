using DelimitedFiles

mutable struct AceData
    #=
    Structure to contain information found in ace files

    IZ: List of 16 nuclide ZA identifiers. Needed for  S(α, β).
    AW: List of 16 atomic weight ratios to go with IZ. Needed for  S(α, β).
    NXS: List of 16 integers defining features of the ACE data file.
    JXS: List of 32 integers locating different parts of the ACE data file.
    XSS: List of floats containing the majority of the nuclear data.
    =#
    IZ::Vector{Int64}
    AW::Vector{Float64}
    NXS::Vector{Int64}
    JXS::Vector{Int64}
    XSS::Vector{Float64}
end

function read_ace(filename::String)
    #=
    Creates an AceData struct based on the path to an ACE data file.

    :param filename: String for the path to the ACE data file.
    :return acedata: AceData struct containing all the information in the file.
    =#
    file_lines = readlines(filename)
    pos = 1

    if length(split(file_lines[1])) == 3
        pos = pos + 2
    end
    HZ, AW, TZ, HD = split(file_lines[pos], " ", keepempty = false)
    HK = file_lines[pos+1][1:73]
    HM = file_lines[pos+1][74:80]
    pos = pos + 2
    IZAW = split(join(file_lines[pos:pos+3]))
    IZ = parse.(Int, IZAW[1:2:end])
    AW = parse.(Float64, IZAW[2:2:end])
    pos = pos + 4
    NXS = parse.(Int, split(join(file_lines[pos:pos+1])))
    pos = pos + 2
    JXS = parse.(Int, split(join(file_lines[pos:pos+3])))
    pos = pos + 4
    XSS = parse.(Float64, split(join(file_lines[pos:end])))
    #XSS = vec(readdlm(filename, Float64, skipstart = pos-1))
    # Return data
    acedata = AceData(IZ, AW, NXS, JXS, XSS)
    return acedata
end

function E_grid(acedata::AceData)
    #=
    Returns the full common energy grid from an AceData struct.

    :param acedata: AceData struct containing the relevant information
    :return: commom energy grid in the ACE data.
    =#
    NES = acedata.NXS[3]
    return acedata.XSS[acedata.JXS[1] : acedata.JXS[1] + NES - 1]
end

function E_idx(acedata::AceData, idx::Int)
    #=
    Returns energy at a specified index in the grid.

    :param acedata: AceData struct containing relevant information
    :param idx: index of interset
    =#
    if idx > get_NES(acedata)
        throw(ArgumentError("Index out of grid range."))
    end
    return acedata.XSS[acedata.JXS[1] + idx - 1]
end

function σ_T_idx(acedata::AceData, idx::Int)
    #=
    Returns microscopic total cross section at a specified index in the grid.

    :param acedata: AceData struct containing relevant information
    :param idx: index of interset
    =#
    NES = get_NES(acedata)
    if idx > NES
        throw(ArgumentError("Index out of grid range."))
    end
    return acedata.XSS[NES + acedata.JXS[1] + idx - 1]
end

function σ_a_idx(acedata::AceData, idx::Int)
    #=
    Returns microscopic absorption cross section at a specified index in the
    grid.

    :param acedata: AceData struct containing relevant information
    :param idx: index of interset
    =#
    NES = get_NES(acedata)
    if idx > NES
        throw(ArgumentError("Index out of grid range."))
    end
    return acedata.XSS[2*NES + acedata.JXS[1] + idx - 1]
end

function σ_el_idx(acedata::AceData, idx::Int)
    #=
    Returns microscopic elastic scattering cross section at a specified index in
    the grid.

    :param acedata: AceData struct containing relevant information
    :param idx: index of interset
    =#
    NES = get_NES(acedata)
    if idx > NES
        throw(ArgumentError("Index out of grid range."))
    end
    return acedata.XSS[4*NES + acedata.JXS[1] + idx - 1]
end

function σ_f_idx(acedata::AceData, idx::Int)
    #=
    Returns microscopic total fission cross section at a specified index in the
    grid.

    :param acedata: AceData struct containing relevant information
    :param idx: index of interset
    =#
    F = acedata.JXS[21]
    if F == 0
        return 0
    else
        IE = convert(Int, acedata.XSS[F])
        if idx < IE
            return 0
        else
            NE = convert(Int, acedata.XSS[F + 1])
            if idx > NE + IE - 1
                throw(ArgumentError("Index out of grid range for FIS."))
            end
        end
        return acedata.XSS[F + acedata.JXS[1] + idx - IE]
    end
end

function heat_idx(acedata::AceData, idx::Int)
    #=
    Returns microscopic heating cross section at a specified index in the grid.

    :param acedata: AceData struct containing relevant information
    :param idx: index of interset
    =#
    NES = get_NES(acedata)
    if idx > NES
        throw(ArgumentError("Index out of grid range."))
    end
    return acedata.XSS[5*NES + acedata.JXS[1] + idx - 1]
end

function get_NES(acedata::AceData)
    #=
    Returns the length of the energy grid.

    :param acedata: AceData struct containing relevant information
    =#
    return acedata.NXS[3]
end

function find_idx(
    acedata::AceData,
    E::Float64,
    lo::Int = 1,
    hi::Int = get_NES(acedata))
    #=
    Returns the index of the enegy space directly below or equal to the given
    energy. Uses recursive binary search algorithm.

    :param acedata: AceData struct containing relevant information
    :param E: Energy of interest
    :param lo: Lower bound index of search
    param hi: Upper bound index of search
    =#
    # Base case to terminate search
    if hi - lo == 1
        return lo
    end
    # Split search down middle
    mid = (lo + hi) ÷ 2
    # Recursive cases if mid point is too high or low
    if E_idx(acedata, mid) > E
        return find_idx(acedata, E, lo, mid)
    else
        return find_idx(acedata, E, mid, hi)
    end
end

function nu_prompt(acedata::AceData, E::Float64)

end

function interp_frac_log(
    acedata::AceData,
    E::Number,
    idx::Int)
    #=
    Returns the logarithmic fraction between the adjacent energy points in the
    grid used for interpolation

    :param acedata: AceData struct containing relevant information
    :param E: Energy of interest
    :param idx: Energy grid index directly below the energy of interest
    =#
    E0 = E_idx(acedata, idx)
    E1 = E_idx(acedata, idx + 1)
    return log(E/E0) / log(E1/E0)
end
