using DelimitedFiles

# Unspecified scientific constancs
include("C:\\Users\\warner\\Desktop\\My Stuff\\Personal Learning\\Monte Carlo\\constants.jl")

# Table of masses of nuclides
masses = readdlm("C:\\Users\\warner\\Desktop\\My Stuff\\Personal Learning\\Monte Carlo\\masses.txt")
atomic_weight = Dict()
for i = 1:size(masses)[1]
    atomic_weight[convert(Int64, masses[i, 1])] = masses[i, 2]
end
masses = 0

# Table of natural abundances of nuclides
abundances = readdlm("C:\\Users\\warner\\Desktop\\My Stuff\\Personal Learning\\Monte Carlo\\abundances.txt")

# Table of atomic numbers and symbols
symbols = readdlm("C:\\Users\\warner\\Desktop\\My Stuff\\Personal Learning\\Monte Carlo\\atomic_symbols.txt")


function to_ZA(nuclide_str::String)
    #=
    Converts a string representing a nuclide into an Int64 in ZA notation

    :param nuclide_str: string representing a nuclide. "92235" or "U-235"
    :return ZA: Int64 representing the nuclide in ZA notation. 1000*Z + A
    =#
    if isnumeric(nuclide_str[1])
        ZA =  parse(Int64, nuclide_str)
    else
        hyphen = findfirst("-", nuclide_str)[]
        symb = nuclide_str[1 : hyphen-1]
        Z = findfirst(x -> x == symb, symbols[:,2])
        A_str = nuclide_str[hyphen+1 : end]
        if A_str == "nat"
            A = 0
        else
            A = parse(Int64, A_str)
        end
        ZA = 1000*Z + A
    end
    return ZA
end


function natural_nuclide(ZA::Int)
    #=
    Gather information needed to implement natural isotopes

    :param ZA: Int representing the nuclide in ZA notation. 1000*Z + A
    :return natural_weight: Float, average mass of atom given natural abundances
    :return nuclides: list of Int for ZA of nuclides present
    :fracs: List of Float for fractions of each nuclide naturally occuring
    =#
    Z = ZA ÷ 1000
    indices = findall(x -> x == Z, abundances[:, 1])
    As = abundances[indices, 2]
    nuclides = convert.(Int64, As.+Z*1000)
    fracs = abundances[indices, 4]
    weights = abundances[indices, 3]
    natural_weight = sum(fracs.*weights)
    mass_fracs = fracs.*weights / natural_weight
    return natural_weight, nuclides, fracs
end


# Add natural isotope atomic weights to table
for Z in 1:92
    atomic_weight[Z*1000] = natural_nuclide(Z*1000)[1]
end


mutable struct Material
    #=
    Structure to contain information about a material

    nuclides: List of Int to represent ZA of each nuclide in material
    atomic_dens: List of Float for the number denisity of each nuclide (1/cm^3)
    temperature: Float of temperature of material in K
    =#
    nuclides::Vector{Int64}
    atomic_dens::Vector{Float64}
    temperature::Float64
    color::Vector{Int64}
end


function define_Material(;
    density::Float64 = 0.0,
    number_density::Float64 = 0.0,
    nuclides::Vector{String} = String[],
    mass_fracs::Vector{Float64} = Float64[],
    mass_dens::Vector{Float64} = Float64[],
    atomic_fracs::Vector{Float64} = Float64[],
    atomic_dens::Vector{Float64} = Float64[],
    temperature::Float64 = 293.6,
    normalize_fracs::Bool = true,
    color::Vector{Int64} = [-1, -1, -1])
    #=
    This function allows for a variety of methods to define a material. As long
    as the provided information is sufficient to deduce the namber densities of
    each nuclide, the function will work unless there is conflicting
    information.

    :param density: density of material in g/cm^3
    :param number_density: number density of material in atoms/cm^3
    :param nuclides: List of nuclides present. Format may be "92235" or "U-235".
    Natural isotopes are indicated as "92000" or "U-nat".
    :param mass_fracs: fraction of mass from each nuclide
    :param mass_dens: density of each nuclide in g/cm^3
    :param atomic_fracs: fraction of total atoms from each nuclide
    :param atomic_dens: List of the number denisity of each nuclide (1/cm^3)
    :param temperature: temperature of material in K
    :param normalize_fracs: if false and fractions don't add to 1, the total
    number density or density will be inconsistent with the input. Setting to
    true will maintain nuclide ratios and the input total density or number
    density.
    :return material: Material struct based on the inputs
    =#
     if sum(length.([mass_fracs, mass_dens, atomic_fracs, atomic_dens]).>0) > 1
         throw(ArgumentError("Confliciting information in define_Material. Only
         one of mass_fracs, mass dens, atomic_fracs, or atomic_dens can be
         provided."))
     end

     nuclides = to_ZA.(nuclides)

     if length(mass_fracs) > 0

         if length(mass_fracs) != length(nuclides)
             throw(ArgumentError("Mismatch in define_Material. Number of nuclides
             and number of mass fractions provided must match."))
         end

         if !(density > 0)
             throw(ArgumentError("Incomplete information in define_Material. If
             using mass fractions, total density must be provided."))
         end

         if normalize_fracs
             mass_dens = mass_fracs * density / sum(mass_fracs)
         else
             mass_dens = mass_fracs * density
         end
     end

     if length(mass_dens) > 0

         if length(mass_dens) != length(nuclides)
             throw(ArgumentError("Mismatch in define_Material. Number of nuclides
             and number of mass densities provided must match."))
         end

         atomic_weights = [atomic_weight[nuclide] for nuclide in nuclides]
         atomic_dens = N_A * mass_dens ./ atomic_weights
     end

     if length(atomic_fracs) > 0

         if length(atomic_fracs) != length(nuclides)
             throw(ArgumentError("Mismatch in define_Material. Number of nuclides
             and number of atomic_dens fractions provided must match."))
         end

         if !(number_density > 0) & !(density > 0)
             throw(ArgumentError("Incomplete information in define_Material. If
             using atomic fractions, total number density or density must be provided."))
         end

         if (number_density > 0) & (density > 0)
             throw(ArgumentError("Conflicting information in define_Material. Only
             one of density and number density can be provided."))
         end

         if density > 0
             atomic_weights = [atomic_weight[nuclide] for nuclide in nuclides]
             number_density = sum(atomic_fracs) * N_A * density / sum(atomic_weights.*atomic_fracs)
         end

         if normalize_fracs
             atomic_dens = number_density * atomic_fracs / sum(atomic_fracs)
         else
             atomic_dens = number_density * atomic_fracs
         end
     end

     # Replace natural nuclides with their constituents
     full_nuclides = []
     full_atomic_dens = []
     for i = 1:length(nuclides)
         if nuclides[i] % 1000 == 0
             AW, ZA, fracs = natural_nuclide(nuclides[i])
             for item in ZA
                 push!(full_nuclides, item)
             end
             for item in fracs
                 push!(full_atomic_dens, atomic_dens[i]*item)
             end
         else
             push!(full_nuclides, nuclides[i])
             push!(full_atomic_dens, atomic_dens[i])
         end
     end

    order = sortperm(full_nuclides)

    return Material(full_nuclides[order],
                    full_atomic_dens[order],
                    temperature,
                    color)
end


function density(material::Material)
    #=
    Returns the density of a material in g/cm^3

    :param material: Material struct
    =#
    dens = 0
    for i = 1:length(material.nuclides)
        dens += material.atomic_dens[i] * atomic_weight[material.nuclides[i]] / N_A
    end
    return dens
end


function mass_frac(material::Material, nuclide)
    #=
    Returns the fraction of a material's mass from a given nuclide

    :param material: Material struct
    :param nuclide: Nuclide identifier String or ZA Int
    =#
    if typeof(nuclide) == String
        nuclide = to_ZA(nuclide)
    end
    nuclide_dens = 0
    dens = 0
    for i = 1:length(material.nuclides)
        partial_dens = material.atomic_dens[i] * atomic_weight[material.nuclides[i]] / N_A
        dens += partial_dens
        if material.nuclides[i] == nuclide
            nuclide_dens += partial_dens
        end
    end
    return nuclide_dens / dens
end


function mass_frac_elem(material::Material, element)
    #=
    Returns the fraction of a material's mass from a given element

    :param material: Material struct
    :param element: String or Int to identify element. "U" or 92.
    =#
    if typeof(element) == String
        element = findfirst(x -> x == element, symbols[:,2])
    end
    elem_dens = 0
    dens = 0
    for i = 1:length(material.nuclides)
        partial_dens = material.atomic_dens[i] * atomic_weight[material.nuclides[i]] / N_A
        dens += partial_dens
        if material.nuclides[i] ÷ 1000 == element
            elem_dens += partial_dens
        end
    end
    return elem_dens / dens
end


function mass_dens(material::Material, nuclide)
    #=
    Returns the density of a nuclide in a material in g/cm^3

    :param material: Material struct
    :param nuclide: Nuclide identifier String or ZA Int
    =#
    if typeof(nuclide) == String
        nuclide = to_ZA(nuclide)
    end
    i = findfirst(x -> x == nuclide, material.nuclides)
    if i == nothing
        return 0
    else
        return atomic_weight[nuclide] * material.atomic_dens[i] / N_A
    end
end


function mass_dens_elem(material::Material, element)
    #=
    Returns the density of an element in a material in g/cm^3

    :param material: Material struct
    :param element: String or Int to identify element. "U" or 92.
    =#
    if typeof(element) == String
        element = findfirst(x -> x == element, symbols[:,2])
    end
    elem_dens = 0
    for i = 1:length(material.nuclides)
        if material.nuclides[i] ÷ 1000 == element
            elem_dens += material.atomic_dens[i] * atomic_weight[material.nuclides[i]] / N_A
        end
    end
    return elem_dens
end
