#=
Definition of the Atom structure.
=#


module atom

import ..data as ase_data


"""
Structure for representing a single atom.

### Parameters:

`symbol`: str or int
- Can be a chemical symbol (str) or an atomic number (int).
`position`: vector of 3 floats
- Atomic position.
`tag`: int
- Special purpose tag.
`momentum`: vector of 3 floats
- Momentum for atom.
`mass`: float
- Atomic mass in atomic units.
`magmom`: float or vector of 3 floats
- Magnetic moment.
`charge`: float
- Atomic charge.
"""
mutable struct Atom
    __slots__::Array{String}
    atoms::Any
    index::Union{Missing,Int64}
    data::Dict{String,Any}

    __names__::Dict{String,Any}

    function Atom(symbol::String;
                  position=[0., 0., 0.]::Array{Float64},
                  tag=missing::Union{Missing,Int64},
                  momentum=missing::Union{Missing,Array{Float64}},
                  mass=missing::Union{Missing,Float64},
                  magmom=missing::Union{Missing,Float64},
                  charge=missing::Union{Missing,Float64},
                  atoms=missing::Any,
                  index=missing::Union{Missing,Int64})

        __slots = ["data", "atoms", "index"]
        __data = Dict{String,Any}()
        __names__ = Dict{String,Any}(
        "position" => ("positions", zeros(3)),
        "number" => ("numbers", 0),
        "tag" => ("tags", 0),
        "momentum" => ("momenta", zeros(3)),
        "mass" => ("masses", missing),
        "magmom" => ("initial_magmoms", 0.0),
        "charge" => ("initial_charges", 0.0)
        )
        
        if typeof(symbol) == String
            __data["number"] = ase_data.atomic_numbers[symbol]
        else
            __data["number"] = symbol
        end

        __data["position"] = position
        __data["tag"] = tag
        __data["momentum"] = momentum
        __data["mass"] = mass

        if typeof(magmom) <: Number
            magmom = [magmom]
        end

        __data["magmom"] = magmom
        __data["charge"] = charge

        new(__slots, atoms, index, __data, __names__)
    end
end


"""Custom short-form printing of an Atom object"""
function Base.show(io::IO, a::Atom)
    s = "Atom('$(a.symbol)'"
    if a.atoms === missing
        s = string(s, ", $(a.position))")
    else
        s = string(s, ", index=$(a.index))")
    end
    print(io, s)
end

"""Custom long-form printing of an Atom object"""
function Base.show(io::IO, ::MIME"text/plain", a::Atom)
    s = "Atom('$(a.symbol)', $(a.position)"
    for name in ["tag", "momentum", "mass", "magmom", "charge"]
        value = get_raw(a, name)
        if !(value === missing)
            s = string(s, ", $name=$value")
        end
    end
    if a.atoms === missing
        s = string(s, ")")
    else
        s = string(s, ", index=$(a.index))")
    end
    print(io, s)
end


"""Cut reference to atoms object."""
function cut_reference_to_atoms!(a::Atom)
    for name in keys(a.__names__)
        a.data[name] = a.get_raw(name)
    end
    a.index = missing
    a.atoms = missing
end


"""Get name attribute, return missing if not explicitly set."""
function get_raw(a::Atom, name::String)
    if name == "symbol"
        return ase_data.chemical_symbols[get_raw(a, "number")]
    end

    if a.atoms === missing
        return a.data[name]
    end

    plural = a.__names__[name][1]
    
    if plural in a.atoms.arrays
        return a.atoms.arrays[plural][a.index]
    else
        return missing
    end
end


"""Get name attribute, return default if not explicitly set."""
function get(a::Atom, name::String)
    value = get_raw(a, name)
    if value === missing
        if name == "mass"
            value = ase_data.atomic_masses[a.number]
        else
            value = a.__names__[name][2]
        end
    end
    return value
end


"""Set name attribute to value."""
function set!(a::Atom, name::String, value::Any)
    if name == "symbol"
        name = "number"
        value = ase_data.atomic_numbers[value]
    end
    if a.atoms === missing
        @assert name in keys(a.__names__)
        a.data[name] = value
    else
        plural, default = a.__names__[name]
        if plural in a.atoms.arrays
            array = a.atoms.arrays[plural]
            if name == "magmom" && ndims(array) == 2
                @assert len(value) == 3
            end
            array[a.index] = value
        else
            if name == "magmom" && ndims(value) == 1
                array = zeros((length(a.atoms), 3))
            elseif name == "mass"
                array = ase_atoms.get_masses(a.atoms) # IS THIS MUTABLE?
            else
                if !(typeof(default) <: Array)
                    default = [default]
                end
                # IF SOMETHING BREAKS, THIS CAN BE ONE REASON WHY
                array = zeros((length(a.atoms),) + size(default))
            end
            array[a.index] = value
            ase_atoms.new_array!(a, plural, array)
        end
    end
end


"""Remove name attribute"""
function delete(a::Atom, name::String)
    if atoms == missing && !(name in ["number", "symbol", "position"])
        a.data[name] = missing
    end
end


function Base.setproperty!(a::Atom, p::Symbol, value::Any)
    if p === :scaled_position
        pos = atoms.cell.cartesian_positions(value)
        a.position = pos

    #atomproperty setter
    elseif p in [:symbol, :number, :position, :tag, :momentum, :mass, :magmom, :charge]
        name = String(p)
        set!(a, name, value)

    # abcproperty setter
    elseif p in [:a, :b, :c]
        if p === :a
            index = 1
        else
            p === :b ? index = 2 : index = 3
        end
        @assert a.scaled_position != missing
        a.scaled_position[index] = value

    # xyzproperty setter
    elseif p in [:x, :y, :z]
        if p === :x
            index = 1
        else
            p === :y ? index = 2 : index = 3
        end
        a.position[index] = value

    # Default. Note that Julia doesn't allow adding new fields the same way Python does
    else
        setfield!(a, p, value)
    
    end
end

function Base.getproperty(a::Atom, p::Symbol)
    if p === :scaled_position
        if a.atoms != missing
            pos = a.position
            spos = a.atoms.cell.scaled_positions([pos])
            return spos[1]
        else
            return missing
        end

    # atomproperty getter
    elseif p in [:symbol, :number, :position, :tag, :momentum, :mass, :magmom, :charge]
        name = String(p)
        return get(a, name)
    
    # abcproperty getter
    elseif p in [:a, :b, :c]
        if p === :a
            index = 1
        else
            p === :b ? index = 2 : index = 3
        end
        if a.scaled_position != missing
            return a.scaled_position[index]
        else
            return missing
        end

    # xyzproperty getter
    elseif p in [:x, :y, :z]
        if p === :x
            index = 1
        else
            p === :y ? index = 2 : index = 3
        end
        return a.position[index]

    # Default
    else
        return getfield(a, p)

    end
end

end # end module
