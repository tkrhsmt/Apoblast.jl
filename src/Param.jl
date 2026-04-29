module Param

export Model

using SymPy
using PyCall

"""
Model struct represents a mathematical model with specified coordinates and fields.
The constructor takes tuples of strings for coordinates and fields, converts them to symbolic variables and functions, and initializes the Model struct.

# Parameters
- `coords`: A tuple of symbolic variables representing the coordinates (e.g., spatial dimensions).
- `fields`: A tuple of symbolic functions representing the fields defined over the coordinates (e.g., physical quantities like velocity, pressure, etc.).

"""
struct Model
    coords::Tuple{Vararg{Sym}}
    fields::Tuple{Vararg{Sym}}

    function Model(coords::Tuple{Vararg{String}}, fields::Tuple{Vararg{String}})

        coords_sym = Tuple(map(c -> Sym(c), coords))
        fields_sym = Tuple(map(f -> SymFunction(f)(coords_sym...), fields))

        new(coords_sym, fields_sym)
    end
end

end
