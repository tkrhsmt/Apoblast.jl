module Param

export Model, Library

using SymPy
using PyCall
using ..Utils

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

        # ex. coords = ("x", "y") -> coords_sym = (x, y)
        coords_sym = Tuple(map(c -> Sym(c), coords))
        # ex. fields = ("u", "v") -> fields_sym = (u(x, y), v(x, y))
        fields_sym = Tuple(map(f -> SymFunction(f)(coords_sym...), fields))

        new(coords_sym, fields_sym)
    end
end

"""
Library struct represents a collection of terms (symbolic expressions) that are validated against a given Model. The constructor takes a Model instance and a tuple of symbolic expressions, validates the expressions using the Utils.validate_terms function, and initializes the Library struct.

# Parameters
- `model`: An instance of the Model struct containing defined coordinates and fields.
- `term`: A tuple of symbolic expressions that are validated to ensure they only contain symbols from the model's coordinates and fields.
"""
struct Library
    term::Tuple{Vararg{Sym}}

    function Library(model::Model, term::Tuple{Vararg{Sym}})

        # validate that the terms only contain symbols from model.coords and model.fields
        Utils.validate_terms(model, term)

        new(term)
    end
end

end
