module Param

export Model, Library, Transformation

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

"""
Transformation struct represents a transformation of coordinates and fields in a given Model. The constructor takes a Model instance, tuples of new coordinates and fields, optional parameters, and an optional Jacobian inverse. It creates replacement dictionaries for the coordinates and fields, computes the Jacobian inverse if not provided, and initializes the Transformation struct.

# Parameters
- `model`: An instance of the Model struct containing defined coordinates and fields.
- `coords_replace`: A tuple of symbolic variables representing the new coordinates after transformation.
- `fields_replace`: A tuple of symbolic functions representing the new fields after transformation.
- `parameter`: An optional tuple of symbolic variables representing parameters that may be involved in the transformation (default is an empty tuple).
- `jacobian_inverse`: An optional symbolic expression representing the inverse of the Jacobian matrix of the transformation. If not provided, it will be computed based on the new coordinates and the original coordinates.
"""
struct Transformation
    coords_replace::Tuple{Vararg{Sym}}
    fields_replace::Tuple{Vararg{Sym}}
    coords_replace_dict::Dict{Sym,Sym}
    fields_replace_dict::Dict{Sym,Sym}
    parameter::Tuple{Vararg{Sym}}
    jacobian_inverse::Array{Sym}

    function Transformation(
        model::Model,
        coords_replace::Tuple{Vararg{Sym}},
        fields_replace::Tuple{Vararg{Sym}};
        parameter::Tuple{Vararg{Sym}}=(),
        jacobian_inverse::Union{Nothing, Sym}=nothing,
    )

        # create replacement dictionaries for coordinates and fields
        coords_replace_dict = Dict(zip(model.coords, coords_replace))
        fields_replace_dict = Dict(zip(model.fields, fields_replace))

        # if the inverse of the Jacobian matrix is not provided, compute it
        if jacobian_inverse === nothing
            # compute the Jacobian matrix of the transformation
            jacobian_matrix = SymPy.Sym[coords_replace...].jacobian(SymPy.Sym[model.coords...])
            # compute the inverse of the Jacobian matrix
            jacobian_inverse = SymPy.simplify(SymPy.inv(jacobian_matrix))
        end

        new(coords_replace, fields_replace, coords_replace_dict, fields_replace_dict, parameter, jacobian_inverse)

    end
end

end
