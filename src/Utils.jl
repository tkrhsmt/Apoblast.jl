module Utils

export validate_terms

using SymPy
using PyCall

"""
    validate_terms(model, terms)

Validates that the given terms only contain symbols from the model's coordinates and fields.

# Arguments
- `model`: An instance of the `Model` struct containing defined coordinates and fields.
- `terms`: A tuple of symbolic expressions to validate.
"""
function validate_terms(model, terms)
    allowed_coords = Set(model.coords)
    allowed_fields = Set(model.fields)

    for expr in terms

        # check for invalid symbols
        invalid_symbols = setdiff(Set(free_symbols(expr)), allowed_coords)
        isempty(invalid_symbols) || throw(ArgumentError(
            "term $(expr) contains symbols outside model.coords: $(collect(invalid_symbols))",
        ))

        # check for invalid field-like functions
        invalid_field = find_invalid_fieldlike_function(expr, model.coords, allowed_fields)
        invalid_field === nothing || throw(ArgumentError(
            "term $(expr) contains a field-like function not defined in model.fields: $(invalid_field)",
        ))
    end

    return nothing
end

"""
    find_invalid_fieldlike_function(expr, coords, allowed_fields)

Recursively checks if the expression contains any function that looks like a field (i.e., a function of the coordinates) but is not in the allowed fields.

# Arguments
- `expr`: The symbolic expression to check.
- `coords`: The coordinates of the model, used to identify field-like functions.
- `allowed_fields`: The set of allowed field symbols defined in the model.
"""
function find_invalid_fieldlike_function(expr, coords, allowed_fields)
    if expr isa Sym
        if Bool(expr.is_Function) && expr.args == coords && !(expr in allowed_fields)
            return expr
        end

        for arg in expr.args
            invalid_field = find_invalid_fieldlike_function(arg, coords, allowed_fields)
            invalid_field === nothing || return invalid_field
        end
    elseif expr isa Tuple
        for arg in expr
            invalid_field = find_invalid_fieldlike_function(arg, coords, allowed_fields)
            invalid_field === nothing || return invalid_field
        end
    end

    return nothing
end


end
