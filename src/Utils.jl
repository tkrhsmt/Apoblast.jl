module Utils

export validate_terms, safe_print, listorder

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

"""
    safe_print(print_progress, str)

Prints the given string if `print_progress` is true.
"""
function safe_print(print_progress::Bool, str::String)
    print_progress && println(str)
end

"""
    listorder(monoterms, order)
Generates a list of all combinations of the given monoterms up to the specified order. For example, if `monoterms` is `[x, y]` and `order` is 2, it will generate `[1, x, y, x^2, x*y, y^2]`.

# Arguments
- `monoterms`: A vector of symbolic variables representing the monoterms to combine.
- `order`: An integer specifying the maximum order of combinations to generate.
"""
function listorder(monoterms::Vector, order::Int=2)
    isempty(monoterms) && return [1]

    return vcat([
        [(monoterms[1]^i) * t
         for t in listorder(monoterms[2:end], order - i)]
        for i in 0:order
    ]...)
end

"""
    listorder(monoterm, coords, order)
Generates a list of all combinations of derivatives of the given monoterm with respect to the specified coordinates up to the specified order. For example, if `monoterm` is `u`, `coords` is `[x, y]`, and `order` is 2, it will generate a list of derivatives like `[u, diff(u, x), diff(u, y), diff(u, x, x), diff(u, x, y), diff(u, y, y)].

# Arguments
- `monoterm`: A symbolic expression representing the monoterm to differentiate.
- `coords`: A vector of symbolic variables representing the coordinates with respect to which the derivatives will be taken.
- `order`: An integer specifying the maximum order of derivatives to generate.
"""
function listorder(monoterm::Sym, coords::Vector, order::Int=2)
    isempty(coords) && return monoterm

    return vcat([
        [diff(t, coords[1], i)
         for t in listorder(monoterm, coords[2:end], order - i)]
        for i in 0:order
    ]...)
end


end
