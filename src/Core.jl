module Core

using SymPy
using PyCall
using ..Utils
using ..Param

function _is_derivative(expr::Sym)
    class_key = expr.class_key()
    return length(class_key) >= 3 && string(class_key[3]) == "Derivative"
end

function _is_plain_function_application(model::Model, expr::Sym)
    Bool(expr.is_Atom) && return false
    _is_derivative(expr) && return false
    expr in model.fields && return false
    return true
end

function _apply_derivative_rule(model::Model, trans::Transformation, expr::Sym)
    differentiated = apply_trans(model, trans, expr.args[1])

    result = differentiated
    for variable_count in expr.args[2:end]
        coord = variable_count.args[1]
        count = Int(variable_count.args[2])
        coord_index = findfirst(isequal(coord), model.coords)
        coord_index === nothing && return expr

        for _ in 1:count
            result = sum(
                trans.jacobian_inverse[coord_index, j] * diff(result, model.coords[j]) for
                j in eachindex(model.coords)
            )
        end
    end

    return result
end

"""
    apply_trans(model::Model, trans::Transformation, term::Sym)

Applies the given transformation to the specified term (symbolic expression) based on the rules defined in the Transformation struct. The function handles different cases for plain function applications, coordinate replacements, field replacements, and derivatives, recursively applying the transformation as needed.

# Parameters
- `model`: An instance of the Model struct containing defined coordinates and fields.
- `trans`: An instance of the Transformation struct containing the transformation rules.
- `term`: A symbolic expression (Sym) to which the transformation will be applied.
"""
function apply_trans(model::Model, trans::Transformation, term::Sym)

    # If s = f(e_1, ..., e_k) for an ordinary scalar-valued function f,
    # recursively apply the transformation to each argument first.
    if _is_plain_function_application(model, term)
        transformed_args = map(arg -> apply_trans(model, trans, arg), term.args)
        return term.func(transformed_args...)
    end

    # If s is a coordinate x_i, apply the coordinate replacement rule.
    for (coord, replacement) in zip(model.coords, trans.coords_replace)
        term == coord && return replacement
    end

    # If s is a field u_i, apply the field replacement rule.
    for (field, replacement) in zip(model.fields, trans.fields_replace)
        term == field && return replacement
    end

    # If s is a derivative, recursively transform the inner expression and
    # then apply the Jacobian-based derivative replacement rule.
    _is_derivative(term) && return _apply_derivative_rule(model, trans, term)

    # If no rule matches, return the original expression.
    return term
end


end
