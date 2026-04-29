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

# ==============================================================================

function _classify(term::Sym)

    # expand the term to get a sum of products form
    term_ex = expand(term)
    # get the coefficients of each term in the expanded expression
    term_dict = term_ex.as_coefficients_dict()

    return Dict(k => term_dict[k] for k in keys(term_dict))
end

"""
    coeff_matrix(terms_after::Vector{<:Sym}, terms_before::Vector{<:Sym})

Given two lists of symbolic expressions, `terms_after` and `terms_before`, this function computes the transformation matrix `T` that relates the two sets of terms, as well as a coefficient matrix `D` for any "leak" terms that appear in `terms_after` but not in `terms_before`. The function returns the transformation matrix `T`, the coefficient matrix `D`, and a list of leak terms `g`.

# Parameters
- `terms_after`: A vector of symbolic expressions representing the terms after transformation.
- `terms_before`: A vector of symbolic expressions representing the terms before transformation.
"""
function coeff_matrix(terms_after::Vector{<:Sym}, terms_before::Vector{<:Sym})

    # Get the number of terms (should be the same for both after and before)
    terms_len_after, terms_len_before = length(terms_after), length(terms_before)

    # Classify the terms after transformation to get their coefficients in terms of the original terms
    z = []
    for term in terms_after
        push!(z, _classify(term))
    end

    # Collect the unique terms from the original terms that appear in the classified terms
    g = []
    for term_dict in z
        for (key, _) in term_dict
            if !(key in terms_before)
                push!(g, key)
            end
        end
    end
    g = Set(g)
    g_len = length(g)

    # Construct the transformation matrix T
    T = zeros(Sym, terms_len_after, terms_len_before)
    for j in 1:terms_len_before
        for i in 1:terms_len_after
            T[i, j] = get(z[i], terms_before[j], Sym(0))
        end
    end

    # Construct the coefficient matrix D for the leak terms g
    D = zeros(Sym, terms_len_after, g_len)
    for j in 1:g_len
        for i in 1:terms_len_after
            D[i, j] = get(z[i], g[j], Sym(0))
        end
    end

    return T, D, g
end

# ==============================================================================
"""
    trans_matrix(model::Model, trans::Transformation, terms_before::Vector{<:Sym})

Applies the given transformation to a list of symbolic expressions (terms) and computes the transformation matrix `T`, the coefficient matrix `D` for any leak terms, and the list of leak terms `g`. The function first applies the transformation to each term in `terms_before` to get `terms_after`, and then uses the `coeff_matrix` function to compute the matrices and leak terms.

# Parameters
- `model`: An instance of the Model struct containing defined coordinates and fields.
- `trans`: An instance of the Transformation struct containing the transformation rules.
- `terms_before`: A vector of symbolic expressions representing the terms before transformation.
"""
function trans_matrix(model::Model, trans::Transformation, terms_before::Vector{<:Sym})

    # Apply the transformation to each term in terms_before to get terms_after
    terms_len = length(terms_before)
    terms_after = []
    for term in terms_before
        push!(terms_after, apply_trans(model, trans, term))
    end
    terms_after = Sym.(terms_after)

    # Compute the transformation matrix T, the coefficient matrix D for the leak terms, and the list of leak terms g
    T, D, g = coeff_matrix(terms_after, terms_before)

    return T, D, g
end

end
