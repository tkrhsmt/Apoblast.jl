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
    z = Dict{Any,Any}[]
    for term in terms_after
        push!(z, _classify(term))
    end

    # Collect the unique terms from the original terms that appear in the classified terms
    g = Sym[]
    for term_dict in z
        for (key, _) in term_dict
            if key isa Sym && !(key in terms_before) && !(key in g)
                push!(g, key)
            end
        end
    end
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
    terms_after = []
    for term in terms_before
        push!(terms_after, apply_trans(model, trans, term))
    end
    terms_after = Sym.(terms_after)

    # Compute the transformation matrix T, the coefficient matrix D for the leak terms, and the list of leak terms g
    T, D, g = coeff_matrix(terms_after, terms_before)

    return T, D, g
end

# ==============================================================================
"""
    infinitesimal_trans_matrix(model::Model, trans::Transformation, parameter::Sym, terms_before::Vector{<:Sym})

Applies an infinitesimal transformation to a list of symbolic expressions (terms) and computes the transformation matrix `T`, the coefficient matrix `D` for any leak terms, and the list of leak terms `g`. The function first applies the transformation to each term in `terms_before`, differentiates the transformed terms with respect to the specified parameter, simplifies the result, and then uses the `coeff_matrix` function to compute the matrices and leak terms.

# Parameters
- `model`: An instance of the Model struct containing defined coordinates and fields.
- `trans`: An instance of the Transformation struct containing the transformation rules.
- `parameter`: A symbolic variable representing the parameter with respect to which the differentiation will be performed.
- `terms_before`: A vector of symbolic expressions representing the terms before transformation.
"""
function infinitesimal_trans_matrix(model::Model, trans::Transformation, parameter::Sym, terms_before::Vector{<:Sym})

    terms_after = Sym[]
    for term in terms_before
        transformed = apply_trans(model, trans, term)
        differentiated = diff(transformed, parameter)
        push!(terms_after, simplify(differentiated(trans.parameter_fixed)))
    end

    T, D, g = coeff_matrix(terms_after, terms_before)

    return T, D, g
end

# ==============================================================================

function make_K_matrix!(K, T, X, D)
    K_row, K_col = size(K)
    K_row_count = 1

    X_row, X_col = size(X)
    T_row, T_col = size(T)
    D_row, D_col = size(D)

    # R T - X R = 0 -> K Rv = 0
    # R : unknown coefficient matrix, Rv : vectorized form of R
    for a in 1:X_row
        for b in 1:T_col

            # R T -> K Rv
            for c in 1:T_row
                K[K_row_count, c+(a-1)*T_col] += T[c, b]
            end

            # - X R -> K Rv
            for d in 1:X_col
                K[K_row_count, b+(d-1)*T_col] -= X[a, d]
            end

            K_row_count += 1
        end
    end

    # R D = 0 -> K Rv = 0
    for a in 1:X_row
        for b in 1:D_col

            # R D -> K Rv
            for c in 1:D_row
                K[K_row_count, c+(a-1)*D_row] += D[c, b]
            end

            K_row_count += 1
        end
    end

    return nothing
end

function sympy_eye(n::Int)

    tmp = zeros(Sym, n, n)
    for i in 1:n
        tmp[i, i] = Sym(1)
    end
    return tmp
end

function _apply_constraint!(Kb, L, L⁺, Li, f̃, f, model::Model, comp_trans::Transformation, trans_builder, span_error::AbstractString; check_rank::Bool=false)
    T̃, D̃, _ = trans_builder(model, comp_trans, f̃)
    X = L * T̃ * L⁺
    check_rank && @assert rank(X) == length(Li) "Transformed Li must remain linearly independent"
    @assert all(iszero, L * D̃) span_error

    T, D, g = trans_builder(model, comp_trans, f)

    K_row = length(Li) * (length(f) + length(g))
    K_col = length(Li) * length(f)
    K = zeros(Sym, K_row, K_col)
    make_K_matrix!(K, T, X, D)

    K = K * Kb
    N = hcat(K.nullspace()...)
    return Kb * N
end

function coeff_basis(model::Model, library::Library, trans::Tuple{Vararg{Transformation}}, Li::Vector{<:Sym}, f̃::Vector{<:Sym})

    L, _, gd = coeff_matrix(Li, f̃)
    @assert gd == Sym[] "Li must not include any leak terms that are not in f̃"
    L⁺ = L.pinv()
    f = [library.term...]

    Kb = sympy_eye(length(Li) * length(f))
    for comp_trans in trans

        if comp_trans.parameter == ()
            Kb = _apply_constraint!(Kb, L, L⁺, Li, f̃, f, model, comp_trans, trans_matrix, "Transformed Li must stay within span(f̃)"; check_rank=true)
        else
            for parameter in comp_trans.parameter
                builder = (model, trans, terms) -> infinitesimal_trans_matrix(model, trans, parameter, terms)
                Kb = _apply_constraint!(Kb, L, L⁺, Li, f̃, f, model, comp_trans, builder, "Infinitesimal action on Li must stay within span(f̃)")
            end
        end

    end

    return Kb
end

function collect_follower(model::Model, library::Library, Li::Vector{<:Sym}, f̃::Vector{<:Sym}, trans::Transformation...)
    Kb = coeff_basis(model, library, trans, Li, f̃)
    Kb_row, Kb_col = size(Kb)

    output = []
    f = [library.term...]
    f_length = length(f)
    for i in 1:Kb_col
        v = Kb[:, i]
        m = permutedims(reshape(v, f_length, Kb_row ÷ f_length))
        push!(output, m * f)
    end
    return output
end

end
