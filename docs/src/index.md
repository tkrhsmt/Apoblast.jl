```@meta
CurrentModule = Apoblast
```

# Apoblast

Documentation for [Apoblast](https://github.com/tkrhsmt/Apoblast.jl) (A Priori Operator-Based Linear Algorithm for Symmetric Terms).
This package is a Julia implementation of the [python package Apoblast](https://github.com/JyJyJcr/apoblast).

## Installation

```julia
using Pkg
Pkg.add(path="https://github.com/tkrhsmt/Apoblast.jl")
```

## Quick start

We can use Apoblast to collect followers of a transformation. For example, consider the following code:

```julia
using Apoblast
using SymPy

# Define a model
model = Model(("x", "y"), ())
x, y = model.coords

# Define a library
lib = Library(model, Tuple(listorder([x, y], 10)))

# Define a transformation
trans_1 = Transformation(model, (x, -y), ())
θ = Sym(2) / Sym(5) * Sym(pi)
trans_2 = Transformation(model, (x*cos(θ) + y*sin(θ), y*cos(θ) - x*sin(θ)), ())

# Define a LHS
Li = [Sym(1)]
f̃ = [Sym(1)]

# Collect followers
result = collect_follower(model, lib, Li, f̃, trans_1, trans_2)

# Print the result
for list in result
    display(list)
end
```

The output of the above code is as follows:

```math
    \left[\begin{smallmatrix}1\end{smallmatrix}\right], \left[\begin{smallmatrix}x^{2} + y^{2}\end{smallmatrix}\right], \left[\begin{smallmatrix}x^{4} + 2 x^{2} y^{2} + y^{4}\end{smallmatrix}\right], \left[\begin{smallmatrix}x^{5} - 10 x^{3} y^{2} + 5 x y^{4}\end{smallmatrix}\right], \cdots
```
