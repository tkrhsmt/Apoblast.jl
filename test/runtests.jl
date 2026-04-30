using Apoblast
using Test
using SymPy

@testset "Apoblast.jl / Model" begin
    model = Apoblast.Model(("x", "y"), ("u", "v"))
    @test model.coords == (Sym("x"), Sym("y"))
    @test model.fields == (SymFunction("u")(Sym("x"), Sym("y")), SymFunction("v")(Sym("x"), Sym("y")))
end

@testset "Apoblast.jl / Library" begin
    model = Model(("x", "y"), ("u", "v"))
    x, y = model.coords
    u, v = model.fields

    @test Apoblast.Param.Library(model, (u, v, u + v, diff(u, x), sin(v))).term ==
          (u, v, u + v, diff(u, x), sin(v))

    z = Sym("z")
    w = SymFunction("w")(x, y)

    @test_throws ArgumentError Apoblast.Param.Library(model, (u + z,))
    @test_throws ArgumentError Apoblast.Param.Library(model, (w,))
end

@testset "Apoblast.jl / Transformation" begin
    model = Model(("x", "y"), ("u", "v"))
    x, y = model.coords
    u, v = model.fields

    reflection = Transformation(model, (x, -y), (u, -v))
    @test reflection.coords_replace == (x, -y)
    @test reflection.fields_replace == (u, -v)
    @test reflection.coords_replace_dict == Dict(x => x, y => -y)
    @test reflection.fields_replace_dict == Dict(u => u, v => -v)
    @test reflection.parameter == ()
    @test reflection.parameter_fixed == Dict()
    @test reflection.jacobian_inverse == Sym[1 0; 0 -1]

    a = Sym("a")
    shifted = Transformation(model, (x + a, -y), (u + a, -v); parameter=((a, Sym(0)),))
    @test shifted.parameter == (a,)
    @test shifted.parameter_fixed == Dict(a => Sym(0))
    @test shifted.coords_replace_dict == Dict(x => x + a, y => -y)
    @test shifted.fields_replace_dict == Dict(u => u + a, v => -v)
    @test shifted.jacobian_inverse == Sym[1 0; 0 -1]
end

@testset "Apoblast.jl / apply_trans" begin
    model = Model(("x", "y"), ("u", "v"))
    x, y = model.coords
    u, v = model.fields

    reflection = Transformation(model, (x, -y), (u, -v))
    @test Apoblast.Core.apply_trans(model, reflection, x) == x
    @test Apoblast.Core.apply_trans(model, reflection, y) == -y
    @test Apoblast.Core.apply_trans(model, reflection, v) == -v
    @test Apoblast.Core.apply_trans(model, reflection, sin(x + v)) == sin(x - v)
    @test Apoblast.Core.apply_trans(model, reflection, diff(u, y)) == -diff(u, y)
    @test Apoblast.Core.apply_trans(model, reflection, diff(u, y, 2)) == diff(u, y, 2)

    shifted = Transformation(model, (x, y), (u + x, v))
    @test Apoblast.Core.apply_trans(model, shifted, diff(u, x)) == diff(u, x) + 1
end

@testset "Apoblast.jl / coeff_matrix" begin
    x = Sym("x")
    y = Sym("y")
    z = Sym("z")

    T, D, g = Apoblast.Core.coeff_matrix(Sym[2*x+3*y, x-z], Sym[x, y])

    @test T == Sym[2 3; 1 0]
    @test D == Sym[0; -1;;]
    @test g == Sym[z]
end

@testset "Apoblast.jl / trans_matrix" begin
    model = Model(("x", "y"), ("u", "v"))
    x, y = model.coords
    u, v = model.fields

    reflection = Transformation(model, (x, -y), (u, -v))
    T, D, g = Apoblast.Core.trans_matrix(model, reflection, Sym[v, diff(u, y)])
    @test T == Sym[-1 0; 0 -1]
    @test size(D) == (2, 0)
    @test g == Sym[]

    shifted = Transformation(model, (x, y), (u + x, v))
    T, D, g = Apoblast.Core.trans_matrix(model, shifted, Sym[diff(u, x)])
    @test T == Sym[1;;]
    @test D == Sym[1;;]
    @test g == Sym[Sym(1)]
end

@testset "Apoblast.jl / infinitesimal_trans_matrix" begin
    model = Model(("x", "y"), ("u", "v"))
    x, y = model.coords
    u, v = model.fields
    θ = Sym("θ")

    rotation = Transformation(
        model,
        (x * cos(θ) + y * sin(θ), y * cos(θ) - x * sin(θ)),
        (u * cos(θ) + v * sin(θ), v * cos(θ) - u * sin(θ));
        parameter=((θ, Sym(0)),),
    )

    T, D, g = Apoblast.Core.infinitesimal_trans_matrix(model, rotation, θ, Sym[u, diff(u, x)])

    @test T == Sym[0 0; 0 0]
    @test sum(D[1, j] * g[j] for j in eachindex(g)) == v
    @test sum(D[2, j] * g[j] for j in eachindex(g)) == diff(v, x) - diff(u, y)
end
