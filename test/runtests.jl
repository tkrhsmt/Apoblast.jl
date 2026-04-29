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
    @test reflection.jacobian_inverse == Sym[1 0; 0 -1]

    a = Sym("a")
    shifted = Transformation(model, (x + a, -y), (u + a, -v); parameter = (a,))
    @test shifted.parameter == (a,)
    @test shifted.coords_replace_dict == Dict(x => x + a, y => -y)
    @test shifted.fields_replace_dict == Dict(u => u + a, v => -v)
    @test shifted.jacobian_inverse == Sym[1 0; 0 -1]
end
