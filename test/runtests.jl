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
