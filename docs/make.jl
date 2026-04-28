using Apoblast
using Documenter

DocMeta.setdocmeta!(Apoblast, :DocTestSetup, :(using Apoblast); recursive=true)

makedocs(;
    modules=[Apoblast],
    authors="T. Hashimoto",
    sitename="Apoblast.jl",
    format=Documenter.HTML(;
        canonical="https://tkrhsmt.github.io/Apoblast.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/tkrhsmt/Apoblast.jl",
    devbranch="main",
)
