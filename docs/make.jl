using Hyperelastics
using Documenter

DocMeta.setdocmeta!(Hyperelastics, :DocTestSetup, :(using Hyperelastics); recursive=true)

makedocs(;
    modules=[Hyperelastics],
    authors="Carson Farmer",
    repo="https://github.com/cfarm6/Hyperelastics.jl/blob/{commit}{path}#{line}",
    sitename="Hyperelastics.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://TRACER-LULab.github.io/Hyperelastics.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/cfarm6/Hyperelastics.jl",
    devbranch="main",
)
