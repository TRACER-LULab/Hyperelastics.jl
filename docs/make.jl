using Hyperelastics
using Documenter
using DocumenterCitations

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "paper.bib");
    style=:numeric
)

DocMeta.setdocmeta!(Hyperelastics, :DocTestSetup, :(using Hyperelastics); recursive=true)

makedocs(bib;
    modules=[Hyperelastics],
    authors="Carson Farmer <59753859+cfarm6@users.noreply.github.com> and contributors",
    repo="https://github.com/cfarm6/Hyperelastics.jl/blob/{commit}{path}#{line}",
    sitename="Hyperelastics.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://cfarm6.github.io/Hyperelastics.jl",
        edit_link="main",
        assets=String[]
    ),
    pages=[
        "Home" => "index.md",
        "API" => "API.md",
        "Example" => "example.md",
    ]
)

deploydocs(;
    repo="github.com/cfarm6/Hyperelastics.jl",
    devbranch="main"
)
