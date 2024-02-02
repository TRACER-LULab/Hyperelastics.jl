using Hyperelastics
using Documenter
# using DocumenterCitations

# bib = CitationBibliography(joinpath(@__DIR__, "src", "paper.bib"); style = :numeric)

DocMeta.setdocmeta!(Hyperelastics, :DocTestSetup, :(using Hyperelastics); recursive = true)

makedocs(;
    # plugins = [bib],
    modules = [Hyperelastics],
    authors = "Carson Farmer <59753859+cfarm6@users.noreply.github.com> and contributors",
    repo = "https://github.com/TRACER-LULab/Hyperelastics.jl/blob/{commit}{path}#{line}",
    sitename = "Hyperelastics.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://TRACER-LULab.github.io/Hyperelastics.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = ["Home" => "index.md", "API" => "API.md", "Example" => "example.md"],
)

deploydocs(; repo = "github.com/TRACER-LULab/Hyperelastics.jl", devbranch = "main")
