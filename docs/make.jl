using Documenter, DiseaseOutbreak

makedocs(;
    modules=[DiseaseOutbreak],
    format=Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages=[
        "Home" => "index.md",
        "Model" => "models.md"
    ],
    repo="https://github.com/chelseatrotter/DiseaseOutbreak.jl/blob/{commit}{path}#L{line}",
    sitename="DiseaseOutbreak.jl",
    authors="Sen Research Group, Saunak Sen, Hyeonju Kim, Gregory Farage, Chelsea Trotter",
    assets=String[],
)

deploydocs(;
    repo="github.com/chelseatrotter/DiseaseOutbreak.jl",
)
