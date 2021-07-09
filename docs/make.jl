using Crystallography
using CrystallographyBase
using Documenter

DocMeta.setdocmeta!(
    Crystallography,
    :DocTestSetup,
    :(using Crystallography);
    recursive = true,
)

makedocs(;
    modules = [Crystallography, CrystallographyBase],
    authors = "Qi Zhang <singularitti@outlook.com>",
    repo = "https://github.com/MineralsCloud/Crystallography.jl/blob/{commit}{path}#{line}",
    sitename = "Crystallography.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://mineralscloud.github.io/Crystallography.jl",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => Any["Installation"=>"installation.md", "Development"=>"develop.md"],
        "API by modules" => Any["`Crystallography` module"=>"api.md",],
    ],
)

deploydocs(; repo = "github.com/MineralsCloud/Crystallography.jl")
