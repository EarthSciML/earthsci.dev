using Documenter

makedocs(;
    authors = "EarthSciML Authors and Contributors",
    repo = "https://github.com/EarthSciML/earthsci.dev/blob/{commit}{path}#{line}",
    sitename = "EarthSciML",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "earthsci.dev",
        assets = ["assets/favicon.ico"],
        repolink = "https://github.com/EarthSciML/earthsci.dev"
    ),
    pages = [
        "Home" => "index.md",
        "Getting Started" => [
            "Overview" => "quickstart/getting_started.md",
            "Preparing Your Environment" => "quickstart/prepare_environment.md",
            "Using ModelingToolkit" => "quickstart/using_mtk.md",
            "Using EarthSciML" => "quickstart/using_earthsciml.md"
        ],
        "Examples" => [
            "Overview" => "examples/overview.md",
            "Scenario Analysis" => "examples/scenario_analysis.md",
            "Optimization" => "examples/optimization.md"
        ],
        "EarthSci EDU" => [
            "Overview" => "edu/overview.md",
            "Atmospheric Chemistry" => "edu/atmos_chem.md"
        ],
        "Contributing" => "contributing.md",
        "Libraries" => [
            "Overview" => "libraries/overview.md",
            "🔗 EarthSciMLBase.jl" => "libraries/base.md",
            "🔗 GasChem.jl" => "libraries/gaschem.md",
            "🔗 Aerosol.jl" => "libraries/aerosol.md",
            "🔗 AtmosphericDeposition.jl" => "libraries/deposition.md",
            "🔗 EnvironmentalTransport.jl" => "libraries/transport.md",
            "🔗 EarthSciData.jl" => "libraries/data.md"
        ]
    ]
)

deploydocs(;
    repo = "github.com/EarthSciML/earthsci.dev",
    devbranch = "main"
)
