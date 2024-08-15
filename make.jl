using Documenter

makedocs(;
    authors="EarthSciML Authors and Contributors",
    repo="https://github.com/EarthSciML/earthsci.dev/blob/{commit}{path}#{line}",
    sitename="EarthSciML",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="earthsci.dev",
        assets=String[],
        repolink="https://github.com/EarthSciML/earthsci.dev",
    ),
    pages=[
        "Home" => "index.md",
        "Quick Start" => [
            "Installation" => "quickstart/getting_started.md",
#            "1D Simulation" => "example.md",
            "3D Simulation" => "quickstart/3d_sim.md",
        ],
        "Libraries" => [
            "Overview" => "libraries/overview.md",
            "🔗 EarthSciMLBase.jl" => "libraries/base.md",
            "🔗 GasChem.jl" => "libraries/gaschem.md",
            "🔗 Aerosol.jl" => "libraries/aersol.md",
            "🔗 AtmosphericDeposition.jl" => "libraries/deposition.md",
            "🔗 EnvironmentalTransport.jl" => "libraries/transport.md",
            "🔗 EarthSciData.jl" => "libraries/data.md",
        ],
        "Contributing" => "contributing.md",
    ],
)

deploydocs(;
    repo="github.com/EarthSciML/earthsci.dev",
    devbranch="main",
)
