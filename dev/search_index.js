var documenterSearchIndex = {"docs":
[{"location":"libraries/data/","page":"🔗 EarthSciData.jl","title":"🔗 EarthSciData.jl","text":"","category":"page"},{"location":"libraries/data/#title:-\"Redirecting...\"","page":"🔗 EarthSciData.jl","title":"title: \"Redirecting...\"","text":"","category":"section"},{"location":"libraries/data/","page":"🔗 EarthSciData.jl","title":"🔗 EarthSciData.jl","text":"<html>\n<head>\n    <meta http-equiv=\"refresh\" content=\"0; url=https://data.earthsci.dev/\" />\n</head>\n<body>\n    <p>If you are not redirected automatically, follow this <a href=\"https://data.earthsci.dev/\">link</a>.</p>\n</body>\n</html>","category":"page"},{"location":"libraries/overview/#Available-Software-Libraries","page":"Overview","title":"Available Software Libraries","text":"","category":"section"},{"location":"contributing/","page":"Contributing","title":"Contributing","text":"+++ title = \"Contributing\" subtitle = \"\" hascode = true date = Date(2023, 5, 5) rss = \"How to contribute\"","category":"page"},{"location":"contributing/","page":"Contributing","title":"Contributing","text":"tags = [\"syntax\", \"code\"] +++","category":"page"},{"location":"contributing/#Contributing-to-EarthSciML","page":"Contributing","title":"Contributing to EarthSciML","text":"","category":"section"},{"location":"contributing/","page":"Contributing","title":"Contributing","text":"See the SciML contributor's guide here.","category":"page"},{"location":"libraries/aerosol/","page":"🔗 Aerosol.jl","title":"🔗 Aerosol.jl","text":"","category":"page"},{"location":"libraries/aerosol/#title:-\"Redirecting...\"","page":"🔗 Aerosol.jl","title":"title: \"Redirecting...\"","text":"","category":"section"},{"location":"libraries/aerosol/","page":"🔗 Aerosol.jl","title":"🔗 Aerosol.jl","text":"<html>\n<head>\n    <meta http-equiv=\"refresh\" content=\"0; url=https://aerosol.earthsci.dev/\" />\n</head>\n<body>\n    <p>If you are not redirected automatically, follow this <a href=\"https://aerosol.earthsci.dev/\">link</a>.</p>\n</body>\n</html>","category":"page"},{"location":"libraries/transport/","page":"🔗 EnvironmentalTransport.jl","title":"🔗 EnvironmentalTransport.jl","text":"","category":"page"},{"location":"libraries/transport/#title:-\"Redirecting...\"","page":"🔗 EnvironmentalTransport.jl","title":"title: \"Redirecting...\"","text":"","category":"section"},{"location":"libraries/transport/","page":"🔗 EnvironmentalTransport.jl","title":"🔗 EnvironmentalTransport.jl","text":"<html>\n<head>\n    <meta http-equiv=\"refresh\" content=\"0; url=https://transport.earthsci.dev/\" />\n</head>\n<body>\n    <p>If you are not redirected automatically, follow this <a href=\"https://transport.earthsci.dev/\">link</a>.</p>\n</body>\n</html>","category":"page"},{"location":"libraries/base/","page":"🔗 EarthSciMLBase.jl","title":"🔗 EarthSciMLBase.jl","text":"","category":"page"},{"location":"libraries/base/#title:-\"Redirecting...\"","page":"🔗 EarthSciMLBase.jl","title":"title: \"Redirecting...\"","text":"","category":"section"},{"location":"libraries/base/","page":"🔗 EarthSciMLBase.jl","title":"🔗 EarthSciMLBase.jl","text":"<html>\n<head>\n    <meta http-equiv=\"refresh\" content=\"0; url=https://base.earthsci.dev/\" />\n</head>\n<body>\n    <p>If you are not redirected automatically, follow this <a href=\"https://base.earthsci.dev/\">link</a>.</p>\n</body>\n</html>","category":"page"},{"location":"libraries/deposition/","page":"🔗 AtmosphericDeposition.jl","title":"🔗 AtmosphericDeposition.jl","text":"","category":"page"},{"location":"libraries/deposition/#title:-\"Redirecting...\"","page":"🔗 AtmosphericDeposition.jl","title":"title: \"Redirecting...\"","text":"","category":"section"},{"location":"libraries/deposition/","page":"🔗 AtmosphericDeposition.jl","title":"🔗 AtmosphericDeposition.jl","text":"<html>\n<head>\n    <meta http-equiv=\"refresh\" content=\"0; url=https://gaschem.earthsci.dev/\" />\n</head>\n<body>\n    <p>If you are not redirected automatically, follow this <a href=\"https://gaschem.earthsci.dev/\">link</a>.</p>\n</body>\n</html>","category":"page"},{"location":"quickstart/3d_sim/#3D-Simulation","page":"3D Simulation","title":"3D Simulation","text":"","category":"section"},{"location":"quickstart/3d_sim/","page":"3D Simulation","title":"3D Simulation","text":"using EarthSciMLBase, EarthSciData, GasChem, EnvironmentalTransport\nusing ModelingToolkit, DifferentialEquations\nusing DomainSets, Unitful, Dates\nusing NCDatasets, Plots\nusing SciMLBase: DiscreteCallback\nusing ProgressMeter\n\n@parameters t [unit = u\"s\", description = \"Time\"]\n@parameters lat = 40\n@parameters lon = -97\n@parameters lev = 1\nemis = NEI2016MonthlyEmis(\"mrggrid_withbeis_withrwc\", t, lon, lat, lev; dtype=Float64)\n\ngeosfp = GEOSFP(\"4x5\", t; dtype = Float64,\n    coord_defaults = Dict(:lon => 0.0, :lat => 0.0, :lev => 1.0))\n\nstarttime = datetime2unix(DateTime(2016, 5, 1, 0, 0))\nendtime = datetime2unix(DateTime(2016, 5, 1, 5, 0))\n\ndomain = DomainInfo(\n    [partialderivatives_δxyδlonlat,\n        partialderivatives_δPδlev_geosfp(geosfp)],\n    constIC(16.0, t ∈ Interval(starttime, endtime)),\n    constBC(16.0,\n        lon ∈ Interval(deg2rad(-130.0), deg2rad(-60.0)),\n        lat ∈ Interval(deg2rad(9.75), deg2rad(60.0)),\n        lev ∈ Interval(1, 15)),\n    dtype = Float64)\n\nchem = SuperFast(t)\nphotolysis = FastJX(t)\n\noutfile = (\"RUNNER_TEMP\" ∈ keys(ENV) ? ENV[\"RUNNER_TEMP\"] : tempname()) * \"out.nc\" # This is just a location to save the output.\noutput = NetCDFOutputter(outfile, 3600.0)\n\ncsys = couple(chem, photolysis, geosfp, emis, domain, output)\n\nadv = AdvectionOperator(600.0, l94_stencil)\n\ncsys = couple(csys, adv)\n\nfunction pbar(start, finish)\n    p = Progress(Int(round(finish-start)))\n    DiscreteCallback(\n        (_, _, _) -> true, \n        (integrator) -> update!(p, Int(round(integrator.t-start)));\n        save_positions = (false, false),\n    )\nend\n#csys = couple(csys, pbar(starttime, endtime))\n\n\nsim = Simulator(csys, [deg2rad(15), deg2rad(10), 1])\nst = SimulatorStrangThreads(Tsit5(), SSPRK22(), 600.0)\n\n@time run!(sim, st, save_on=false, save_start=false, save_end=false,\n    initialize_save=false)\n\nds = NCDataset(outfile, \"r\")\n\nanim = @animate for i ∈ 1:size(ds[\"SuperFast₊O3\"])[4]\n    plot(\n        heatmap(ds[\"SuperFast₊O3\"][:, :, 1, i]', title=\"Ground-Level\"),\n        heatmap(ds[\"SuperFast₊O3\"][:, 2, :, i]', title=\"Vertical Cross-Section\"),\n    )\nend\ngif(anim, fps = 15)","category":"page"},{"location":"quickstart/example/","page":"Example","title":"Example","text":"+++ title = \"EarthSciML Example\" subtitle = \"\" hascode = true date = Date(2023, 7, 1) rss = \"EarthSciML Example\"","category":"page"},{"location":"quickstart/example/","page":"Example","title":"Example","text":"tags = [\"syntax\", \"code\"] +++","category":"page"},{"location":"quickstart/example/#Example","page":"Example","title":"Example","text":"","category":"section"},{"location":"quickstart/example/","page":"Example","title":"Example","text":"\\toc","category":"page"},{"location":"quickstart/example/#Figure-1","page":"Example","title":"Figure 1","text":"","category":"section"},{"location":"quickstart/example/","page":"Example","title":"Example","text":"\nusing EarthSciData, EarthSciMLBase, GasChem,\n    DomainSets, ModelingToolkit, MethodOfLines, \n    DifferentialEquations, Dates, Distributions,\n    Latexify, Plots, SciMLBase,\n    CairoMakie, GraphMakie, MetaGraphsNext\n\n@parameters t lev lon lat\nmodel = SuperFast(t)\nls = latexify(model.rxn_sys) # TODO: Change to model.rxn_sys\nprint(replace(ls.s, raw\"\\require{mhchem}\" => \"\")) # hide","category":"page"},{"location":"quickstart/example/","page":"Example","title":"Example","text":"This is what ls looks like:","category":"page"},{"location":"quickstart/example/","page":"Example","title":"Example","text":"\\textoutput{./code/ex1}","category":"page"},{"location":"quickstart/example/#Figure-2","page":"Example","title":"Figure 2","text":"","category":"section"},{"location":"quickstart/example/","page":"Example","title":"Example","text":"# Set start and end times.\nstart, finish = Dates.datetime2unix.([\n    Dates.DateTime(2022, 5, 1),\n    Dates.DateTime(2022, 5, 3),\n])\n\n# Create a model by combining SuperFast chemistry and FastJX photolysis.\nmodel_ode = SuperFast(t) + FastJX(t)\n\n# Run the simulation.\nprob = ODEProblem(structural_simplify(get_mtk(model_ode)), [], (start, finish), [])\nsol = solve(prob, TRBDF2(), saveat=1800.0)\n\n# Plot result.\nticks = Dates.datetime2unix.([Dates.DateTime(2022, 5, 1), Dates.DateTime(2022, 5, 2), \n    Dates.DateTime(2022, 5, 3)])\ntickstrings = Dates.format.(Dates.unix2datetime.(ticks), \"m-d-yy\")\nPlots.plot(sol,ylims=(0,30),xlabel=\"Date\", ylabel=\"Concentration (ppb)\", \n    ylabelfontsize=9, xlabelfontsize=9, \n    xticks=(ticks, tickstrings), legend=:outertopright, size=(500, 310))\nPlots.savefig(joinpath(@OUTPUT, \"ode.svg\")) # hide","category":"page"},{"location":"quickstart/example/","page":"Example","title":"Example","text":"\\textoutput{./code/ex2}","category":"page"},{"location":"quickstart/example/","page":"Example","title":"Example","text":"\\fig{ode}","category":"page"},{"location":"quickstart/example/#Figure-3","page":"Example","title":"Figure 3","text":"","category":"section"},{"location":"quickstart/example/","page":"Example","title":"Example","text":"struct Emissions <: EarthSciMLODESystem\n    sys\n    function Emissions(t, μ_lon, σ)\n        @variables NO(t) ISOP(t)\n        dist = MvNormal([start,μ_lon],[3600.0,σ])\n        @parameters lon\n        D = Differential(t)\n        new(ODESystem([\n            D(NO) ~ pdf(dist, [t, lon]) * 50, \n            D(ISOP) ~ pdf(dist, [t, lon]) * 50,\n        ], t, name=:emissions))\n    end\nend\nBase.:(+)(e::Emissions, b::SuperFast) = operator_compose(b, e)\nBase.:(+)(b::SuperFast, e::Emissions) = e + b","category":"page"},{"location":"quickstart/example/","page":"Example","title":"Example","text":"\\textoutput{./code/ex3}","category":"page"},{"location":"quickstart/example/#Supplemental-code-to-skip-unit-enforcement","page":"Example","title":"Supplemental code to skip unit enforcement","text":"","category":"section"},{"location":"quickstart/example/","page":"Example","title":"Example","text":"function ModelingToolkit.check_units(eqs...) # Skip unit enforcement for now\n    nothing\n    #ModelingToolkit.validate(eqs...) || @info \"Some equations had invalid units. See warnings for details.\"\nend\n\n# Skip returning the observed variables (i.e. variables that are not state variables)\n# because it is currently extremely slow to do so for this demo. \nSciMLBase.observed(sol::SciMLBase.AbstractTimeseriesSolution, sym, i::Colon) = zeros(Float64, length(sol.t))","category":"page"},{"location":"quickstart/example/","page":"Example","title":"Example","text":"\\textoutput{./code/ex4}","category":"page"},{"location":"quickstart/example/#Figure-4","page":"Example","title":"Figure 4","text":"","category":"section"},{"location":"quickstart/example/","page":"Example","title":"Example","text":"domain = DomainInfo(\n    partialderivatives_lonlat2xymeters,\n    constIC(1.0f0, t ∈ Interval(start, finish)),\n    zerogradBC(lon ∈ Interval(-130f0, -100f0)))\n\ngeos = GEOSFP(\"0.25x0.3125_NA\", t; coord_defaults = Dict(:lat => 34.0, :lev => 1))\n\nmodel = SuperFast(t) + FastJX(t) + domain +\n    Emissions(t, -118.2, 0.2)+ Advection() + geos\n\ng = graph(model)\n\nf, ax, p = graphplot(g; ilabels=labels(g))\nhidedecorations!(ax); hidespines!(ax); ax.aspect = DataAspect()\n\nsave(joinpath(@OUTPUT, \"graph.svg\"), f) # hide","category":"page"},{"location":"quickstart/example/","page":"Example","title":"Example","text":"\\output{./code/ex5a} \\fig{graph}","category":"page"},{"location":"quickstart/example/","page":"Example","title":"Example","text":"discretization = MOLFiniteDifference([lon => 50], t, approx_order=2)\nprob = discretize(get_mtk(model), discretization)\nsol = solve(prob, TRBDF2(), saveat=3600.0)","category":"page"},{"location":"quickstart/example/","page":"Example","title":"Example","text":"\\output{./code/ex5b}","category":"page"},{"location":"quickstart/example/#Supplemental-plotting-code","page":"Example","title":"Supplemental plotting code","text":"","category":"section"},{"location":"quickstart/example/","page":"Example","title":"Example","text":"discrete_lon = sol[lon]\ndiscrete_t = sol[t]\n\n@variables superfast₊O3(..) superfast₊NO(..) superfast₊ISOP(..) superfast₊NO2(..)\nsol_isop = sol[superfast₊ISOP(t, lon)]\nsol_o3 = sol[superfast₊O3(t, lon)]\nsol_no = sol[superfast₊NO(t, lon)]\nsol_no2 = sol[superfast₊NO2(t, lon)]\n\nusing LaTeXStrings\nPlots.plot(\n    Plots.heatmap(discrete_t, discrete_lon[1:end], sol_isop[:, 1:end]', \n        xticks=(ticks, tickstrings), \n        ylabel=\"Longitude (°)\", title=\"Isoprene\", titlefontsize=12,\n        margin=0Plots.mm,\n        size=(600, 400), dpi=300),\n    Plots.heatmap(discrete_t, discrete_lon[1:end], sol_no[:, 1:end]', \n        xticks=(ticks, tickstrings), title=\"NO\", titlefontsize=12),\n    Plots.heatmap(discrete_t, discrete_lon[1:end], sol_no2[:, 1:end]', \n        xticks=(ticks, tickstrings), xlabel=\"Date\", title=L\"\\textrm{NO_2}\"),\n    Plots.heatmap(discrete_t, discrete_lon[1:end], sol_o3[:, 1:end]',colorbar_title=\"Concentration (ppb)\",\n        xticks=(ticks, tickstrings), title=L\"\\textrm{O_3}\", titlefontsize=12),\n)\nPlots.savefig(joinpath(@OUTPUT, \"pde.svg\"))","category":"page"},{"location":"quickstart/example/","page":"Example","title":"Example","text":"\\textoutput{./code/ex6}","category":"page"},{"location":"quickstart/example/","page":"Example","title":"Example","text":"\\fig{pde}","category":"page"},{"location":"libraries/gaschem/","page":"🔗 GasChem.jl","title":"🔗 GasChem.jl","text":"","category":"page"},{"location":"libraries/gaschem/#title:-\"Redirecting...\"","page":"🔗 GasChem.jl","title":"title: \"Redirecting...\"","text":"","category":"section"},{"location":"libraries/gaschem/","page":"🔗 GasChem.jl","title":"🔗 GasChem.jl","text":"<html>\n<head>\n    <meta http-equiv=\"refresh\" content=\"0; url=https://gaschem.earthsci.dev/\" />\n</head>\n<body>\n    <p>If you are not redirected automatically, follow this <a href=\"https://gaschem.earthsci.dev/\">link</a>.</p>\n</body>\n</html>","category":"page"},{"location":"#Earth-Science-Machine-Learning","page":"Home","title":"Earth Science Machine Learning","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A next-generation software framework for geoscientific modeling and analysis","category":"page"},{"location":"#Available-libraries","page":"Home","title":"Available libraries","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"EarthSciMLBase.jl\nEarthSciData.jl\nAerosol.jl\nAtmosphericDeposition.jl\nGasChem.jl\nEnvironmentalTransport.jl","category":"page"},{"location":"#Reproducibility","page":"Home","title":"Reproducibility","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"<details><summary>The documentation of this EarthSciML package was built using these direct dependencies,</summary>","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg # hide\nPkg.status() # hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"</details>","category":"page"},{"location":"","page":"Home","title":"Home","text":"<details><summary>and using this machine and Julia version.</summary>","category":"page"},{"location":"","page":"Home","title":"Home","text":"using InteractiveUtils # hide\nversioninfo() # hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"</details>","category":"page"},{"location":"","page":"Home","title":"Home","text":"<details><summary>A more complete overview of all dependencies and their versions is also provided.</summary>","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg # hide\nPkg.status(;mode = PKGMODE_MANIFEST) # hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"</details>","category":"page"},{"location":"","page":"Home","title":"Home","text":"You can also download the \n<a href=\"","category":"page"},{"location":"","page":"Home","title":"Home","text":"using TOML\nusing Markdown\nversion = TOML.parse(read(\"../Project.toml\",String))[\"version\"]\nname = TOML.parse(read(\"../Project.toml\",String))[\"name\"]\nlink = Markdown.MD(\"https://github.com/EarthSciML/\"*name*\"/tree/gh-pages/v\"*version*\"/assets/Manifest.toml\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"\">manifest</a> file and the\n<a href=\"","category":"page"},{"location":"","page":"Home","title":"Home","text":"using TOML\nusing Markdown\nversion = TOML.parse(read(\"../Project.toml\",String))[\"version\"]\nname = TOML.parse(read(\"../Project.toml\",String))[\"name\"]\nlink = Markdown.MD(\"https://github.com/EarthSciML/\"*name*\"/tree/gh-pages/v\"*version*\"/assets/Project.toml\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"\">project</a> file.","category":"page"},{"location":"quickstart/getting_started/#Getting-started","page":"Installation","title":"Getting started","text":"","category":"section"},{"location":"quickstart/getting_started/","page":"Installation","title":"Installation","text":"First, install the package:","category":"page"},{"location":"quickstart/getting_started/","page":"Installation","title":"Installation","text":"] add EarthSciMLBase","category":"page"},{"location":"quickstart/getting_started/","page":"Installation","title":"Installation","text":"Then, load it:","category":"page"},{"location":"quickstart/getting_started/","page":"Installation","title":"Installation","text":"> using EarthSciMLBase","category":"page"}]
}
