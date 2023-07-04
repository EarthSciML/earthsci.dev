+++
title = "EarthSciML Example"
subtitle = ""
hascode = true
date = Date(2023, 7, 1)
rss = "EarthSciML Example"

tags = ["syntax", "code"]
+++

# Example

\toc

## Figure 1
```julia:./code/ex1

using EarthSciData, EarthSciMLBase, GasChem,
    DomainSets, ModelingToolkit, MethodOfLines, 
    DifferentialEquations, Dates, Distributions,
    Latexify, Plots, SciMLBase

@parameters t lev lon lat
model = SuperFast(t)
ls = latexify(model.sys) # TODO: Change to model.rxn_sys
print(ls.s) # hide
```

This is what `ls` looks like:

\textoutput{./code/ex1}

## Figure 2
```julia:./code/ex2
# Set start and end times.
start, finish = Dates.datetime2unix.([
    Dates.DateTime(2022, 5, 1),
    Dates.DateTime(2022, 5, 3),
])

# Create a model by combining SuperFast chemistry and FastJX photolysis.
model_ode = SuperFast(t) + FastJX(t)

# Run the simulation.
prob = ODEProblem(structural_simplify(get_mtk(model_ode)), [], (start, finish), [])
sol = solve(prob, TRBDF2(), saveat=1800.0)

# Plot result.
ticks = Dates.datetime2unix.([Dates.DateTime(2022, 5, 1), Dates.DateTime(2022, 5, 2), Dates.DateTime(2022, 5, 3)])
tickstrings = Dates.format.(Dates.unix2datetime.(ticks), "m-d-yy")
plot(sol,ylims=(0,15),xlabel="Date", ylabel="Concentration (ppb)", ylabelfontsize=9, xlabelfontsize=9, 
    xticks=(ticks, tickstrings), legend=:outertopright, size=(500, 310))
savefig(joinpath(@OUTPUT, "ode.svg")) # hide
```

\fig{ode}

## Figure 3

```julia:./code/ex3
struct Emissions <: EarthSciMLODESystem
    sys
    function Emissions(t, μ_lon, σ)
        @variables NO(t) ISOP(t)
        dist = MvNormal([start,μ_lon],[3600.0,σ])
        @parameters lon
        D = Differential(t)
        new(ODESystem([
            D(NO) ~ pdf(dist, [t, lon]) * 50, 
            D(ISOP) ~ pdf(dist, [t, lon]) * 50,
        ], t, name=:emissions))
    end
end
Base.:(+)(e::Emissions, b::SuperFast) = operator_compose(b, e)
Base.:(+)(b::SuperFast, e::Emissions) = e + b
```

## Supplemental code to skip unit enforcement
```julia:./code/ex4
function ModelingToolkit.check_units(eqs...) # Skip unit enforcement for now
    ModelingToolkit.validate(eqs...) || @info "Some equations had invalid units. See warnings for details."
end

# Skip returning the observed variables (i.e. variables that are not state variables)
# because it is currently extremely slow to do so for this demo. 
SciMLBase.observed(sol::SciMLBase.AbstractTimeseriesSolution, sym, i::Colon) = zeros(Float64, length(sol.t))
```
 
## Figure 4
```julia:./code/ex5
domain = DomainInfo(
    partialderivatives_lonlat2xymeters,
    constIC(1.0f0, t ∈ Interval(start, finish)),
    zerogradBC(lon ∈ Interval(-130f0, -100f0)))

geos = GEOSFP("0.25x0.3125_NA", t; coord_defaults = Dict(:lat => 34.0, :lev => 1))

model = SuperFast(t) + FastJX(t) + domain +
    Emissions(t, -118.2, 0.2)+ Advection() + geos

discretization = MOLFiniteDifference([lon => 50], t, approx_order=2)
prob = discretize(get_mtk(model), discretization)
sol = solve(prob, TRBDF2(), saveat=3600.0)
```

## Supplemental plotting code
```julia:./code/ex6
discrete_lon = sol[lon]
discrete_t = sol[t]

@variables superfast₊O3(..) superfast₊NO(..) superfast₊ISOP(..) superfast₊NO2(..)
sol_isop = sol[superfast₊ISOP(t, lon)]
sol_o3 = sol[superfast₊O3(t, lon)]
sol_no = sol[superfast₊NO(t, lon)]
sol_no2 = sol[superfast₊NO2(t, lon)]

using LaTeXStrings
plot(
    heatmap(discrete_t, discrete_lon[1:end], sol_isop[:, 1:end]', 
        xticks=(ticks, tickstrings), 
        ylabel="Longitude (°)", title="Isoprene", titlefontsize=12,
        margin=0Plots.mm,
        size=(600, 400), dpi=300),
    heatmap(discrete_t, discrete_lon[1:end], sol_no[:, 1:end]', 
        xticks=(ticks, tickstrings), title="NO", titlefontsize=12),
    heatmap(discrete_t, discrete_lon[1:end], sol_no2[:, 1:end]', 
        xticks=(ticks, tickstrings), xlabel="Date", title=L"\textrm{NO_2}"),
    heatmap(discrete_t, discrete_lon[1:end], sol_o3[:, 1:end]',colorbar_title="Concentration (ppb)",
        xticks=(ticks, tickstrings), title=L"\textrm{O_3}", titlefontsize=12),
)
savefig(joinpath(@OUTPUT, "pde.svg"))
```