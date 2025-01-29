# Using EarthSciML

As we saw in [Using ModelingToolkit](@ref), keeping track of variables and equations when working with complex models can be a challenge.
EarthSciML is designed to streamline this process in the case of geoscientific model by providing utilities for coupling model components together and running them, and also by providing reference implementations of standard model components that can be used as building blocks to create larger-scale geoscientific models.

## Model Components

EarthSciML model components are ModelingToolkit `ODESystems` which contain metadata which specifies how each component should be coupled to other components.

So, for example, we can create an instance of the [SuperFast](https://gaschem.earthsci.dev/stable/superfast/) gas-phase chemistry model model from the [GasChem](https://gaschem.earthsci.dev/stable/) package like this:

```@example using_earthsciml
using GasChem

chem = SuperFast()
```

It's just an ODESystem, so we can run it like any other ModelingToolkit model:

```@example using_earthsciml
using ModelingToolkit, OrdinaryDiffEq, Plots

sol = solve(ODEProblem(structural_simplify(chem)), tspan=(0,2*3600))

plot(sol, legend=:topright, xlabel="Time (s)", ylabel="Concentration (ppb)",
    title="SuperFast Chemistry")
```

## Coupling Model Components

EarthSciML model components differ from ordinary ModelingToolkit models in that they contain metadata that specifies how they should be coupled to other components.
(In particular, this is done by specifying `:coupletype` metadata and the `EarthSciMLBase.couple2` function, which is covered in more detail [here](https://base.earthsci.dev/dev/composition/).)
This allows us to use the `EarthSciMLBase.couple` function to couple model components together into a single model that can be run as a single unit, for example as shown below:

```@example using_earthsciml
using EarthSciMLBase

chem_phot = couple(
    SuperFast(),
    FastJX()
)
```

The above example couples the SuperFast with the [Fast-JX](https://gaschem.earthsci.dev/stable/api/#GasChem.FastJX-Tuple{}) photolysis model, and when we run them together as a coupled model, we can see a characteristic diurnal cycle in the O₃ concentration, were concentrations are lower at night and higher during the day owing to photochemistry:

```@example using_earthsciml
chem_phot_sys = convert(ODESystem, chem_phot)
sol = solve(ODEProblem(chem_phot_sys), tspan=(0,3*24*3600))

plot(sol.t, sol[chem_phot_sys.SuperFast₊O3], legend=:topright, xlabel="Time (s)", 
    ylabel="O₃ Concentration (ppb)", label=:none, 
    title="SuperFast Chemistry + Fast-JX Photolysis")
```

# External Data

We can also bring in external data to provide forcing to our model using [EarthSciData](https://data.earthsci.dev).
For example, below we load [GEOS-FP](https://data.earthsci.dev/stable/geosfp/) meteorological data and pollutant emissions from the U.S. [National Emissions Inventory](https://data.earthsci.dev/stable/nei2016/).

```@example using_earthsciml
using EarthSciData
using Dates

domain = DomainInfo(
    DateTime(2016, 5, 1),
    DateTime(2016, 5, 3);
    lonrange = deg2rad(-115):deg2rad(2.5):deg2rad(-68.75),
    latrange = deg2rad(25):deg2rad(2):deg2rad(53.7),
    levrange = 1:15,
    dtype = Float32)

emis = NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", domain; stream=false)
emis = EarthSciMLBase.copy_with_change(emis, discrete_events=[]) # Workaround for bug.

equations(emis)[1:5]
```

```@example using_earthsciml
geosfp = GEOSFP("0.5x0.625_NA", domain; stream=false)
geosfp = EarthSciMLBase.copy_with_change(geosfp, discrete_events=[]) # Workaround for bug.

equations(geosfp)[1:5]
```

In order to load data products, we will usually need to specify a spatiotemporal domain as demonstrated above.
Simulations that don't include a spatial component will use data from the centerpoint of the provided spatial domain.

We can add our data providers to our model from above.
We'll also add in other processes important to a simplified model of air quality, such as [deposition](https://deposition.earthsci.dev/stable/).
The coupling process will automatically make the appropriate connections between the model components to build a unified model.

```@example using_earthsciml
using AtmosphericDeposition

model = couple(
    SuperFast(),
    FastJX(),
    #DrydepositionG(), Not currently working
    Wetdeposition(),
    emis,
    geosfp,
)
```
## Box Model Simulation

Once we specify our model, we can try running it in a zero-dimensional "box model" format. 
This can allow us study the dynamics or diagnose any programs without using a lot of computing power.

```@example using_earthsciml
model_sys = convert(ODESystem, model)
sol = solve(ODEProblem(model_sys), Rosenbrock23(), tspan=get_tspan(domain))

plot(unix2datetime.(sol.t), sol[model_sys.SuperFast₊O3],
    ylabel="O₃ Concentration (ppb)", 
    xlabel="Time", label=:none, 
    title="Air Quality Box Model")
```

## 3D Model Simulation

Once we're satisfied with our box model, we can move on to a full 3D simulation.
To do this, we'll add in any spatial operators that we'll need, for example [advection](https://transport.earthsci.dev/stable/advection/) from the [EnvironmentalTransport](https://transport.earthsci.dev/stable/) package.
We can also add in a model component that will [write the model results to the disk](https://data.earthsci.dev/dev/api/#EarthSciData.NetCDFOutputter),
which can be useful for large simulations where the results may not all fit in the memory of your computer.
Finally, we also need to add our spatiotemporal domain to the model.

```@example using_earthsciml
using EnvironmentalTransport

dt = 300.0 # Operator splitting timestep
advection = AdvectionOperator(dt, upwind1_stencil, ZeroGradBC())

outfile = ("RUNNER_TEMP" ∈ keys(ENV) ? ENV["RUNNER_TEMP"] : tempname()) * "out.nc" # This is just a location to save the output.
output = NetCDFOutputter(outfile, 3600.0) # Save output every 3600 sec. of simulation time,

model_3d = couple(model, advection, output, domain)
```

### Graph Visualization

Once we've created our model, we can also check a graph representation of its components.
This shows each of the components in the model as nodes and couplings between the components as edges.

```@example using_earthsciml
using MetaGraphsNext
using CairoMakie, GraphMakie

g = graph(model_3d)
f = Figure(backgroundcolor = :transparent)
ax = Axis(f[1, 1], backgroundcolor = :transparent)
p = graphplot!(ax, g; ilabels=labels(g), edge_color=:gray)
CairoMakie.xlims!(ax, -2.4, 1.5)
CairoMakie.ylims!(ax, -2.3, 2.5)
hidedecorations!(ax); hidespines!(ax)
f
```

So for example in the figure above we can see that `WetDeposition` is coupled to SuperFast (because it is depositing the SuperFast chemical species) and to GEOS-FP (because it gets information such as Temperature and cloud fraction from GEOS-FP).

### Running 3D simulation

Running a 3D simulation is similar to running the box model simulation, except we need to specify a [SolverStrategy](https://base.earthsci.dev/stable/api/#EarthSciMLBase.SolverStrategy) which will tell the solver how to combine the non-spatial and spatial components of the model.
Here we will use [SolverStrangSerial](https://base.earthsci.dev/stable/api/#EarthSciMLBase.SolverStrangSerial).

We then use our solver strategy as an argument to `ODEProblem`.
We also specify some options in `solve` show a progress bar (because this simulation can take a while) and to prevent the output data from being saved in memory because we have already specified that we want to save it to disk instead.

```@example using_earthsciml
st = SolverStrangSerial(Rosenbrock23(), dt)

prob = ODEProblem(model_3d, st)

using ProgressLogging # Needed for progress bar. Use `TerminalLoggers` if in a terminal.

sol = solve(prob, SSPRK22(); dt=dt, progress=true, progress_steps=1,
    save_on=false, save_start=false, save_end=false, initialize_save=false)
```

### Visualizing 3D model results

Because we have saved our results to disk instead of in the `sol` variable, we need to visualize differently than we have above. To do so, we can use the [NCDatasets](https://juliageo.org/NCDatasets.jl/stable/) package.

```@example using_earthsciml
using NCDatasets

ds = NCDataset(outfile, "r");

cross_section_lat = size(ds["SuperFast₊O3"], 2) ÷ 2
anim = @animate for i ∈ 1:size(ds["SuperFast₊O3"])[4]
    Plots.plot(
        Plots.heatmap(ds["SuperFast₊O3"][:, :, 1, i]', title="Ground-Level"),
        Plots.heatmap(ds["SuperFast₊O3"][:, cross_section_lat, :, i]', 
            title="Vertical Cross-Section"),
        size=(1200, 400)
    )
end
gif(anim, fps = 5)
```

## Next Steps

Now that you have a basic understanding of how to use EarthSciML, you can take a look at the model components and features available in the different [Available Software Libraries](@ref) that we provide.