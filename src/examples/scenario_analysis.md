# Scenario Analysis

A common application in geoscientific modeling is conducting scenario analysis: comparing the results of a model under different sets of assumptions or conditions. In this example, we will demonstrate how to use EarthSciML to conduct a simple scenario analysis, comparing the results of a model of atmospheric chemistry under two different sets of emissions.

First, we create a base model containing all the components that will be shared between the two scenarios. This model is similar to the one described in detail in [Using EarthSciML](@ref), including a gas-phase chemistry model, a photolysis model, emissions data, and meteorological data:

```@example scenario_analysis
using EarthSciMLBase, GasChem, AtmosphericDeposition, EarthSciData
using EnvironmentalTransport, ModelingToolkit, OrdinaryDiffEq
using Dates, Plots, NCDatasets
using ProgressLogging # Needed for progress bar. Use `TerminalLoggers` if in a terminal.

domain = DomainInfo(
    DateTime(2016, 5, 1),
    DateTime(2016, 5, 3);
    lonrange = deg2rad(-115):deg2rad(2.5):deg2rad(-68.75),
    latrange = deg2rad(25):deg2rad(2):deg2rad(53.7),
    levrange = 1:15,
    dtype = Float32)

geosfp = GEOSFP("0.5x0.625_NA", domain; stream=false)
geosfp = EarthSciMLBase.copy_with_change(geosfp, discrete_events=[]) # Workaround for bug.

emis = NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", domain; stream=false)
emis = EarthSciMLBase.copy_with_change(emis, discrete_events=[]) # Workaround for bug.

dt = 300.0 # Operator splitting timestep

base_model = couple(
    SuperFast(),
    FastJX(),
    #DrydepositionG(), Not currently working
    Wetdeposition(),
    AdvectionOperator(dt, upwind1_stencil, ZeroGradBC()),
    emis,
    geosfp,
    domain
)
```

Once we've created a base model, we can then add in the components that will be unique to each scenario. In this case, we will create two scenarios: a "business as usual" scenario and a "reduced emissions" scenario. For our business as usual scenario, the only component we need to add is one that writes the model output to a file called `bau_output.nc`.

```@example scenario_analysis
bau_outfile = ("RUNNER_TEMP" ∈ keys(ENV) ? ENV["RUNNER_TEMP"] : tempname()) * "bau_output.nc" # This is just a location to save the output.
bau_model = couple(base_model, NetCDFOutputter(bau_outfile, 3600.0))
```

In this case, we want to explore the impacts of a 30% reduction in emissions from the "onroad" sector.
To do this, we initialize a new emissions component only containing the emissions from the onroad sector, and we specify a scale factor of -0.3, which when added to the total emissions from all sectors that are already present in the base model represents the total emissions minus 30% of the onroad emissions.
(Additional information about the emissions data we're using is available [here](https://data.earthsci.dev/dev/api/#EarthSciData.NEI2016MonthlyEmis-Tuple{AbstractString,%20EarthSciMLBase.DomainInfo}).)
We also create a component that writes the model output to a file called `scenario_output.nc`, and then we couple both of these components to the base model to create the scenario model.

```@example scenario_analysis
@named scenario_emis = NEI2016MonthlyEmis("onroad", domain; scale=-0.3, stream=false)
scenario_emis = EarthSciMLBase.copy_with_change(scenario_emis, discrete_events=[]) # Workaround for bug.


scenario_outfile = ("RUNNER_TEMP" ∈ keys(ENV) ? ENV["RUNNER_TEMP"] : tempname()) * "scenario_output.nc"
scenario_model = couple(base_model, scenario_emis, NetCDFOutputter(scenario_outfile, 3600.0))
```

Once we've configured both of out models, we can run them both as described in [Using EarthSciML](@ref) above:

```@example scenario_analysis
st = SolverStrangSerial(Rosenbrock23(), dt)

bau_prob = ODEProblem(bau_model, st)
sol = solve(bau_prob, SSPRK22(); dt=dt, progress=true, progress_steps=1,
    save_on=false, save_start=false, save_end=false, initialize_save=false)


scenario_prob = ODEProblem(scenario_model, st)
sol = solve(scenario_prob, SSPRK22(); dt=dt, progress=true, progress_steps=1,
    save_on=false, save_start=false, save_end=false, initialize_save=false)
nothing # hide
```

One the simulations are complete, we can compare the results of the two scenarios. 
In this case, we'll create an animation that shows the difference in ground-level ozone concentrations between the two scenarios at each time step:

```@example scenario_analysis
bau_data = NCDataset(bau_outfile, "r");
scenario_data = NCDataset(scenario_outfile, "r");

anim = @animate for i ∈ 1:size(bau_data["SuperFast₊O3"])[4]
    Plots.plot(
        Plots.heatmap(bau_data["SuperFast₊O3"][:, :, 1, i]', title="Business as Usual"),
        Plots.heatmap(scenario_data["SuperFast₊O3"][:, :, 1, i]', title="Scenario"),
        Plots.heatmap(scenario_data["SuperFast₊O3"][:, :, 1, i]' - 
            bau_data["SuperFast₊O3"][:, :, 1, i]', title="Difference"),
        size=(1900, 400), layout=(1, 3)
    )
end
gif(anim, fps = 5)
```

You can see that in this case, the scenario with reduced emissions actually results in an increase on ozone concentrations during the short simulation time that we run.
Can you think of any way to use the models above to explore why this might be happening?