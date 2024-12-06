# 3D Simulation

```@example 3d_sim
using EarthSciMLBase, EarthSciData, GasChem, EnvironmentalTransport
using AtmosphericDeposition
using ModelingToolkit
using OrdinaryDiffEq
using Dates
using NCDatasets, Plots
using MetaGraphsNext
using CairoMakie, GraphMakie
using ProgressLogging # Needed for progress bar. Use `TerminalLoggers` if in a terminal.

domain = DomainInfo(
    DateTime(2016, 5, 1),
    DateTime(2016, 5, 1, 6);
    lonrange = deg2rad(-115):deg2rad(2.5):deg2rad(-68.75),
    latrange = deg2rad(25):deg2rad(2):deg2rad(53.7),
    levrange = 1:15,
    dtype = Float64)

emis = NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", domain; stream=false)
emis = EarthSciMLBase.copy_with_change(emis, discrete_events=[]) # Workaround for bug.

geosfp = GEOSFP("0.5x0.625_NA", domain; stream=false)
geosfp = EarthSciMLBase.copy_with_change(geosfp, discrete_events=[]) # Workaround for bug.

domain = EarthSciMLBase.add_partial_derivative_func(domain, partialderivatives_δPδlev_geosfp(geosfp))

chem = SuperFast()
photolysis = FastJX()
#drydep = DrydepositionG()
wetdep = Wetdeposition()

dt = 300.0 # Splitting timestep
adv = AdvectionOperator(dt, upwind1_stencil, ZeroGradBC())

outfile = ("RUNNER_TEMP" ∈ keys(ENV) ? ENV["RUNNER_TEMP"] : tempname()) * "out.nc" # This is just a location to save the output.
output = NetCDFOutputter(outfile, 3600.0)

csys = couple(chem, photolysis, geosfp, emis, domain, output, adv,  wetdep) # drydep,

g = graph(csys)
f, ax, p = graphplot(g; ilabels=labels(g))
hidedecorations!(ax); hidespines!(ax)
f

st = SolverStrangSerial(Rosenbrock23(autodiff=false), dt)

prob = ODEProblem(csys, st; tspan=BigFloat.(get_tspan(domain)))

@time sol = solve(prob, SSPRK22(); dt=dt, progress=true, progress_steps=1,
    save_on=false, save_start=false, save_end=false, initialize_save=false)

ds = NCDataset(outfile, "r");

anim = @animate for i ∈ 1:size(ds["SuperFast₊O3"])[4]
    Plots.plot(
        Plots.heatmap(ds["SuperFast₊O3"][:, :, 1, i]', title="Ground-Level"),
        Plots.heatmap(ds["SuperFast₊O3"][:, 2, :, i]', title="Vertical Cross-Section"),
        size=(1200, 400)
    )
end
gif(anim, fps = 5)
```