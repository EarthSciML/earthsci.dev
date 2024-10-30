# 3D Simulation

```@example 3d_sim
using EarthSciMLBase, EarthSciData, GasChem, EnvironmentalTransport
using ModelingToolkit, DifferentialEquations
using Dates
using NCDatasets, Plots
using ProgressLogging # Needed for progress bar. Use `TerminalLoggers` if in a terminal.

domain = DomainInfo(
    DateTime(2016, 5, 1),
    DateTime(2016, 5, 1, 2);
    lonrange = deg2rad(-115):deg2rad(5):deg2rad(-68.75),
    latrange = deg2rad(25):deg2rad(4):deg2rad(53.7),
    levrange = 1:15,
    dtype = Float64)

emis = NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", domain)

geosfp = GEOSFP("0.5x0.625_NA", domain)

domain = EarthSciMLBase.add_partial_derivative_func(domain, partialderivatives_δPδlev_geosfp(geosfp))

chem = SuperFast()
photolysis = FastJX()

dt = 60.0 # Splitting timestep
adv = AdvectionOperator(dt, upwind1_stencil, ZeroGradBC())

outfile = ("RUNNER_TEMP" ∈ keys(ENV) ? ENV["RUNNER_TEMP"] : tempname()) * "out.nc" # This is just a location to save the output.
output = NetCDFOutputter(outfile, 3600.0)

csys = couple(chem, photolysis, geosfp, adv, domain, output) # emis, 

st = SolverStrangSerial(Rosenbrock23(), dt)
prob = ODEProblem(csys, st)
sol = solve(prob, SSPRK22(); dt=dt, progress=true, progress_steps=1,
    save_on=false, save_start=false, save_end=false, initialize_save=false, abstol=1e-8, reltol=1e-8)

ds = NCDataset(outfile, "r")

anim = @animate for i ∈ 1:size(ds["SuperFast₊O3"])[4]
    plot(
        heatmap(ds["SuperFast₊O3"][:, :, 1, i]', title="Ground-Level"),
        heatmap(ds["SuperFast₊O3"][:, 2, :, i]', title="Vertical Cross-Section"),
    )
end
gif(anim, fps = 15)
```