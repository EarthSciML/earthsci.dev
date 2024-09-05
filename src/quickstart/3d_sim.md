# 3D Simulation

```@example 3d_sim
using EarthSciMLBase, EarthSciData, GasChem, EnvironmentalTransport
using ModelingToolkit, DifferentialEquations
using ModelingToolkit: t
using DomainSets, DynamicQuantities, Dates
using NCDatasets, Plots
using SciMLBase: DiscreteCallback
using ProgressMeter

@parameters lon = deg2rad(-97) [unit=u"rad"]
@parameters lat = deg2rad(40) [unit=u"rad"]
@parameters lev = 1
emis = NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", lon, lat, lev; dtype=Float64)

geosfp, geosfp_updater = GEOSFP("0.5x0.625_NA"; dtype = Float64,
    coord_defaults = Dict(:lon => deg2rad(-97), :lat => deg2rad(40), :lev => 1.0))

starttime = datetime2unix(DateTime(2016, 5, 1, 0, 0))
endtime = datetime2unix(DateTime(2016, 5, 5, 0, 0))

domain = DomainInfo(
    [partialderivatives_δxyδlonlat,
        partialderivatives_δPδlev_geosfp(geosfp)],
    constIC(16.0, t ∈ Interval(starttime, endtime)),
    constBC(16.0,
        lon ∈ Interval(deg2rad(-115), deg2rad(-68.75)),
        lat ∈ Interval(deg2rad(25), deg2rad(53.7)),
        lev ∈ Interval(1, 15)),
    dtype = Float64)

chem = SuperFast()
photolysis = FastJX()

outfile = ("RUNNER_TEMP" ∈ keys(ENV) ? ENV["RUNNER_TEMP"] : tempname()) * "out.nc" # This is just a location to save the output.
output = NetCDFOutputter(outfile, 3600.0)

# This simulation doesn't work unless we get rid of this coupling.
function EarthSciMLBase.couple2(c::GasChem.SuperFastCoupler, g::EarthSciData.GEOSFPCoupler)
    c, g = c.sys, g.sys
    ConnectorSystem([], c, g)
end


csys = couple(chem, photolysis, geosfp, geosfp_updater, emis, domain, output)

adv = AdvectionOperator(100.0, upwind1_stencil, ZeroGradBC())

csys = couple(csys, adv)

function pbar(start, finish)
    p = Progress(Int(round(finish-start)))
    DiscreteCallback(
        (_, _, _) -> true, 
        (integrator) -> update!(p, Int(round(integrator.t-start)));
        save_positions = (false, false),
    )
end
#csys = couple(csys, pbar(starttime, endtime))


sim = Simulator(csys, [deg2rad(2), deg2rad(2), 1])
st = SimulatorStrangThreads(Rosenbrock23(), SSPRK22(), 300.0)

u0 = init_u(sim);
u0 .+= rand(size(u0)) * 1e-6
@time run!(sim, st, u0, save_on=false, save_start=false, save_end=false,
    initialize_save=false)

ds = NCDataset(outfile, "r")

anim = @animate for i ∈ 1:size(ds["SuperFast₊O3"])[4]
    plot(
        heatmap(ds["SuperFast₊O3"][:, :, 1, i]', title="Ground-Level"),
        heatmap(ds["SuperFast₊O3"][:, 2, :, i]', title="Vertical Cross-Section"),
    )
end
gif(anim, fps = 15)
```
