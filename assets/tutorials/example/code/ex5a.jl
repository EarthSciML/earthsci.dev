# This file was generated, do not modify it. # hide
domain = DomainInfo(
    partialderivatives_lonlat2xymeters,
    constIC(1.0f0, t ∈ Interval(start, finish)),
    zerogradBC(lon ∈ Interval(-130f0, -100f0)))

geos = GEOSFP("0.25x0.3125_NA", t; coord_defaults = Dict(:lat => 34.0, :lev => 1))

model = SuperFast(t) + FastJX(t) + domain +
    Emissions(t, -118.2, 0.2)+ Advection() + geos

g = graph(model)

f, ax, p = graphplot(g; ilabels=labels(g))
hidedecorations!(ax); hidespines!(ax); ax.aspect = DataAspect()

save(joinpath(@OUTPUT, "graph.svg"), f) # hide