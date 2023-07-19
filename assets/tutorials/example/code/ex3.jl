# This file was generated, do not modify it. # hide
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