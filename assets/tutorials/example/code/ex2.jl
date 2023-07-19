# This file was generated, do not modify it. # hide
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
ticks = Dates.datetime2unix.([Dates.DateTime(2022, 5, 1), Dates.DateTime(2022, 5, 2), 
    Dates.DateTime(2022, 5, 3)])
tickstrings = Dates.format.(Dates.unix2datetime.(ticks), "m-d-yy")
Plots.plot(sol,ylims=(0,30),xlabel="Date", ylabel="Concentration (ppb)", 
    ylabelfontsize=9, xlabelfontsize=9, 
    xticks=(ticks, tickstrings), legend=:outertopright, size=(500, 310))
Plots.savefig(joinpath(@OUTPUT, "ode.svg")) # hide