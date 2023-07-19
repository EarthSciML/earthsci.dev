# This file was generated, do not modify it. # hide
discrete_lon = sol[lon]
discrete_t = sol[t]

@variables superfast₊O3(..) superfast₊NO(..) superfast₊ISOP(..) superfast₊NO2(..)
sol_isop = sol[superfast₊ISOP(t, lon)]
sol_o3 = sol[superfast₊O3(t, lon)]
sol_no = sol[superfast₊NO(t, lon)]
sol_no2 = sol[superfast₊NO2(t, lon)]

using LaTeXStrings
Plots.plot(
    Plots.heatmap(discrete_t, discrete_lon[1:end], sol_isop[:, 1:end]', 
        xticks=(ticks, tickstrings), 
        ylabel="Longitude (°)", title="Isoprene", titlefontsize=12,
        margin=0Plots.mm,
        size=(600, 400), dpi=300),
    Plots.heatmap(discrete_t, discrete_lon[1:end], sol_no[:, 1:end]', 
        xticks=(ticks, tickstrings), title="NO", titlefontsize=12),
    Plots.heatmap(discrete_t, discrete_lon[1:end], sol_no2[:, 1:end]', 
        xticks=(ticks, tickstrings), xlabel="Date", title=L"\textrm{NO_2}"),
    Plots.heatmap(discrete_t, discrete_lon[1:end], sol_o3[:, 1:end]',colorbar_title="Concentration (ppb)",
        xticks=(ticks, tickstrings), title=L"\textrm{O_3}", titlefontsize=12),
)
Plots.savefig(joinpath(@OUTPUT, "pde.svg"))