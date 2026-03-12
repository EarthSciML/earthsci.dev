# Climate Sensitivity Example

This document demonstrates how to use the **ClimateSensitivity** module to perform a simple ensemble-based climate parameter estimation and projection workflow. The example loads scenario forcing inputs, generates pseudo-observations from a known “truth,” assimilates those observations, and then compares prior and posterior temperature projections.

The workflow has four main steps. First, we define the time range, model settings, and prior parameter specification. Second, we generate synthetic observations of temperature and ocean heat uptake that act as our calibration targets. Third, we run the ensemble data assimilation procedure to estimate posterior parameter samples. Finally, we propagate both prior and posterior ensembles forward in time and visualize how the uncertainty in projected temperature evolves.

## Load Packages and Module

We begin by importing the necessary libraries. 

```@example climate_sensitivity
using Random
using CSV
using DataFrames
using Statistics
using Plots

include(joinpath(@__DIR__, "..", "ClimateSensitivity.jl"))
using .ClimateSensitivity

```
## Configuration

Next, we define the full simulation years and the subset of years used for assimilation. We also specify the ensemble size, random seed, and other configuration parameters such as observational noise levels, process noise, and optimization settings.

```@example climate_sensitivity

years_full  = collect(1850:2100)
years_assim = collect(1850:2050)

Nens = parse(Int, get(ENV, "NENS", "500"))
seed = parse(Int, get(ENV, "SEED", "1234"))

params = (
    years        = years_full,
    assim_years  = years_assim,
    dt           = 1.0,
    co2_C0_ppm   = 278.0,
    phi_q        = 0.2,
    sigma_q      = 0.27,
    sigma_Tobs   = 0.10,
    sigma_Qobs   = 0.1,
    maxiters     = 500,
    accept_J_max = 1e2,
)

ps = make_prior_spec(; params...)
rng = MersenneTwister(seed)
```

## Load Scenario Forcing Inputs

Here we load the scenario input data from scenario_inputs.csv. This file provides the time series of atmospheric CO₂ concentration and SO₂ emissions used to force the climate model.

```@example climate_sensitivity

data_path = joinpath(@__DIR__, "data", "scenario_inputs.csv")
if !isfile(data_path)
    error("The data file 'scenario_inputs.csv' was not found in the data folder.")
end

@info "Using scenario file: $(data_path)"

df = CSV.read(data_path, DataFrame)

function align_series(df::DataFrame, years::AbstractVector{<:Integer}, col::Symbol)
    lookup = Dict(Int.(df.year) .=> Float64.(df[!, col]))
    return [lookup[y] for y in years]
end

co2_full = align_series(df, years_full, :co2_ppm)
so2_full = align_series(df, years_full, :so2_mt_per_yr)
```

## Generate Pseudo-Observations (Truth)

Here we generate synthetic observations from a known parameter vector. In this example, the “truth” is chosen as the prior mean. The model is then run forward using that truth to create pseudo-observations of temperature and ocean heat uptake.

```@example climate_sensitivity


x_true = ClimateSensitivity.prior_mean(ps)
T10, T20, λ, γ, C1, C2, ε, f1, f3, f1s, f2s, C0s, q = ClimateSensitivity.unpack_x(x_true, ps)
F2x = ClimateSensitivity.forcing_doubling_CO2(; f1_CO2 = f1, f3_CO2 = f3, C0_CO2 = params.co2_C0_ppm)

@info "Generating pseudo-observations with correct TCR..."
T_obs, Q_obs_ZJ, F_assim_true = generate_pseudo_obs(
    rng,
    ps,
    co2_full,
    so2_full;
    x_true = x_true,
    params...,
)
```
## Restrict Observations to the Modern Period

In this example, we only retain pseudo-observations from 2020 onward. Earlier years are set to NaN so that they do not influence the assimilation.

```@example climate_sensitivity

for (i, y) in enumerate(years_assim)
    if y < 2020
        T_obs[i] = NaN
        Q_obs_ZJ[i] = NaN
    end
end

```

## Estimate the Posterior

Starting from the prior specification, the algorithm updates the ensemble members so that they better match the pseudo-observations while remaining consistent with the prior assumptions.

```@example climate_sensitivity

@info "Running RML-EDA with Nens=$(Nens)..."

X_post, J_post = ClimateSensitivity.rml_eda_parallel(
    rng,
    Nens,
    ps,
    years_assim,
    F_assim_true,
    T_obs,
    Q_obs_ZJ;
    perturb_obs = true,
    verbose = true,
    nworkers = Threads.nthreads(),
    params...,
)
```

## Run Ensembles (Prior & Posterior)

Now that the posterior parameter ensemble has been obtained, we run the climate model forward for both the prior and posterior ensembles. This allows us to compare the projected temperature response before and after assimilating the pseudo-observations.

```@example climate_sensitivity

@info "Running forecast ensembles (2020–2100)..."

X_prior = hcat([sample_prior(rng, ps) for _ in 1:Nens]...)

Y_prior, ECS_prior, TCR_prior =
    forecast_ensemble(ps, X_prior, co2_full, so2_full; params...)

Y_post, ECS_post, TCR_post =
    forecast_ensemble(ps, X_post, co2_full, so2_full; params...)
```

##  Plot Temperature Projections

Finally, we summarize the prior and posterior temperature ensembles and plot the median and 5–95% uncertainty ranges for the temperature projections from 2020 to 2100.

```@example climate_sensitivity

function quantiles_over_ensemble(years::AbstractVector, Y::AbstractMatrix)
    p05 = [quantile(view(Y, i, :), 0.05) for i in axes(Y, 1)]
    p50 = [quantile(view(Y, i, :), 0.50) for i in axes(Y, 1)]
    p95 = [quantile(view(Y, i, :), 0.95) for i in axes(Y, 1)]
    return DataFrame(year = years, p05 = p05, p50 = p50, p95 = p95)
end

q_prior = quantiles_over_ensemble(years_full, Y_prior)
q_post  = quantiles_over_ensemble(years_full, Y_post)

obs_df = DataFrame(year = years_assim, T_obs = T_obs, Q_obs_ZJ = Q_obs_ZJ)

df_et = DataFrame(
    sample_type = vcat(fill("prior", Nens), fill("posterior", Nens)),
    ECS = vcat(ECS_prior, ECS_post),
    TCR = vcat(TCR_prior, TCR_post),
)

q_prior_2020 = filter(row -> row.year >= 2020, q_prior)
q_post_2020  = filter(row -> row.year >= 2020, q_post)
obs_2020     = filter(row -> row.year >= 2020 && !isnan(row.T_obs), obs_df)

prior_lower = q_prior_2020.p50 .- q_prior_2020.p05
prior_upper = q_prior_2020.p95 .- q_prior_2020.p50
post_lower  = q_post_2020.p50 .- q_post_2020.p05
post_upper  = q_post_2020.p95 .- q_post_2020.p50

p = plot(
    q_prior_2020.year,
    q_prior_2020.p50;
    ribbon = (prior_lower, prior_upper),
    label = "Prior median",
    fillalpha = 0.20,
    color = :thistle,
    linestyle = :dash,
    linewidth = 2,
    xlabel = "Year",
    ylabel = "Global Surface Temperature Anomaly (°C)",
    title = "\nTemperature Projections (2020–2100)",
    legend = :topleft,
    xlims = (2020, 2100),
    ylims = (0.0, 5.0),
    size = (900, 500),
    left_margin = 8Plots.mm,
    right_margin = 8Plots.mm,
    top_margin = 8Plots.mm,
    bottom_margin = 8Plots.mm
)

plot!(
    p,
    q_post_2020.year,
    q_post_2020.p50;
    ribbon = (post_lower, post_upper),
    label = "Posterior median",
    fillalpha = 0.25,
    color = :steelblue,
    linewidth = 2,
)

plot!(
    p,
    obs_2020.year,
    obs_2020.T_obs;
    label = "Truth (Pseudo-obs)",
    color = :black,
    linestyle = :dash,
    linewidth = 2,
)

vline!(p, [2050]; label = "Forecast start", color = :gray, linestyle = :dot, linewidth = 2)
annotate!(p, 2051, 0.6, text("Forecast starts", 10, :gray))

p
```