module ClimateSensitivity

"""
The `ClimateSensitivity` module defines a two-layer energy-balance climate model together with:

- a prior distribution over uncertain parameters,
- forcing functions for CO2 and SO2,
- functions for fitting model parameters to observation,
- utilities for ensemble forecasts and output.

The main control vector is:

x = [T10, T20, λ, γ, C1, C2, ε, f1_CO2, f3_CO2, f1_SO2, f2_SO2, C0_SO2, q1, ..., qn]

where:

- T10 = initial upper-layer temperature [K]
- T20 = initial deep-layer temperature [K]
- λ = climate feedback parameter [W K^-1]
- γ = heat exchange between layers [W K^-1]
- C1 = upper-layer heat capacity [J K^-1]
- C2 = deep-layer heat capacity [J K^-1]
- ε = ocean heat uptake efficacy [-]
- f1_CO2 = CO2 forcing parameter [W m^-2], 
- f3_CO2 = CO2 forcing parameter [W m^-2 (ppm CO₂)^-1/2]
- f1_SO2 = SO2 forcing parameter [W]
- f2_SO2 = SO2 forcing parameter [W (MtSO₂ / yr)^-1]
- C0_SO2 = SO2 forcing parameter [MtSO₂ / yr]
- q1 ... qn = yearly stochastic forcing perturbations

The two-layer model is:

dT1/dt = [F(t) + q(t) - λ*T1 - ε*γ*(T1 - T2)] / C1

dT2/dt = [γ*(T1 - T2)] / C2

"""

export PriorSpec,
       make_prior_spec, sample_prior, prior_mean,
       unpack_x, Q_to_ZJ,
       forcing_CO2, forcing_SO2, forcing_CO2_SO2, forcing_total,
       integrate_two_layer, compute_ECS_TCR, forcing_doubling_CO2,
       ar1_covariance, sample_ar1,
       rml_eda, rml_eda_parallel, quantiles_over_ensemble,
       generate_pseudo_obs, forecast_ensemble,
       write_ensemble_csvs

using LinearAlgebra
using Random
using Base.Threads
using Statistics
using Optim
using ForwardDiff
using CSV
using DataFrames
using ModelingToolkit
using DifferentialEquations
using DataInterpolations

"""
This struct stores:

- μ = prior mean vector
- P = prior covariance matrix
- L = Cholesky factor of `P`, so that P = L * L'
- idx_q = index range of the yearly forcing perturbation terms q

"""
struct PriorSpec
    μ::Vector{Float64}
    P::Matrix{Float64}
    L::LowerTriangular{Float64,Matrix{Float64}}
    idx_q::UnitRange{Int}
end

"""
Builds and return a PriorSpec object containing the prior mean vector,
prior covariance matrix, Cholesky factor, and index range for the yearly
model errors, whose covariance is specified with an AR(1) structure.
"""
function make_prior_spec(; assim_years::AbstractVector{<:Integer},
                           sigma_q::Real = 0.27,
                           phi_q::Real = 0.2,
                           kwargs...)
    n = length(assim_years)

    # Prior means
    T10     = 0.0
    T20     = 0.0
    λ       = 1.258
    γ       = 0.7
    C1      = 8.0
    C2      = 100.0
    ε       = 1.58
    f1_CO2  = 4.58
    f3_CO2  = 0.086
    f1_SO2  = -0.0047
    f2_SO2  = -0.0096
    C0_SO2  = 170.6

    μ_base = [
        T10, T20,
        λ, γ, C1, C2, ε,
        f1_CO2, f3_CO2,
        f1_SO2, f2_SO2, C0_SO2,
    ]
    μ = vcat(μ_base, zeros(n))

    # Prior standard deviations
    σ_T10     = 0.2
    σ_T20     = 0.2
    σ_λ       = 0.38
    σ_γ       = 0.21
    σ_C1      = 2.4
    σ_C2      = 30.0
    σ_ε       = 0.128
    σ_f1_CO2  = 0.519
    σ_f3_CO2  = 0.026
    σ_f1_SO2  = 0.0014
    σ_f2_SO2  = 0.0029
    σ_C0_SO2  = 51.2

    σ_base = [
        σ_T10, σ_T20,
        σ_λ, σ_γ, σ_C1, σ_C2, σ_ε,
        σ_f1_CO2, σ_f3_CO2,
        σ_f1_SO2, σ_f2_SO2, σ_C0_SO2,
    ]

    # covariance matrix
    P = zeros(length(μ), length(μ))
    for i in 1:length(σ_base)
        P[i, i] = σ_base[i]^2
    end

    C_q = Matrix{Float64}(undef, n, n)
    for i in 1:n, j in 1:n
        C_q[i, j] = (sigma_q^2) * (phi_q^abs(i - j))
    end

    idx_q = (length(μ_base) + 1):(length(μ_base) + n)
    P[idx_q, idx_q] .= Symmetric(C_q)

    # Cholesky factor
    L = cholesky(Symmetric(P) + 1e-10I).L
    return PriorSpec(μ, P, L, idx_q)
end

# Returns the prior-mean state vector.
prior_mean(ps::PriorSpec) = ps.μ

# Draws one sample from the prior distribution.
function sample_prior(rng::AbstractRNG, ps::PriorSpec)
    z = randn(rng, length(ps.μ))
    return ps.μ .+ (ps.L * z)
end

# Unpacks the state vector into model parameters and yearly forcing perturbations.
function unpack_x(x::AbstractVector, ps::PriorSpec)
    T10 = x[1]
    T20 = x[2]
    λ   = x[3]
    γ   = x[4]
    C1  = x[5]
    C2  = x[6]
    ε   = x[7]
    f1_C = x[8]
    f3_C = x[9]
    f1_S = x[10]
    f2_S = x[11]
    C0_S = x[12]
    q = x[ps.idx_q]
    return T10, T20, λ, γ, C1, C2, ε, f1_C, f3_C, f1_S, f2_S, C0_S, q
end

# Constructs the covariance matrix for an AR(1)-style yearly forcing perturbation process.
function ar1_covariance(n::Int, φ::Float64, σ::Float64)
    C = Matrix{Float64}(undef, n, n)
    for i in 1:n, j in 1:n
        C[i, j] = (σ^2) * (abs(i - j) == 0 ? 1.0 : φ^abs(i - j))
    end
    return Symmetric(C)
end

# Samples one AR(1) perturbation time series.
function sample_ar1(rng::AbstractRNG, n::Int, φ::Float64, σ::Float64)
    q = zeros(n)
    q[1] = σ * randn(rng)
    α = sqrt(max(0.0, 1 - φ^2)) * σ
    for t in 1:(n - 1)
        q[t + 1] = φ * q[t] + α * randn(rng)
    end
    return q
end

# --- Forcing ---

# Computes radiative forcing from atmospheric CO₂ concentration.
@inline function forcing_CO2(C_ppm::Real; f1, f2 = 0.0, f3, C0_ppm)
    (C_ppm <= 0) && return NaN
    return f1 * log(C_ppm / C0_ppm) + f2 * (C_ppm - C0_ppm) + f3 * (sqrt(C_ppm) - sqrt(C0_ppm))
end

# Computes radiative forcing from SO₂ emissions.
@inline function forcing_SO2(E_mtpy::Real; f1, f2, C0_mtpy)
    val = 1 + E_mtpy / C0_mtpy
    (val <= 0) && return NaN
    return f1 * log(val) + f2 * E_mtpy
end

# Combines CO₂ and SO₂ forcing histories into a single annual forcing time series.
function forcing_CO2_SO2(
    years::AbstractVector{<:Integer},
    co2_ppm::AbstractVector{<:Real},
    so2_mtpy::AbstractVector{<:Real};
    f1_CO2,
    f3_CO2,
    C0_CO2,
    f1_SO2,
    f2_SO2,
    C0_SO2,
)
    length(years) == length(co2_ppm) == length(so2_mtpy) ||
        throw(ArgumentError("years, co2_ppm, and so2_mtpy must have the same length."))

    F = Vector{Float64}(undef, length(co2_ppm))
    for i in eachindex(years)
        F[i] =
            forcing_CO2(co2_ppm[i]; f1 = f1_CO2, f3 = f3_CO2, C0_ppm = C0_CO2) +
            forcing_SO2(so2_mtpy[i]; f1 = f1_SO2, f2 = f2_SO2, C0_mtpy = C0_SO2)
    end
    return F
end

forcing_total = forcing_CO2_SO2

# Computes the forcing for a doubling of CO₂ from its reference concentration.
function forcing_doubling_CO2(; f1_CO2, f3_CO2, C0_CO2)
    return forcing_CO2(2 * C0_CO2; f1 = f1_CO2, f3 = f3_CO2, C0_ppm = C0_CO2)
end

@register_symbolic eval_forcing(t, itp::Any)
eval_forcing(t, itp) = itp(t)

# Builds two-layer energy-balance model used for time integration.
function build_mtk_system()
    @parameters t
    @parameters λ γ C1 C2 ε
    @parameters itp_F::Any itp_q::Any

    @variables T1(t) T2(t)
    D = Differential(t)

    F_t = eval_forcing(t, itp_F)
    q_t = eval_forcing(t, itp_q)

    eqs = [
        D(T1) ~ (F_t + q_t - λ * T1 - ε * γ * (T1 - T2)) / C1,
        D(T2) ~ (γ * (T1 - T2)) / C2,
    ]

    return ODESystem(
        eqs,
        t,
        [T1, T2],
        [λ, γ, C1, C2, ε, itp_F, itp_q];
        name = :TwoLayerEnergyBalance,
    )
end

const _SYS = build_mtk_system()
const _SIM_SYS = structural_simplify(_SYS)

# Integrates the two-layer energy-balance model over the supplied yearly forcing history.
function integrate_two_layer(
    years::AbstractVector{<:Integer},
    F::AbstractVector{<:Real};
    T10::Real,
    T20::Real,
    Q0::Real,
    λ::Real,
    γ::Real,
    C1::Real,
    C2::Real,
    ε::Real,
    q::AbstractVector,
    dt::Real = 1.0,
    kwargs...,
)
    n = length(years)
    length(F) == n || throw(ArgumentError("F length $(length(F)) must match years length $n."))
    length(q) == n || throw(ArgumentError("q length $(length(q)) must match years length $n."))

    Fv = Float64.(F)
    qv = Float64.(q)

    if any(!isfinite, Fv)
        error("F contains NaN/Inf.")
    end
    if any(!isfinite, qv)
        error("q contains NaN/Inf.")
    end
    if C1 <= 1e-2 || C2 <= 1e-2
        return fill(NaN, n), fill(NaN, n), fill(NaN, n)
    end

    t_points = Float64.(years)
    t_itp = vcat(t_points, t_points[end] + Float64(dt))

    itp_F = ConstantInterpolation(vcat(Fv, Fv[end]), t_itp)
    itp_q = ConstantInterpolation(vcat(qv, qv[end]), t_itp)

    u0 = [
        _SIM_SYS.T1 => Float64(T10),
        _SIM_SYS.T2 => Float64(T20),
    ]

    p_vals = [
        _SIM_SYS.λ     => Float64(λ),
        _SIM_SYS.γ     => Float64(γ),
        _SIM_SYS.C1    => Float64(C1),
        _SIM_SYS.C2    => Float64(C2),
        _SIM_SYS.ε     => Float64(ε),
        _SIM_SYS.itp_F => itp_F,
        _SIM_SYS.itp_q => itp_q,
    ]

    op = merge(Dict(u0), Dict(p_vals))
    prob = ODEProblem(_SIM_SYS, op, (t_points[1], t_points[end]))

    sol = solve(prob, Rodas4(); saveat = t_points, abstol = 1e-8, reltol = 1e-8)

    if sol.retcode != SciMLBase.ReturnCode.Success
        return fill(NaN, n), fill(NaN, n), fill(NaN, n)
    end

    T1_sol = sol[_SIM_SYS.T1]
    T2_sol = sol[_SIM_SYS.T2]

    if any(!isfinite, T1_sol) || any(!isfinite, T2_sol)
        return fill(NaN, n), fill(NaN, n), fill(NaN, n)
    end

    Q = fill(Float64(Q0), n)
    return T1_sol, T2_sol, Q
end

# Converts heat uptake from W·yr m⁻² to zettajoules.
Q_to_ZJ(Q_Wyr_m2) = Q_Wyr_m2 * (365.0 * 24 * 3600) * 5.10072e14 / 1e21

# Compute ECS, TCR, and doubled-CO₂ forcing from the current parameter set.
function compute_ECS_TCR(; λ, γ, f1_CO2, f3_CO2, C0_CO2)
    F2x = forcing_doubling_CO2(; f1_CO2 = f1_CO2, f3_CO2 = f3_CO2, C0_CO2 = C0_CO2)
    return F2x / λ, F2x / (λ + γ), F2x
end


# --- Optimization ---

# Evaluates the prior-deviation term of the objective function.
@inline function J_guess(x::AbstractVector, χ::AbstractVector, ps::PriorSpec)
    d = x .- χ
    z = ps.L \ d
    return 0.5 * dot(z, z)
end

# Evaluates the observation-deviation term for surface temperature and heat uptake.
function J_obs(
    T1_sim,
    Q_sim,
    T_obs,
    Q_obs_ZJ;
    sigma_Tobs::Real = 0.10,
    sigma_Qobs::Real = 10.0,
    kwargs...,
)
    length(T1_sim) == length(T_obs) || return Inf

    rT = T1_sim .- T_obs
    rQ = Q_to_ZJ.(Q_sim) .- Q_obs_ZJ
    wT = 1.0 / (sigma_Tobs^2)
    wQ = 1.0 / (sigma_Qobs^2)

    valid_T = .!isnan.(rT)
    valid_Q = .!isnan.(rQ)

    return 0.5 * (wT * sum(abs2, rT[valid_T]) + wQ * sum(abs2, rQ[valid_Q]))
end

# Evaluates the full objective function for one ensemble member.
function cost_function(
    x::AbstractVector,
    χ::AbstractVector,
    years::AbstractVector{<:Integer},
    F::AbstractVector{<:Real},
    T_obs::AbstractVector{<:Real},
    Q_obs_ZJ::AbstractVector{<:Real},
    ps::PriorSpec;
    dt::Real = 1.0,
    sigma_Tobs::Real = 0.10,
    sigma_Qobs::Real = 10.0,
    kwargs...,
)
    T10, T20, λ, γ, C1, C2, ε, f1_C, f3_C, f1_S, f2_S, C0_S, q = unpack_x(x, ps)

    T1, T2, Q = integrate_two_layer(
        years,
        F;
        T10 = T10,
        T20 = T20,
        Q0 = 0.0,
        λ = λ,
        γ = γ,
        C1 = C1,
        C2 = C2,
        ε = ε,
        q = q,
        dt = dt,
    )

    if any(isnan, T1) || any(isnan, Q)
        return Inf
    end

    Jg = J_guess(x, χ, ps)
    Jo = J_obs(T1, Q, T_obs, Q_obs_ZJ; sigma_Tobs = sigma_Tobs, sigma_Qobs = sigma_Qobs)

    J = Jg + Jo
    return (isnan(J) || isinf(J)) ? Inf : J
end

# Optimizes to find the best-fitting parameter set for one ensemble member.
function rml_member(
    rng::AbstractRNG,
    ps::PriorSpec,
    years::AbstractVector{<:Integer},
    F::AbstractVector{<:Real},
    T_obs::AbstractVector{<:Real},
    Q_obs_ZJ::AbstractVector{<:Real};
    perturb_obs::Bool = true,
    dt::Real = 1.0,
    sigma_Tobs::Real = 0.10,
    sigma_Qobs::Real = 10.0,
    maxiters::Int = 200,
    kwargs...,
)
    χ = sample_prior(rng, ps)
    T_obs_i = collect(Float64, T_obs)
    Q_obs_i = collect(Float64, Q_obs_ZJ)

    if perturb_obs
        vt = .!isnan.(T_obs_i)
        vq = .!isnan.(Q_obs_i)
        T_obs_i[vt] .+= sigma_Tobs .* randn(rng, count(vt))
        Q_obs_i[vq] .+= sigma_Qobs .* randn(rng, count(vq))
    end

    obj(x) = cost_function(
        x,
        χ,
        years,
        F,
        T_obs_i,
        Q_obs_i,
        ps;
        dt = dt,
        sigma_Tobs = sigma_Tobs,
        sigma_Qobs = sigma_Qobs,
    )

    function g!(G, x)
        try
            ForwardDiff.gradient!(G, obj, x)
        catch
            fill!(G, 0.0)
        end
    end

    options = Optim.Options(iterations = maxiters, show_trace = false, g_tol = 1e-8)

    try
        res = Optim.optimize(obj, g!, χ, Optim.LBFGS(), options)
        return Optim.minimizer(res), Optim.minimum(res), Optim.converged(res)
    catch
        return χ, Inf, false
    end
end

# Runs the optimization for ensemble members in parallel.
function rml_eda_parallel(
    rng::AbstractRNG,
    N::Int,
    ps::PriorSpec,
    years::AbstractVector{<:Integer},
    F::AbstractVector{<:Real},
    T_obs::AbstractVector{<:Real},
    Q_obs_ZJ::AbstractVector{<:Real};
    perturb_obs::Bool = true,
    verbose::Bool = true,
    nworkers::Int = Threads.nthreads(),
    accept_J_max::Real = 1e2,
    dt::Real = 1.0,
    sigma_Tobs::Real = 0.10,
    sigma_Qobs::Real = 10.0,
    maxiters::Int = 200,
    kwargs...,
)
    d = length(ps.μ)
    X_post = Matrix{Float64}(undef, d, N)
    J_post = Vector{Float64}(undef, N)

    try
        LinearAlgebra.BLAS.set_num_threads(1)
    catch
    end

    stop_flag = Threads.Atomic{Int}(0)
    attempts = Threads.Atomic{Int}(0)
    ch = Channel{Tuple{Vector{Float64}, Float64}}(N + nworkers)

    base_seed = rand(rng, UInt)

    function worker(worker_id::Int)
        local_rng = MersenneTwister(hash((base_seed, worker_id)))
        while stop_flag[] == 0
            Threads.atomic_add!(attempts, 1)

            x, J, conv = rml_member(
                local_rng,
                ps,
                years,
                F,
                T_obs,
                Q_obs_ZJ;
                perturb_obs = perturb_obs,
                dt = dt,
                sigma_Tobs = sigma_Tobs,
                sigma_Qobs = sigma_Qobs,
                maxiters = maxiters,
            )

            if conv && isfinite(J) && J <= accept_J_max
                if stop_flag[] == 0
                    put!(ch, (x, J))
                end
            end
        end
        return nothing
    end

    if verbose
        println("Starting parallel assimilation of $N members with $nworkers workers...")
    end

    tasks = [Threads.@spawn worker(w) for w in 1:nworkers]

    accepted = 0
    while accepted < N
        x, J = take!(ch)
        accepted += 1
        X_post[:, accepted] = x
        J_post[accepted] = J

        if verbose
            print("\rProgress: $accepted / $N  (attempts=$(attempts[]))")
        end
    end

    stop_flag[] = 1

    for t in tasks
        wait(t)
    end

    if verbose
        println("\nDone.")
    end

    return X_post, J_post
end

# --- Output ---

# Generates synthetic observations.
function generate_pseudo_obs(
    rng::AbstractRNG,
    ps::PriorSpec,
    co2,
    so2;
    x_true = nothing,
    years,
    assim_years,
    dt::Real = 1.0,
    co2_C0_ppm::Real = 278.0,
    phi_q::Real = 0.2,
    sigma_q::Real = 0.27,
    kwargs...,
)
    μ = isnothing(x_true) ? prior_mean(ps) : x_true
    T10, T20, λ, γ, C1, C2, ε, f1c, f3c, f1s, f2s, C0s, _ = unpack_x(μ, ps)

    F_full = forcing_CO2_SO2(
        years,
        co2,
        so2;
        f1_CO2 = f1c,
        f3_CO2 = f3c,
        C0_CO2 = co2_C0_ppm,
        f1_SO2 = f1s,
        f2_SO2 = f2s,
        C0_SO2 = C0s,
    )

    i0 = findfirst(==(assim_years[1]), years)
    i1 = findfirst(==(assim_years[end]), years)
    isnothing(i0) && error("assim_years[1] was not found in years.")
    isnothing(i1) && error("assim_years[end] was not found in years.")

    F_assim = F_full[i0:i1]
    q_assim = sample_ar1(rng, length(assim_years), Float64(phi_q), Float64(sigma_q))

    T1, _, Q = integrate_two_layer(
        assim_years,
        F_assim;
        T10 = T10,
        T20 = T20,
        Q0 = 0.0,
        λ = λ,
        γ = γ,
        C1 = C1,
        C2 = C2,
        ε = ε,
        q = q_assim,
        dt = dt,
    )

    return copy(T1), [Q_to_ZJ(q) for q in Q], F_assim
end

# Evaluate the explicit two-layer model right-hand side for fast ensemble solves.
@inline function two_layer_rhs!(du, u, p, t)
    T1 = u[1]
    T2 = u[2]
    F = p.itp_F(t)

    du[1] = (F - p.λ * T1 - p.ε * p.γ * (T1 - T2)) / p.C1
    du[2] = (p.γ * (T1 - T2)) / p.C2
    return nothing
end

# Run forecasts for all ensemble members.
function forecast_ensemble(
    ps::PriorSpec,
    X_ens,
    co2,
    so2;
    years,
    dt::Real = 1.0,
    co2_C0_ppm::Real = 278.0,
    alg = Tsit5(),
    ensemblealg = SciMLBase.EnsembleThreads(),
    abstol = 1e-8,
    reltol = 1e-8,
    batch_size = size(X_ens, 2),
    kwargs...,
)
    t_pts = Float64.(years)
    n_t = length(t_pts)
    n_ens = size(X_ens, 2)
    t_itp = vcat(t_pts, t_pts[end] + Float64(dt))

    itp_F0 = DataInterpolations.ConstantInterpolation(vcat(0.0, 0.0), [t_pts[1], t_pts[1] + Float64(dt)])

    p0 = (λ = 1.0, γ = 0.7, C1 = 8.0, C2 = 100.0, ε = 1.58, itp_F = itp_F0)
    prob0 = ODEProblem(two_layer_rhs!, [0.0, 0.0], (t_pts[1], t_pts[end]), p0)

    prob_func = function (prob, i, repeat)
        x = @view X_ens[:, i]
        T10, T20, λ, γ, C1, C2, ε, f1c, f3c, f1s, f2s, C0s, _ = unpack_x(x, ps)

        F = forcing_CO2_SO2(
            years,
            co2,
            so2;
            f1_CO2 = f1c,
            f3_CO2 = f3c,
            C0_CO2 = co2_C0_ppm,
            f1_SO2 = f1s,
            f2_SO2 = f2s,
            C0_SO2 = C0s,
        )

        itp_F = DataInterpolations.ConstantInterpolation(vcat(F, F[end]), t_itp)
        p = (λ = λ, γ = γ, C1 = C1, C2 = C2, ε = ε, itp_F = itp_F)

        return remake(prob; u0 = [T10, T20], p = p)
    end

    output_func = function (sol, i)
        if sol.retcode != SciMLBase.ReturnCode.Success
            return (fill(NaN, n_t), false)
        end
        T1 = sol[1, :]
        if any(!isfinite, T1)
            return (fill(NaN, n_t), false)
        end
        return (T1, false)
    end

    reduction = function (u, batch, I)
        for (k, idx) in enumerate(I)
            @inbounds u[:, idx] .= batch[k]
        end
        return (u, false)
    end

    u_init = Matrix{Float64}(undef, n_t, n_ens)

    ensprob = SciMLBase.EnsembleProblem(
        prob0;
        prob_func = prob_func,
        output_func = output_func,
        reduction = reduction,
        u_init = u_init,
        safetycopy = false,
    )

    sim = solve(
        ensprob,
        alg,
        ensemblealg;
        trajectories = n_ens,
        batch_size = batch_size,
        saveat = t_pts,
        dense = false,
        save_everystep = false,
        abstol = abstol,
        reltol = reltol,
    )

    Y_T1 = sim.u

    ECS = zeros(n_ens)
    TCR = zeros(n_ens)
    Threads.@threads for i in 1:n_ens
        x = @view X_ens[:, i]
        _, _, λ, γ, _, _, _, f1c, f3c, _, _, _, _ = unpack_x(x, ps)
        ecs, tcr, _ = compute_ECS_TCR(; λ = λ, γ = γ, f1_CO2 = f1c, f3_CO2 = f3c, C0_CO2 = co2_C0_ppm)
        ECS[i] = ecs
        TCR[i] = tcr
    end

    return Y_T1, ECS, TCR
end

# Computes yearly ensemble quantiles and return them as a table.
function quantiles_over_ensemble(years, Y; probs = [0.05, 0.50, 0.95], prefix = "p")
    out = DataFrame(year = years)
    for p in probs
        q = [
            begin
                vals = collect(filter(!isnan, view(Y, t, :)))
                isempty(vals) ? NaN : quantile(vals, p)
            end
            for t in 1:length(years)
        ]
        col = Symbol(prefix * lpad(string(Int(round(p * 100))), 2, '0'))
        out[!, col] = q
    end
    return out
end

# Writes prior and posterior ensemble temperature trajectories, along with their yearly quantile summaries, to CSV files.
function write_ensemble_csvs(outdir, years, Yprior, Ypost)
    mkpath(outdir)

    CSV.write(
        joinpath(outdir, "prior_temperature_ensemble.csv"),
        DataFrame([Symbol("ens_$i") => Yprior[:, i] for i in 1:size(Yprior, 2)]),
    )
    CSV.write(
        joinpath(outdir, "temperature_quantiles_prior.csv"),
        quantiles_over_ensemble(years, Yprior),
    )
    CSV.write(
        joinpath(outdir, "posterior_temperature_ensemble.csv"),
        DataFrame([Symbol("ens_$i") => Ypost[:, i] for i in 1:size(Ypost, 2)]),
    )
    CSV.write(
        joinpath(outdir, "temperature_quantiles_posterior.csv"),
        quantiles_over_ensemble(years, Ypost),
    )
end

end # module
