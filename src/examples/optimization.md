# Analysis and Optimization

There are many types of analysis and optimization that one may want to perform using a geoscientific model, for example calculating the relationship between an input and an output, or finding the value of an input or parameter that is consistent with a measured observation.
A main goal of EarthSciML is to streamline and support these types of analysis.

## Gradient Calculation

Many types of analysis and optimization require knowledge of the relationship between one or more of a model's outputs and one or more of its inputs or parameters.
In other words, we often want to know if we change a model input or parameter, how will the model output change, or what change in the inputs or parameters is necessary to result in a desired change in the output.
These types of relationships are referred to in mathematics as [derivatives](https://en.wikipedia.org/wiki/Derivative), [gradients](https://en.wikipedia.org/wiki/Gradient), or [Jacobians](https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant), depending on the number of inputs and outputs involved.
In this example, we are going to calculate the gradient of an model's average output air pollution concentration over time with respect to parameters in the model, which will tell us how we should change those parameters to change the model output in a certain way.

To start, we will define an air quality model similar to the one described in [Using EarthSciML](@ref):

```@example optimization
using EarthSciMLBase, GasChem, AtmosphericDeposition, EarthSciData
using EnvironmentalTransport, ModelingToolkit, OrdinaryDiffEq
using DiffEqCallbacks
using ForwardDiff
using SymbolicIndexingInterface
using ModelingToolkit: t
using Dates, Plots, NCDatasets, Statistics, DynamicQuantities
using ProgressLogging # Needed for progress bar. Use `TerminalLoggers` if in a terminal.
using LinearSolve
using DiffEqCallbacks

domain = DomainInfo(
    DateTime(2016, 5, 1),
    DateTime(2016, 5, 1, 1);
    lonrange = deg2rad(-115):deg2rad(2.5):deg2rad(-68.75),
    latrange = deg2rad(25):deg2rad(2):deg2rad(53.7),
    levrange = 1:2,
    dtype = Float64)

geosfp = GEOSFP("0.5x0.625_NA", domain; stream=false)
geosfp = EarthSciMLBase.copy_with_change(geosfp, discrete_events=[]) # Workaround for bug.

emis = NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", domain; stream=false)
emis = EarthSciMLBase.copy_with_change(emis, discrete_events=[]) # Workaround for bug.

model_base = couple(
    SuperFast(),
    FastJX(),
    #DrydepositionG(), Not currently working
    Wetdeposition(),
    AdvectionOperator(NaN, upwind1_stencil, ZeroGradBC()),
    #emis, # Not currently working with automatic differentiation
    geosfp,
    domain
)
```

In this example, we want to understand how we can "nudge" or adjust each of the variables in the model while it is running to achieve a desired output.
So the first thing that we need to do is to create a new model component that includes the nudging factors.
It gets slightly complicated, but what we want to do is to add a new term to each differential equation in our system that adjusts the rate of change of the variable in question by the current value of that variable times a "nudge factor".
For more information about the specifics of what we're doing here, you can refer to the documentation about [creating and composing model components](https://base.earthsci.dev/dev/composition/) for more information.

```@example optimization
struct NudgeCoupler
    sys
end
function Nudge(; name=:nudge)
    vars = []
    params = []
    for i in 1:13
        n = Symbol(:nudge_, i)
        push!(vars, only(@variables $n(t) = 0.0 [unit = u"s^-1", description = "Nudge for species $i"]))
        cn = Symbol(:nudge_c, i)
        push!(params, only(@parameters $(cn) = 0.0 [unit = u"s^-1", description = "Nudge constant for species $i"]))
    end
    eqs = vars .~ params
    ODESystem(eqs, t, [vars...], [params...]; name=name, metadata=Dict(:coupletype => NudgeCoupler))
end

function EarthSciMLBase.couple2(c::GasChem.SuperFastCoupler, n::NudgeCoupler)
    c, n = c.sys, n.sys
    operator_compose(c, n, Dict(
        c.O3 => n.nudge_1 => c.O3,
        c.OH => n.nudge_2 => c.OH,
        c.HO2 => n.nudge_3 => c.HO2,
        c.H2O => n.nudge_4 => c.H2O,
        c.NO => n.nudge_5 => c.NO,
        c.NO2 => n.nudge_6 => c.NO2,
        c.CH3O2 => n.nudge_7 => c.CH3O2,
        c.CH2O => n.nudge_8 => c.CH2O,
        c.CO => n.nudge_9 => c.CO,
        c.CH3OOH => n.nudge_10 => c.CH3OOH,
        c.ISOP => n.nudge_11 => c.ISOP,
        c.H2O2 => n.nudge_12 => c.H2O2,
        c.HNO3 => n.nudge_13 => c.HNO3,
    ))
end

model = couple(model_base, Nudge())
```

Once we've created our nudging model component and coupled it into our base model, we're almost ready to calculate our gradient.
Before we do that, though, we need to get a few preliminaries out of the way, including extracting the nudging parameters from our model so that we can use them later and figuring out the location in the model results where the NO2 concentration is going to be stored, as that is the output that we're interested in adjusting.

```@example optimization
model_sys = convert(ODESystem, model)
nudge_params = parameters(model_sys)[[only(findall((x)->x==Symbol(:nudge₊nudge_c, i), Symbol.(parameters(model_sys)))) for i in 1:13]]

usize = size(EarthSciMLBase.init_u(model_sys, domain))
iNO2 = only(findall((x) -> x==Symbol("SuperFast₊NO2(t)"), Symbol.(unknowns(model_sys))))

st = SolverIMEX(stiff_sparse=false)
prob = ODEProblem{false}(model, st, callback=PositiveDomain(save=false))
```
The last thing that we need to set up is an objective function: what do we want to calculate the gradient with respect to?
In this case, let's say that we know that our model should always output an NO2 concentration of 42 ppb, so we're interested in finding nudging factors that minimize the difference between the model's predicted NO2 concentration and 42 ppb at all times.

To operationalize this goal, we create a function that we'll call `loss` that takes a vector of nudging factors as input, runs the model with those nudging factors, and calculates the mean squared difference between the model's predicted NO2 concentration and 42 ppb.

```@example optimization
function loss(nudge_vals)
    the_answer = 42.0
    new_params = remake_buffer(model_sys, prob.p, Dict(nudge_params .=> nudge_vals))
    newprob = remake(prob, p=new_params)
    sol = solve(newprob, KenCarp5(linsolve=LUFactorization()); progress=true, progress_steps=1, saveat=3600)
    mean([mean((reshape(ui, usize...)[iNO2] .- the_answer).^2) for ui in sol.u])
end
```
Since there are 13 chemical species in our model, we have 13 nudging factors, so we can run the `loss` function with a vector of 13 zeros to see what the loss is with no nudging factors applied.

```@example optimization
loss(zeros(13))
```

We can also see how our loss value changes if we nudge the derivatives of all the model variables down by 10%:

```@example optimization
loss(zeros(13) .- 0.1)
```

As you can see, the loss is higher in this new simulation, which means that our change made things worse!
So what values of the nudging factors should we use to make the loss as small as possible?

The first step to get that answer is to calculate the gradient of the loss function with respect to the nudging factors, so we know which direction to change each nudging factor to make the loss smaller.

We can use automatic differentiation to efficiently calculate the gradient:

```@example optimization
grad = ForwardDiff.gradient(loss, zeros(13))

bar(["O3", "OH", "HO2", "H2O", "NO", "NO2", "CH3O2", "CH2O", "CO", 
    "CH3OOH", "ISOP", "H2O2", "HNO3"], grad, permute=(:x, :y), size=(400, 250),
    label=:none, xlabel="Species", ylabel="Nudging Sensitivity")
```

Now that's a start! We can see that the error in NO2 concentrations is most sensitive to adjustments to NO2 dynamics (unsurprisingly), but also sensitive to adjustments in CH3O2, NO, and O3.

That's all we have for now, but we will add to this example as we build out the EarthSciML optimization and analysis capabilities.
Check back soon!