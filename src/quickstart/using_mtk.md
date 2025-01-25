# Using ModelingToolkit

[ModelingToolkit](https://mtk.sciml.ai/) is a Julia library that provides a high-level interface for defining and manipulating mathematical models. It is the foundation of EarthSciML, and is used to define the physical and chemical processes that are used in EarthSciML models.
Refer to the [documentation](https://mtk.sciml.ai/stable/) for comprehensive information, but we'll go through a quick tutorial here.

The first step is to load the library (see [Installing Software Packages](@ref) for details):

```@example using_mtk
using ModelingToolkit
```
We'll also be working with units here, which will help us check our work and avoid unit conversion errors. 
This requires the `DynamicQuantities` package:

```@example using_mtk
using DynamicQuantities
```
Finally, we will load a two specific variables from the ModelingToolkit library to make our code more concise, specifically `t` which is time in units of seconds, and `D` which is the differential operator with respect to `t`.

```@example using_mtk
using ModelingToolkit: t, D
```

## Variables, Parameters, and Constants

Using ModelingToolkit, we create models by defining systems of differential or algebraic equations, which in turn are functions of *variables*, *parameters*, and *constants*.
We demonstrate how this works by modeling a chemical reaction.

### Variables

In ModelingToolkit, variables are entities that are affected by the dynamics of the system.
Our reaction will have two variables `A` and `B`, which represent the concentrations of two chemical species.
We define them like this:

```@example using_mtk
@variables A(t) = 1 [unit=u"kg/m^3", description="Concentration of A"]
@variables B(t) = 0 [unit=u"kg/m^3", description="Concentration of B"]
nothing #hide
```

There's a lot going on there, let's unpack it. 
`@variables` is something in Julia called a [macro](https://docs.julialang.org/en/v1/manual/metaprogramming/#man-macros).
Feel free to go to the link to see what a macro is, but for our purposes we just need to know that it specifies that we are defining a variable.
The next part, `A(t)` or `B(t)`, species the variable name and that the variable is a function of time `t`.
After that, the `= 1` or `= 0` specifies the default initial value of the variable, which we can override later if we want.
Finally, `[unit=u"kg/m^3", description="Concentration of A"]` specifies the units of the variable and a description of what the variable represents.

### Parameters

In ModelingToolkit, parameters are entities that are fixed for the duration of the simulation, but can be changed between simulations.
Our model will have one parameter `T`, which represents the temperature of the system:

```@example using_mtk
@parameters T = 293.15 [unit=u"K", description="Temperature"]
nothing #hide
```

`@parameters` works the same as `@variables`, but it specifies that we are defining a parameter instead of a variable.
`= 300` 

### Constants

Finally, constants in ModelingToolkit are values that are not expected to ever change.
Our model will have contants T_0, which is a reference temperature, and `k`, which is a rate constant for the reaction:

```@example using_mtk
@constants k = 0.1 [unit=u"1/s", description="Rate constant"]
@constants T_0 = 300 [unit=u"K", description="Initial temperature"]
nothing #hide
```

Now we're ready to define the equations that govern the dynamics of our system.

## Equation Systems

The equations that govern the dynamics of our system can be defined as a system of ordinary differential equations (ODEs).
We can define them like this:

```@example using_mtk
eqs = [
    D(A) ~ -k * exp(T/T_0) * A,
    D(B) ~ k * exp(T/T_0) * A
]
```

Here, `eqs` is the variable we're creating to store our equations, the `[...]` square brackets enclose the list of equations we're creating (which are separated by a comma), `D(A)` and `D(B)` represent the derivatives of the variables `A` and `B` with respect to time, and `~` is where we would use `=` in a normal equation (because `=` already means something else in Julia).
Overall, we have specified that `A` decreases and `B` increases at the rate of `k * exp(T/T0) * A`, which is just a standard first-order chemical reaction.

## Running an ODE Simulation

Now that we have defined our system of equations, we can run a simulation.
First, we set up an [`ODESystem`](https://docs.sciml.ai/ModelingToolkit/stable/systems/ODESystem/), which is a ModelingToolkit object that represents a system of ODEs.

```@example using_mtk
@named sys = ODESystem(eqs, t, [A, B], [T])
```

Again, there's several parts here, but starting from the left:
`@named sys` is saying the both the name of the model and the name of the variable we're assigning it to is `sys`.
Then, `ODESystem(eqs, t, [A, B], [T])` is creating a new `ODESystem` object, where `eqs` is the list of equations we defined earlier, `t` is the time variable, `[A, B]` is the list of variables, and `[T]` is the list of parameters.

Now, almost we're ready to run our simulation. 
First we need the [`OrdinaryDiffEq`](https://diffeq.sciml.ai/stable/) package, which provides the differential equation solvers we need to solve our ODE system.

```@example using_mtk
using OrdinaryDiffEq
```

Now, we do just a few more things:
first, we run the [`structural_simplify`](https://docs.sciml.ai/ModelingToolkit/dev/basics/Composition/#Structural-Simplify) function on our system, which checks if there are any manipulations that can be made to our equation system to make it easier to solve (there aren't any in this case).
Then, we create an [`ODEProblem`](https://docs.sciml.ai/ModelingToolkit/dev/systems/ODESystem/#SciMLBase.ODEProblem-Tuple{ModelingToolkit.AbstractODESystem,%20Vararg{Any}}) from our ODESystem.
`ODEProblem`s are structures are can directly be solved.
Then, we finally can call the [`solve`](https://diffeq.sciml.ai/stable/basics/common_solver_opts/#Solving-the-Problem) function to get the result of our simulation, which we save in the `sol` variable.
When we call the `solve` function, we include the `tspan=(0,10)` argument, which specifies that we want to run the simulation with the start time of 0 seconds and the end time of 10 seconds.

```@example using_mtk
sys_simplified = structural_simplify(sys)
prob = ODEProblem(sys_simplified)
sol = solve(prob, Tsit5(), tspan=(0,10))
```

The result of the simulation is shown above.
Note that it starts with `retcode: Success` which means that the simulation was successful.
If `retcode` is something else, the solution make include some results, but they will not be correct.

## Plotting the Results

There are many software packages in Julia that can be used to make plots.
Here, we're going to use one called [Plots.jl](https://docs.juliaplots.org/stable/):

```@example using_mtk
using Plots

plot(sol, xlabel="time (s)", ylabel="Concentration (kg/m³)")
```

The plot above shows that `A` starts out at a concentration of 1 and decreases over time, while `B` starts out at a concentration of 0 and increases over time at the same rate.

We can also re-run the simulation with different values of the parameters or different initial values for the variables, and compare the result.
For example, in the code below, we increase the initial concentration of `A` to 1.5 and decrease the temperature to 100 K:

```@example using_mtk
prob2 = ODEProblem(sys_simplified, [A=>1.5], (0, 10), [T=>100])
sol2 = solve(prob2, Tsit5())

plot(
    plot(sol, xlabel="time (s)", ylabel="Concentration (kg/m³)", 
        ylim=(0, 1.5), title="T=293.15"),
    plot(sol2, title="T=100", xlabel="", ylabel="", legend=:none),
    size=(1000, 400)
)
```

## Partial Differential Equations

We can also simulate partial differential equations (PDEs) using ModelingToolkit.
For example, we can consider the same reaction system as above, but occuring in a fluid that is moving.
To do this, let's first define the boundaries of our spatial system, where we have x ∈ {0, 1}, y ∈ {0, 1}, and t ∈ {0, 10} as before.

```@example using_mtk
x_min = y_min = t_min = 0.0
x_max = y_max = 1.0
t_max = 10.0
```
Let's split up our domain into 32 grid cells in each direction:

```@example using_mtk
N = 32
```

and calculate the size of each grid cell:

```@example using_mtk
dx = (x_max-x_min)/N
dy = (y_max-y_min)/N
```

We can now create a function which injects emissions into the domain at a specific location:

```@example using_mtk
islocation(x, y) = x > x_max / 2 - dx && x < x_max / 2 + dx && 
    y > y_max / 10 - dx && 
    y < y_max / 10 + dx
emission(x, y) = ifelse(islocation(x, y), 10, 0)
@register_symbolic emission(x, y)
```

(Because our `emission` function includes operations which are not allowed in ModelingToolkit, we have to register it as a symbolic function. See [here](https://docs.sciml.ai/ModelingToolkit/dev/basics/FAQ/#ERROR:-TypeError:-non-boolean-(Num)-used-in-boolean-context?) for more information.)

We'll also add function to specify the movement of the fluid in the domain, which we will specify as traveling in a circle:

```@example using_mtk
θ(x,y) = atan(y.-0.5, x.-0.5)
u(x,y) = -sin(θ(x,y))
v(x,y) = cos(θ(x,y))
```

Now we can re-create our equation system from above, but this time we will add terms to account for the advection of chemicals `A` and `B` in the host fluid.
Because units and partial differential equations don't currently work together in ModelingToolkit, we will also recreate our variables, parameters, and constants without units.

```@example using_mtk
@parameters T=293.15 x y
t_ = ModelingToolkit.t_nounits # We want the time variable without units now.
D_ = ModelingToolkit.D_nounits # We want the differential operator without units now.
@constants k=0.1 T_0=300.0
@variables A(..) B(..)
Dx = Differential(x)
Dy = Differential(y)
advect(var) = -u(x,y)*Dx(var) - v(x,y)*Dy(var)
eqs = [
    D_(A(x,y,t_)) ~ advect(A(x, y, t_)) + emission(x, y) - k*exp(T/T_0)*A(x,y,t_),
    D_(B(x,y,t_)) ~ advect(B(x, y, t_)) + k*exp(T/T_0)*A(x,y,t_),
]
```

As you can see, we need to define the variables and equations slightly differently for PDEs. (See [here](https://docs.sciml.ai/MethodOfLines/stable/) for more information.)

Next, we need to more formally specify our spatial and temporal domain, which requires the `DomainSets` package:

```@example using_mtk
using DomainSets
domain = [
    x ∈ Interval(x_min, x_max),
    y ∈ Interval(y_min, y_max),
    t ∈ Interval(t_min, t_max),
]
```

We also need to specify our boundary conditions:

```@example using_mtk
bcs = [A(x,y,t_min) ~ 0.0,
       A(x_min,y,t_) ~ A(x_max,y,t_),
       A(x,y_min,t_) ~ A(x,y_max,t_),

       B(x,y,t_min) ~ 0.0,
       B(x_min,y,t_) ~ B(x_max,y,t_),
       B(x,y_min,t_) ~ B(x,y_max,t_),
] 
```

Now we can create our `PDESystem`, which is similar to an `ODESystem` but for PDEs.

```@example using_mtk
@named pdesys = PDESystem(eqs, bcs, domain, [x,y,t], [A(x,y,t), B(x,y,t)], [T])
```

One additional difference with PDE systems is that we need to [discretize](https://docs.sciml.ai/MethodOfLines/stable/api/discretization/) them to convert them to an ODE system that we can solve. 
To do this, we need to use the [MethodOfLines](https://docs.sciml.ai/MethodOfLines/stable/) package:

```@example using_mtk
using MethodOfLines
discretization = MOLFiniteDifference([x=>dx, y=>dy], t)
prob = discretize(pdesys,discretization)
```

Then, we can finally run our simulation.
This time, we specify a different ODE solver called [`TRBDF2`](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/#SDIRK-Methods).
We also specify that we want to save the results at every 0.1 seconds:

```@example using_mtk
sol = solve(prob, TRBDF2(), saveat=0.1)
```

Finally, we can plot the results of our simulation:

```@example using_mtk
disc_t = sol[t]
disc_x = sol[x]
disc_y = sol[y]
solA = sol[A(x, y, t)]
solB = sol[B(x, y, t)]

anim = @animate for i in eachindex(disc_t)
    plot(
        heatmap(disc_x, disc_y, solA[:,:,i], xlabel="x", ylabel="y", title="A Concentration",
            clim=(0, maximum(solA))),
        heatmap(disc_x, disc_y, solB[:,:,i], title="B Concentration",
            clim=(0, maximum(solB))),
        size=(1000, 400)
    )
end
gif(anim, fps=8)
```
!!! warning
    This method for simulating PDEs currently doesn't work for large-scale simulations. EarthSciML currently uses a different method for representing spatial operators, which we will harmonize with the method shown here in the future.

## Isn't There an Easier Way?

As you can see, we can use ModelingToolkit to define and simulate complex systems of equations.
However, you can probably also see that as the system we're modeling gets more complex, it gets more complex to keep track of things and make sure everything is specified correctly.
That's where EarthSciML comes in.

EarthSciML is a collection of Julia packages that are designed to simplify and automate the process of specifying geoscientific models using ModelingToolkit.
Read on to see how it works!