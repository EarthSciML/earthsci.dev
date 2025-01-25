# The Earth Science Machine Learning (EarthSciML) Framework

EarthSciML is a next-generation software framework for geoscientific modeling and analysis.
It allows geoscientific models to be specified as systems of equations which can be coupled together and compiled to a runnable model.
In addition to simplifying the process of creating and using geoscientific models, EarthSciML also (will soon) enable the use of GPU-accelerated simulations as well as automatic creation of advanced analysis tools such as model adjoints.

### Simulation

A main goal of EarthSciML is to allow users to easily set up and explore dynamics of geoscientific models.
Below, we demonstrate this by coupling a model of gas-phase atmospheric chemistry with a model of photolysis and then demonstrating how the resulting dyanamics are affected by temperature:

```@example index
using EarthSciMLBase, GasChem
using ModelingToolkit, OrdinaryDiffEq, Plots, SymbolicIndexingInterface

model = convert(ODESystem, couple(SuperFast(), FastJX()))

prob = ODEProblem(model, (), (0,3*24*3600))
T_setter = setp(prob, [model.SuperFast₊T, model.FastJX₊T])

# Create an animation of the O₃ concentration as a function of temperature.
anim = @animate for T in 100:10:400
    T_setter(prob, (T, T))
    sol = solve(prob, Rosenbrock23())
    plot(sol.t, sol[model.SuperFast₊O3], xlabel="Time (s)", ylim=(0, 50),
        ylabel="O₃ Concentration (ppb)", label=:none, title="T=$T K")
end
gif(anim, fps = 10)
```

### Analysis

A second main goal of EarthSciML is to make it easy to perform advanced analysis of model dynamics, particularly the type of analysis that can streamline the integration of geoscientific modeling and machine learning. Below, we demonstrate this by using gradient descent with [automatic differentiation](https://juliadiff.org/) to find a temperature that results in an average ozone concentration of 25 ppb.

```@example index
using Optimization, OptimizationOptimJL
using Statistics

# Set up functions to run the simulation and calculate the error.
function run_simulation(u)
    T = u[1]
    p2 = remake_buffer(model, prob.p, Dict(model.SuperFast₊T => T, model.FastJX₊T => T))
    solve(remake(prob, p=p2), Rosenbrock23())
end
loss(u, p) = (mean(run_simulation(u)[model.SuperFast₊O3]) - 25)^2

optfn = OptimizationFunction(loss, Optimization.AutoForwardDiff())
optprob = OptimizationProblem(optfn, [320.0], ())

# This is some code to save the results for making an animation.
Ts = []
sols = []
function callback(state, loss_val) 
    push!(Ts, state.u[1])
    push!(sols, run_simulation(state.u))
    false
end

# This line actually does the optimization.
sol = solve(optprob, BFGS(), callback=callback)

# The rest of the code makes the animation.
anim = @animate for i in eachindex(Ts)
    o3avg = round(mean(sols[i][model.SuperFast₊O3]), digits=3)
    plot(sols[i].t, sols[i][model.SuperFast₊O3], label=:none, ylim=(0, 35),
        xlabel="Time (s)", ylabel="O₃ Concentration (ppb)",
        title="Step=$i; T=$(round(Ts[i], digits=3)) K; O₃ avg.=$(o3avg) ppb", )
    plot!([sols[i].t[begin], sols[i].t[end]], [o3avg, o3avg], label=:none,
        linecolor=:black, linestyle=:dash)

end
gif(anim, fps = 5)
```


## Where to Start

* If you'd like to know more about the available capabilities, check out our [Examples](@ref) section.
* If you're interested in jumping right in and using the framework on your own computer, check out the [Getting Started](@ref) section.
* If you're interested in learning about the theory behind geoscientific modeling and machine learning, check out [EarthSci EDU](@ref).
* Finally, if you're interested in contributing to the project, check [Contributing](@ref).


## Reproducibility
```@raw html
<details><summary>The documentation of this EarthSciML package was built using these direct dependencies,</summary>
```
```@example
using Pkg # hide
Pkg.status() # hide
```
```@raw html
</details>
```
```@raw html
<details><summary>and using this machine and Julia version.</summary>
```
```@example
using InteractiveUtils # hide
versioninfo() # hide
```
```@raw html
</details>
```
```@raw html
<details><summary>A more complete overview of all dependencies and their versions is also provided.</summary>
```
```@example
using Pkg # hide
Pkg.status(;mode = PKGMODE_MANIFEST) # hide
```
```@raw html
</details>
```
```@raw html
You can also download the 
<a href="
```
```@eval
using TOML
using Markdown
version = TOML.parse(read("../Project.toml",String))["version"]
name = TOML.parse(read("../Project.toml",String))["name"]
link = Markdown.MD("https://github.com/EarthSciML/"*name*"/tree/gh-pages/v"*version*"/assets/Manifest.toml")
```
```@raw html
">manifest</a> file and the
<a href="
```
```@eval
using TOML
using Markdown
version = TOML.parse(read("../Project.toml",String))["version"]
name = TOML.parse(read("../Project.toml",String))["name"]
link = Markdown.MD("https://github.com/EarthSciML/"*name*"/tree/gh-pages/v"*version*"/assets/Project.toml")
```
```@raw html
">project</a> file.
```
