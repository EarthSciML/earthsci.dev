# Preparing Your Environment

## Installing Julia

EarthSciML is based on the Julia programming language. To get started, you'll need to install Julia. You can download the latest version of Julia from the [official website](https://julialang.org/downloads/).

After installing, you can use Julia straight from the command line:

```bash
julia
```

which will give you the [Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/) (Read-Evaluate-Print Loop) that you can start entering commands into directly.
However, most users prefer to use some sort of development environment to work with Julia.
Some common options include:
* [Pluto notebooks](https://plutojl.org/), which are an interactive way to work with Julia code in a notebook format. Installation instructions are [here](https://plutojl.org/#install).
* [Jupyter notebooks](https://jupyter.org/), which are a more traditional notebook format that can be used with Julia. Installation instructions are [here](https://julialang.github.io/IJulia.jl/stable/manual/running/).
* [VS Code](https://code.visualstudio.com/), which is a popular code editor that has a Julia extension. Installation instructions are [here](https://www.julia-vscode.org/).

## Installing Software Packages

EarthSciML is a collection of Julia software packages that are designed to work together and to work with other packages in the Julia ecosystem.
To use any of the packages, you first need to install and load them.
How to do this varies slightly depending on the development environment you're using:

### In the Terminal or in VS Code

After starting Julia, for example by typing `julia` in a terminal or by entering `cmd/ctrl + shift + p` in VS Code and the typing `Julia: Start REPL`, you can install a package by entering the package manager mode by pressing the `]` key, then typing `add` followed by the package name, for example:

```julia
] add EarthSciMLBase
```

Then, you can press the backspace key to exit the package manager mode and the load the package by typing (for example):

```julia-repl
> using EarthSciMLBase
```
### In a Pluto Notebook

In Pluto notebooks, you don't need to worry about installing the package, you can just directly load it in one of the notebook cells, for example `using EarthSciMLBase`.

### In a Jupyter Notebook

It is not straightforward to install packages from within Jupyter notebooks.
Instead, you can install the package from the Julia terminal as described above, then restart the Jupyter kernel.

## Keeping Track of Your Work

In any computing project, it is important to have to organize your work and keep track of different versions of different documents.
There are obviously many different ways to do this, but one common way is to use a version control system like [Git](https://git-scm.com/), and potentially a system for backing up and sharing your code such as [GitHub](https://github.com/).
You can interact with Git and GitHub directly within VS Code, as described [here](https://youtu.be/i_23KUAEtUM?feature=shared), or you could use a standalone app such as [GitHub Desktop](https://desktop.github.com/).

## Next Steps

Now that we've got our work environment set up, let's look at how to use [ModelingToolkit](https://mtk.sciml.ai/), which is the Julia library with EarthSciML is built on top of. 
Or, if you prefer, you can skip directly to [Using EarthSciML](@ref).
