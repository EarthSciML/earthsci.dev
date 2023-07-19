# This file was generated, do not modify it. # hide
using EarthSciData, EarthSciMLBase, GasChem,
    DomainSets, ModelingToolkit, MethodOfLines, 
    DifferentialEquations, Dates, Distributions,
    Latexify, Plots, SciMLBase,
    CairoMakie, GraphMakie, MetaGraphsNext

@parameters t lev lon lat
model = SuperFast(t)
ls = latexify(model.rxn_sys) # TODO: Change to model.rxn_sys
print(replace(ls.s, raw"\require{mhchem}" => "")) # hide