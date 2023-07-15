@def title = "Earth Science Machine Learning"
@def subtitle = "A next-generation software framework for geoscientific modeling and analysis"
@def tags = ["syntax", "code"]

A next-generation software framework for geoscientific modeling and analysis

<!-- \tableofcontents <!-- you can use \toc as well -->

## Available libraries

Earth Science models are important tools for understanding and making predictions about our environment.
The next generation of geoscientific models needs to be integrated across disciplines and to better facilitate convergent research.
However, current model development focuses on tightly-coupled model components for numerical efficiency, requiring development or gatekeeping by a close-knit team.
These two requirements are mutually exclusive and incompatible: development processes that require contributions from diverse individuals require components that are modular and granular, whereas the development of tightly coupled components requires a centrally-managed organization.
We term this challenge the *community fragmentation problem*.

An additional current barrier to the openness of geoscientific resources is that geoscientists developing algorithms performant enough for inclusion in large-scale models are also required to have extensive knowledge of numerical methods for integrating differential equations, and these skills are not taught in many geoscientific curricula.
We term this the *two-domain problem*.

This project aims to solve these two problems by seeding a transition to symbolic equation-based model development, where model components and their interrelationships are specified symbolically as systems of differential-algebraic equations, visually similar to how they would appear in a scientific journal article.

* [EarthSciMLBase.jl](https://base.earthsci.dev)
* [EarthSciData.jl](https://data.earthsci.dev)
* [GasChem.jl](https://gaschem.earthsci.dev)
* [Aerosol.jl](https://aerosol.earthsci.dev)
