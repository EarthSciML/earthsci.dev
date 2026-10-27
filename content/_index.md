---
title: "EarthSci"
description: "Documentation for the EarthSciML software ecosystem: a language-neutral format for Earth-science models and the tools around it."
---

# EarthSci

EarthSci is a family of libraries for describing, discretizing, feeding data
to, and running Earth-science models. They share one idea: a model is a
**document**, not a program. Models are written once in the `.esm` format (an
abstract syntax tree serialized as JSON), and the same document can then be
read, checked, transformed and simulated from Julia, Python, Rust, TypeScript
or Go.

New here? Start with the **[quick start](quickstart/)**.

## The projects

<ul class="project-cards">
<li><strong><a href="/ast/">EarthSciAST</a></strong><br>
The core. The specification of the <code>.esm</code> format, plus bindings in
Julia, TypeScript, Python, Rust and Go to parse, validate, display, edit,
flatten and (in Julia, Python and Rust) simulate documents. Includes the
cross-language API reference.</li>
<li><strong><a href="/discretizations/">EarthSciDiscretizations</a></strong><br>
The standard library of discretization rules, written as <code>.esm</code>
documents: grids, finite-difference and finite-volume rewrite rules,
regridding and reprojection, with conformance goldens.</li>
<li><strong><a href="/models/">EarthSciModels</a></strong><br>
A browsable catalog of <code>.esm</code> model components (chemistry,
physics, and more), each with inline tests.</li>
<li><strong><a href="/io/">EarthSciIO</a></strong><br>
The data-provider library: resolves a data loader's URL, downloads it into a
shared content-addressed cache, and reads it into native-grid arrays. Python,
Julia and Rust.</li>
<li><strong><a href="/astdiff/">EarthSciASTDiff</a></strong><br>
Analytical sparse Jacobians for <code>.esm</code> models, by symbolic
differentiation of the expression AST.</li>
<li><strong><a href="/integrator/">EarthSciIntegrator</a></strong><br>
A specified-once, natively implemented (Julia and Rust) integrator for full 3-D
<code>.esm</code> models, including partitioning and discrete adjoints.
Currently at the specification stage.</li>
</ul>

## How they fit together

- **EarthSciAST** defines the document format and the tools that read and
  evaluate it. Everything else builds on it.
- **EarthSciDiscretizations** and **EarthSciModels** are *data*: libraries of
  `.esm` documents that models import or reference.
- **EarthSciIO** is the one place that touches the network and the file system:
  it supplies data to a running simulation but never solves anything.
  Variable remapping stays in EarthSciAST and regridding stays in
  EarthSciDiscretizations.
- **EarthSciASTDiff** and **EarthSciIntegrator** are numerical machinery on top
  of the AST: Jacobians for stiff solvers, and a planner and integrator for large
  coupled models.

This site is assembled from the documentation of the six repositories under
the [EarthSciML](https://github.com/EarthSciML) organization; it is rebuilt
whenever any of them changes on `main`.
