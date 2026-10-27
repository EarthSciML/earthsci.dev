---
title: "Quick start"
description: "Install an EarthSci library and load, check and simulate your first model."
---

This page takes you from nothing to a simulated model. It uses
[EarthSciAST](/ast/), the core of the ecosystem, because every other project
builds on its `.esm` documents.

## 1. Install

Pick the language you work in. Each of these installs EarthSciAST.

**Julia**

```julia
using Pkg
Pkg.add("EarthSciAST")
```

**Python**

```bash
pip install earthsci-ast
```

**TypeScript / Node.js**

```bash
npm install @earthsciml/ast
```

**Rust**

```toml
[dependencies]
earthsci-ast = "0.1.1"
```

**Go**

```bash
go get github.com/EarthSciML/EarthSciAST/pkg/earthsci-ast-go
```

Simulation is available in the Julia, Python and Rust bindings; TypeScript and
Go read, check and transform documents but do not simulate. See the
[capability table](/ast/) for details. The Python simulation example below may
need the package's `simulate` extra (`pip install "earthsci-ast[simulate]"`).

## 2. A first model

A model is an `.esm` file. This one, `decay.esm`, is a single quantity `N`
that decays at rate `lambda`, with its half-life as an observed value:

```json
{
  "esm": "1.0.0",
  "metadata": { "name": "Decay", "authors": ["EarthSciML"] },
  "models": {
    "Decay": {
      "variables": {
        "N": { "type": "unknown", "units": "mol", "default": 100.0 },
        "lambda": { "type": "parameter", "units": "1/s", "default": 0.1 },
        "half_life": { "type": "unknown", "units": "s" }
      },
      "equations": [
        { "lhs": { "op": "D", "args": ["N"], "wrt": "t" },
          "rhs": { "op": "*", "args": [{ "op": "-", "args": ["lambda"] }, "N"] } },
        { "lhs": "half_life",
          "rhs": { "op": "/", "args": [{ "op": "log", "args": [2] }, "lambda"] } }
      ]
    }
  }
}
```

You rarely need to write this JSON by hand; the point is that it is plain data
every binding understands.

## 3. Load, validate, simulate (Python)

```python
import earthsci_ast

doc = earthsci_ast.load_path("decay.esm")      # parse + schema check
assert list(doc.models) == ["Decay"]

result = earthsci_ast.validate(doc)            # structural rules and units
assert result.is_valid

prob = earthsci_ast.esm_problem(doc, (0.0, 10.0))
sol = earthsci_ast.solve(prob, saveat=[0.0, 5.0, 10.0])
print(sol["N"])                                # amount remaining at each time
```

The same steps in every language, with runnable and CI-checked examples, are in
the [EarthSciAST getting-started guide](/ast/generated/api/guides/getting-started/).
The generated [API reference](/ast/generated/api/) lists every operation with a
tab per language.

## 4. Where next

| I want to… | Go to |
|---|---|
| Understand the `.esm` format | [EarthSciAST](/ast/) |
| Use a ready-made model component | [EarthSciModels](/models/) |
| Discretize a PDE onto a grid | [EarthSciDiscretizations](/discretizations/) |
| Load real data (ERA5 and others) into a simulation | [EarthSciIO](/io/) |
| Get an analytical Jacobian for a stiff solver | [EarthSciASTDiff](/astdiff/) |
| Run large coupled models | [EarthSciIntegrator](/integrator/) |
