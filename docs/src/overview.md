# Overview of BifurcationKit.jl

```@contents
Pages = ["overview.md"]
Depth = 3
```

The general workflow for using the package is as follows:

- Define a problem
- Solve the problem
- Analyze the output

## Defining Problems

Problems are specified via a type interface. The problem types are designed to contain the necessary information to fully define their associated continuation method. For example, a bifurcation problem is defined by

```math
F(u, pars) = 0
```

with some parameters `pars`, some initial guess `u0`, and scalar parameter axis `lens` contained in `pars`. Therefore the `BifurcationProblem` is defined by those components:

```julia
prob = BifurcationProblem(F, u0, pars, lens)
```

Note that the number types in the solution will match the types you designate
in the problem. However complex types are not advised as they mess up the detection of bifurcation points.

## Continuing from the initial guess

Each type of bifurcation problem has its own problem type which allow the solvers
to dispatch to the right methods. The common interface for calling the solvers is:

```julia
br = continuation(prob, alg; kwargs)
```

Into the command, one passes the bifurcation problem that was defined
`prob`, choose an algorithm `alg`, and change the properties of the solver using keyword arguments. 
The solver returns a branch object `br` which holds all the details for the branch.

## Analyzing the branch

The solution type has a common interface, which makes handling the solution similar between the
different types of bifurcation problems. Tools such as interpolations
are seamlessly built into the solution interface to make analysis easy. This
interface is described in the [`ContResult`](@ref).

Plotting functionality is provided by a recipe to `Plots.jl`. To
use plot branches, simply call the `plot(br)` and the plotter will generate
appropriate plots. Plots can be customized using all the keyword arguments
provided by Plots.jl. Please see Plots.jl's documentation for more information.
