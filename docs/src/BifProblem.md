# Bifurcation problems

```@contents
Pages = ["BifProblem.md"]
Depth = 3
```

The idea behind `BifurcationKit` is to compute bifurcation diagrams in memory limited environments where the device can barely hold the current continuation state. We thus disable by default saving all solutions along the branch and all eigenvectors (see [`ContinuationPar`](@ref) to change this behavior). Still, one needs to save a few solution indicators, like for plotting. This is the reason for the function `record_from_solution` (see below).

## Generic bifurcation problem

[`BifurcationProblem`](@ref) is the basic / generic structure for encoding a bifurcation problem ; it holds the following fields:

- the vector field
- an initial guess
- a set of parameters
- a parameter axis

as well as user defined functions for 

- plotting, `plot_solution`
- recording (`record_from_solution`) indicators about the solution when this one is too large to be saved at every continuation step.

### Example

```julia
f(x,p) = @. sin(x * p.a)
u0 = zeros(100_000_000) 
params = (a = 1.0, b = 2.0)

# record a few components / indicators about x 
myRecord(x,p;k...) = (x1 = x[1], max = maximum(x), nrm = norm(x, Inf))

prob = BifurcationProblem(f, u0, p, (@optic _.a);
	record_from_solution = myRecord
	)
```


## Problem modification

In case you want to modify an existing problem, you should use the following method

```@docs
re_make(prob::BifurcationKit.AbstractBifurcationProblem;
		u0 = prob.u0,
		params = prob.params,
		lens = prob.lens,
		record_from_solution = prob.record_from_solution,
		plot_solution = prob.plot_solution,
       J = missing,
       d2F = missing,
		d3F = missing)
```

### Example

```@example 
using BifurcationKit, Setfield
F(x,p) = @. p.a + x^2
# parameters
par = (a = 0., b = 2)
prob = BifurcationProblem(F, zeros(3), par, (@optic _.a))
# change u0
prob2 = BifurcationKit.re_make(prob, u0 = rand(3))
```

