# Bifurcation Problem

The idea behind `BifurcationKit` is to compute bifurcation diagrams in memory limited environments where the device can barely hold the current continuation state. We thus disable by default saving all solutions along the branch and all eigenvectors (see [`ContinuationPar`](@ref) to change this behaviour). Still, one needs to save a few solution indicators, like for plotting. This is the reason for the function `recordFromSolution` (see below).

## Bifurcation Problem

[`BifurcationProblem`](@ref) is the basic structure for a bifurcation problem which holds the following fields:

- the vector field
- an initial guess
- a set of parameters
- a parameter axis

as well as user defined functions for 

- plotting, `plotSolution`
- recording (`recordFromSolution`) indicators about the solution when this one is too large to be saved at every continuation step.

### Example

```julia
f(x,p) = @. sin(x * p.a)
u0 = zeros(100_000_000) 
params = (a = 1.0, b = 2.0)

# record a few components / indicators about x 
myRecord(x,p) = (x1 = x[1], max = maximum(x), nrm = norm(x, Inf))

prob = BifurcationProblem(f, u0, p, (@lens _.a);
	recordSolution = myRecord
	)
```


## Problem modification

In case you want to modify an existing problem, you should use the following method

```@docs
reMake(prob::BifurcationKit.AbstractBifurcationProblem;
		u0 = prob.u0,
		params = prob.params,
		lens::Lens = prob.lens,
		recordFromSolution = prob.recordFromSolution,
		plotSolution = prob.plotSolution,
        J = missing,
        d2F = missing,
		d3F = missing)
```

## Example

```@example 
using BifurcationKit, Setfield
F(x,p) = @. p.a + x^2
# parameters
par = (a = 0., b = 2)
prob = BifurcationProblem(F, zeros(3), par, (@lens _.a))
```

