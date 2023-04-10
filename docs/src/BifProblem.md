# Bifurcation Problem

## Vector field

The structure [`BifFunction`](@ref) holds the vector field. If you just pass a function, everything (jacobian, second derivative, ...) will be evaluated using automatic differentiation.

## Bifurcation Problem

[`BifurcationProblem`](@ref) is the basic structure for a bifurcation problem which holds the following fields. The link [`BifurcationProblem`](@ref) provides more information.

- the vector field
- an initial guess
- a set of parameters
- a parameter axis

as well as user defined functions for 

- plotting, `plotSolution`
- recording (`recordFromSolution`) indicators about the solution when this one is too large to save at every step.

> `BifurcationKit` is designed to be used in very limited memory environments (GPU) where you can barely hold the equilibrium and some of the eigen-elements of the jacobian. Thus, `BifurcationKit` does not record the equilibria unless stated, *i.e.* when `ContinuationPar.saveSolEveryStep > 0`

### Example

```julia
f(x,p) = @. sin(x * p.a)
u0 = zeros(100_000_000) 
params = (a = 1.0, b = 2.0)

# save a few components / indicators of x because it is too costly 
# to keep all solutions in memory
function myRecord(x,p)
	return (x1 = x[1], max = maximum(x), nrm = norm(x, Inf))
end 

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

