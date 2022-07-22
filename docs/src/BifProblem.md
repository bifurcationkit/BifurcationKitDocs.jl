# Bifurcation Problem

## Vector field

The structure [`BifFunction`](@ref) holds the vector field. If you just pass a function, everything (jacobian, second derivative, ...) will be evaluated using automatic differentiation.

## Bifurcation Problem

The structure [`BifurcationProblem`](@ref) is the basic holder for a bifurcation problem which holds 

- the vector field
- an initial guess
- a set of parameters
- a parameter axis

as well as user defined functions for 

- plotting
- recording a indicators about the solution when this oe is too large to save at every step

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

