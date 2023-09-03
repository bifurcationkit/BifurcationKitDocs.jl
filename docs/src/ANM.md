#  Asymptotic numerical method (ANM)

!!! warning "Work in progress"
    Automatic branch switching is being tested, it will be available soon.

!!! warning "Dimensions"
    This is a method for small dimensions, less than several thousands.

> To access this algorithm, you have to use the package [AsymptoticNumericalMethod.jl](https://github.com/bifurcationkit/AsymptoticNumericalMethod.jl)   

The method [^Rubbert],[^Charpentier] seeks a Taylor approximation of the solutions of

$$\mathbf F(X,p)=0\in\mathbb R^n$$

where $X\in\mathbb R^n, p\in\mathbb R$. The solution is found by solving

$$F(x(s),p(s))= 0$$

$$\langle x(s)-x(0),\frac{\partial u}{\partial s}(0)\rangle + (p(s)-p(0))\frac{\partial p}{\partial s}(0) = s$$

where

$$x(s)=\sum\limits_{k=0}^K x_ks^k,\quad p(s)=\sum\limits_{k=0}^K p_ks^k$$

for some user passed $K>0$. It gives

$$F(x(s),p(s)) = \sum\limits_{k=0}^K F_ks^k+h.o.t.$$

from which we deduce the equations $F_k=0$. We then find:

$$F_{1, x_1=Id, p_1=1}(x_k,p_k) = -F_{k,x_k=0,p_k=0}.$$

The validity range of the solution can be estimated by

$$r_K = \left(\frac{\epsilon}{\lvert\lvert F_K\lvert\lvert}\right)^{1/K}$$

This allows to iterate and find a sequence of series which spans the parameter range.

## Implementation

The method is based on the package [TaylorSeries.jl](https://github.com/JuliaDiff/TaylorSeries.jl) which makes it easy to manipulate Taylor series based on Automatic Differentiation.

## Method

See [`AsymptoticNumericalMethod.ANM`](@ref) for more information.

## Example

We provide an example of use. We define a `BifurcationProblem` as usual and pass the continuation algorithm `ANM`.

```@example ANM
using AsymptoticNumericalMethod, Plots, Parameters, Setfield
using LinearAlgebra: norm
using BifurcationKit
const BK = BifurcationKit

norminf(x) = norm(x, Inf)

function F(x, p)
	@unpack α = p
	f = similar(x)

	f[1] = (-2x[1]+x[2]) + α * exp(x[1])
	f[2] = ( x[1]-2x[2]) + α * exp(x[2])

	return f
end

sol0 = zeros(2)
par = (α = 0.0, )
prob = BifurcationProblem(F, sol0, par, (@lens _.α); record_from_solution = (x,p) -> norminf(x))
```

```@example ANM
optanm = ContinuationPar(dsmin = 0.01, dsmax = 0.15, detect_bifurcation = 3, ds= 0.01, newton_options = NewtonPar(tol = 1e-9, verbose = false), n_inversion = 6, max_bisection_steps = 15, max_steps = 15, )

branm = continuation(prob, ANM(20, 1e-8), optanm, normC = norminf, verbosity = 2)
```

You can plot the result as usual:

```@example ANM
plot(branm)
```

You can also show the radius of convergence of each series:

```@example ANM
plot(branm, plotseries = true)
```

Finally, for each series, we ca evaluate the residual norm:

```@example ANM
plot()
for ii in eachindex(branm.polU)
	s = LinRange(-0*branm.radius[ii], branm.radius[ii], 20)
	plot!([branm.polp[ii].(s)], [norminf(F(branm.polU[ii](_s), BK.setparam(prob,branm.polp[ii](_s)))) for _s in s], legend = false, linewidth=5)#, marker=:d)
end
title!("")
```
### References

[^Charpentier]:> Charpentier, Isabelle, Bruno Cochelin, and Komlanvi Lampoh. “Diamanlab - An Interactive Taylor-Based Continuation Tool in MATLAB,” n.d., 12.

[^Rubbert]:> Rubbert, Lennart, Isabelle Charpentier, Simon Henein, and Pierre Renaud. “Higher-Order Continuation Method for the Rigid-Body Kinematic Design of Compliant Mechanisms”,  n.d., 18.
