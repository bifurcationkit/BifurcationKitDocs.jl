# Periodic predator-prey model

```@contents
Pages = ["tutorialsCodim2PO.md"]
Depth = 3
```

The following is a periodically forced predator-prey model studied in [^Kuznetsov] using shooting technics.

$$\tag{E}\left\{\begin{aligned}
\dot{x} & =r\left(1-\frac{x}{K}\right) x-\frac{a x y}{b_0(1+\varepsilon u)+x} \\
\dot{y} & =e \frac{a x y}{b_0(1+\varepsilon u)+x}-d y \\
\dot{u} & =u-2 \pi v-\left(u^2+v^2\right) u \\
\dot{v} & =2 \pi u+v-\left(u^2+v^2\right) v
\end{aligned}\right.$$

This tutorial is useful in that we show how to start periodic orbits continuation from solutions obtained from solving ODEs. We are interested in cases where $(u,v)\neq (0,0)$ for which we have only periodic solutions. Thus, we cannot rely on branching from equilibria to find these periodic orbits.

```@example TUTPPREY
using Revise
using Test, ForwardDiff, Parameters, Setfield, Plots, LinearAlgebra
using BifurcationKit, Test
const BK = BifurcationKit

norminf(x) = norm(x, Inf)

function Pop!(du, X, p, t = 0)
	@unpack r,K,a,ϵ,b0,e,d, = p
	x, y, u, v = X
	p = a * x / (b0 * (1 + ϵ * u) + x)
	du[1] = r * (1 - x/K) * x - p * y
	du[2] = e * p * y - d * y
	s = u^2 + v^2
	du[3] = u-2pi * v - s * u
	du[4] = 2pi * u + v - s * v
	du
end
Pop(u,p) = Pop!(similar(u),u,p,0)

par_pop = ( K = 1., r = 6.28, a = 12.56, b0 = 0.25, e = 1., d = 6.28, ϵ = 0.2, )

z0 = [0.1,0.1,1,0]

prob = BK.BifurcationProblem(Pop, z0, par_pop, (@lens _.b0);
	recordFromSolution = (x, p) -> (x = x[1], y = x[2], u = x[3]))

opts_br = ContinuationPar(pMin = 0., pMax = 20.0, ds = 0.002, dsmax = 0.01, nInversion = 6, detectBifurcation = 3, maxBisectionSteps = 25, nev = 4, maxSteps = 20000)

nothing #hide
```

## Simulation of the model

This is very straightforward thanks to `SciML`.

```@example TUTPPREY
using DifferentialEquations
prob_de = ODEProblem(Pop!, z0, (0,200.), par_pop)
sol = solve(prob_de, Rodas5())
prob_de = ODEProblem(Pop!, sol.u[end], (0,3.), par_pop, reltol = 1e-8, abstol = 1e-10)
sol = solve(prob_de, Rodas5())
plot(sol)
```

## Utility functions

We start with two helper functions that record and plot the periodic orbits. The following works for shooting, collocation and trapezoid methods for computing periodic orbits.

```@example TUTPPREY
argspo = (recordFromSolution = (x, p) -> begin
		xtt = BK.getPeriodicOrbit(p.prob, x, set(getParams(p.prob), BK.getLens(p.prob), p.p))
		return (max = maximum(xtt[1,:]),
				min = minimum(xtt[1,:]),
				period = getPeriod(p.prob, x, set(getParams(p.prob), BK.getLens(p.prob), p.p)))
	end,
	plotSolution = (x, p; k...) -> begin
		xtt = BK.getPeriodicOrbit(p.prob, x, set(getParams(p.prob), BK.getLens(p.prob), p.p))
		plot!(xtt.t, xtt[1,:]; label = "x", k...)
		plot!(xtt.t, xtt[2,:]; label = "y", k...)
		# plot!(br; subplot = 1, putspecialptlegend = false)
	end)
```


## Periodic orbits with trapezoid method

We compute periodic orbits of (E) using the Trapezoid method.
We are now equipped to build a periodic orbit problem from a solution `sol::ODEProblem`.

```@example TUTPPREY
# this is the function which builds probtrap from sol
probtrap, ci = BK.generateCIProblem(PeriodicOrbitTrapProblem(M = 150),
	prob, sol, 2.)

opts_po_cont = setproperties(opts_br, maxSteps = 50, tolStability = 1e-8)
brpo_fold = continuation(probtrap, ci, PALC(), opts_po_cont;
	verbosity = 3, plot = true,
	argspo...
	)

scene = plot(brpo_fold)
```

We continue w.r.t. to $\epsilon$ and find a period-doubling bifurcation.

```@example TUTPPREY
prob2 = @set probtrap.prob_vf.lens = @lens _.ϵ
brpo_pd = continuation(prob2, ci, PALC(), opts_po_cont;
	verbosity = 3, plot = true,
	argspo...
	)
scene = plot(brpo_pd)
```

## Periodic orbits with parallel standard shooting

We are now ready to build a periodic orbit problem from a solution `sol::ODEProblem`.

```@example TUTPPREY
probsh, cish = generateCIProblem( ShootingProblem(M=3),
	prob, prob_de, sol, 2.; alg = Rodas5())

opts_po_cont = setproperties(opts_br, maxSteps = 50, tolStability = 1e-3)
br_fold_sh = continuation(probsh, cish, PALC(tangent = Bordered()), opts_po_cont;
	verbosity = 3, plot = true,
	argspo...
)

scene = plot(br_fold_sh)
```

We continue w.r.t. to $\epsilon$ and find a period-doubling bifurcation.

```@example TUTPPREY
probsh2 = @set probsh.lens = @lens _.ϵ
brpo_pd_sh = continuation(probsh2, cish, PALC(), opts_po_cont;
	verbosity = 3, plot = true,
	argspo...
	)
scene = plot(brpo_pd_sh)
```

## Periodic orbits with orthogonal collocation

We do the same as in the previous section but using orthogonal collocation. This is the most reliable and precise method for ODE. When the dimension of the ODE is large, it becomes prohibitive.

```@example TUTPPREY
# this is the function which builds probcoll from sol
probcoll, ci = generateCIProblem(PeriodicOrbitOCollProblem(26, 3; updateSectionEveryStep = 0),
	prob, sol, 2.)

opts_po_cont = setproperties(opts_br, maxSteps = 50, tolStability = 1e-8)
brpo_fold = continuation(probcoll, ci, PALC(), opts_po_cont;
	verbosity = 3, plot = true,
	argspo...
	)
scene = plot(brpo_fold)
```

We continue w.r.t. to $\epsilon$ and find a period-doubling bifurcation.

```@example TUTPPREY
prob2 = @set probcoll.prob_vf.lens = @lens _.ϵ
brpo_pd = continuation(prob2, ci, PALC(), ContinuationPar(opts_po_cont, dsmax = 5e-3);
	verbosity = 3, plot = true,
	argspo...
	)

scene = plot(brpo_pd)
```

## References

[^Kuznetsov]:> Yu.A. Kuznetsov, S. Muratori, and S. Rinaldi. Bifurcations and chaos in a periodic predator-prey model, Internat. J. Bifur. Chaos Appl. Sci. Engr., 2 (1992), 117-128.
