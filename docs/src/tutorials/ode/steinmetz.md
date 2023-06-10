# Steinmetz-Larter model

```@contents
Pages = ["steinmetz.md"]
Depth = 3
```

> This example is found in the MatCont ecosystem.

The Steinmetz-Larter model is studied in the MatCont ecosystem because it is a simple example where a Chenciner bifurcation occurs.

$$\tag{E}\left\{\begin{array}{l}
\dot{A}=-k_1 A B X-k_3 A B Y+k_7-k_{-7} A, \\
\dot{B}=-k_1 A B X-k_3 A B Y+k_8, \\
\dot{X}=k_1 A B X-2 k_2 X^2+2 k_3 A B Y-k_4 X+k_6, \\
\dot{Y}=-k_3 A B Y+2 k_2 X^2-k_5 Y,
\end{array}\right.$$

This tutorial is also useful in that we show how to start periodic orbits continuation from solutions obtained from solving ODEs. Being a relatively advanced tutorial, we will not give too many details.

We start by coding the bifurcation problem.

```@example STEINMETZ
using Revise
using Test, ForwardDiff, Parameters, Plots, LinearAlgebra
using BifurcationKit, Test
const BK = BifurcationKit

norminf(x) = norm(x, Inf)

function SL!(du, u, p, t = 0)
	@unpack k1, k2, k3, k4, k5, k6, k7, k₋₇, k8 = p
	A,B,X,Y = u

	du[1]= -k1*A*B*X - k3*A*B*Y + k7 - k₋₇*A
	du[2] = -k1*A*B*X - k3*A*B*Y + k8
	du[3] = k1*A*B*X - 2*k2*X^2 + 2*k3*A*B*Y - k4*X+k6
	du[4] = -k3*A*B*Y + 2*k2*X^2 - k5*Y
	du
end
SL(u,p,t=0) = SL!(similar(u),u,p,t)

z0 = rand(4)
par_sl = (k1=0.1631021, k2=1250., k3=0.046875, k4=20., k5=1.104, k6=0.001, k₋₇=0.1175, k7=1.5, k8=0.75)
prob = BK.BifurcationProblem(SL, z0, par_sl, (@lens _.k8);)

function recordFromSolution(x, p) 
	xtt = BK.getPeriodicOrbit(p.prob, x, set(getParams(p.prob), BK.getLens(p.prob), p.p))
	return (max = maximum(xtt[1,:]),
			min = minimum(xtt[1,:]),
			period = getPeriod(p.prob, x, set(getParams(p.prob), BK.getLens(p.prob), p.p)))
end

function plotSolution(x, p; k...)
	xtt = BK.getPeriodicOrbit(p.prob, x, set(getParams(p.prob), BK.getLens(p.prob), p.p))
	plot!(xtt.t, xtt[:,:]'; label = "", k...)
end

argspo = (recordFromSolution = recordFromSolution,
	plotSolution = plotSolution)

nothing #hide
```

We obtain some trajectories as seeds for computing periodic orbits.

```@example STEINMETZ
using DifferentialEquations
alg_ode = Rodas5P()
prob_de = ODEProblem(SL!, z0, (0,136.), par_sl)
sol = solve(prob_de, alg_ode)
prob_de = ODEProblem(SL!, sol.u[end], (0,30.), sol.prob.p, reltol = 1e-11, abstol = 1e-13)
sol = solve(prob_de, alg_ode)
plot(sol)
```

We generate a shooting problem form the computed trajectories and continue the periodic orbits as function of $k_8$

```@example STEINMETZ
probsh, cish = generateCIProblem( ShootingProblem(M=4), prob, prob_de, sol, 16.; reltol = 1e-10, abstol = 1e-12, parallel = true)

solpo = newton(probsh, cish, NewtonPar(verbose = true))

_sol = BK.getPeriodicOrbit(probsh, solpo.u, sol.prob.p)
	plot(_sol.t, _sol[:,:]')

opts_br = ContinuationPar(pMin = 0., pMax = 20.0, ds = 0.002, dsmax = 0.05, nInversion = 8, detectBifurcation = 3, maxBisectionSteps = 25, nev = 4)
opts_po_cont = setproperties(opts_br, maxSteps = 60, saveEigenvectors = true, tolStability = 1e-3)
@set! opts_po_cont.newtonOptions.verbose = false
@set! opts_po_cont.newtonOptions.maxIter = 10
br_sh = continuation(deepcopy(probsh), cish, PALC(tangent = Bordered()), opts_po_cont;
	#verbosity = 3, plot = true,
	callbackN = BK.cbMaxNorm(10),
	argspo...)
```

## Curve of Fold points of periodic orbits

```@example STEINMETZ
opts_posh_fold = ContinuationPar(br_sh.contparams, detectBifurcation = 2, maxSteps = 35, pMax = 1.9, plotEveryStep = 10, dsmax = 4e-2, ds = 1e-2)
@set! opts_posh_fold.newtonOptions.tol = 1e-12
@set! opts_posh_fold.newtonOptions.verbose = true
fold_po_sh = @time continuation(br_sh, 2, (@lens _.k7), opts_posh_fold;
		#verbosity = 3, plot = true,
		detectCodim2Bifurcation = 2,
		startWithEigen = false,
		usehessian = false,
		jacobian_ma = :minaug,
		normC = norminf,
		callbackN = BK.cbMaxNorm(1e1),
		bdlinsolver = BorderingBLS(solver = DefaultLS(), checkPrecision = false),
		)
plot(fold_po_sh)
```

## Curve of NS points of periodic orbits
```@example STEINMETZ
opts_posh_ns = ContinuationPar(br_sh.contparams, detectBifurcation = 0, maxSteps = 35, pMax = 1.9, plotEveryStep = 10, dsmax = 4e-2, ds = 1e-2)
@set! opts_posh_ns.newtonOptions.tol = 1e-12
ns_po_sh = continuation(br_sh, 1, (@lens _.k7), opts_posh_ns;
		verbosity = 2, plot = false,
		detectCodim2Bifurcation = 2,
		startWithEigen = false,
		usehessian = false,
		jacobian_ma = :minaug,
		normC = norminf,
		callbackN = BK.cbMaxNorm(1e1),
		)
```

```@example STEINMETZ
plot(ns_po_sh, fold_po_sh, branchlabel = ["NS","Fold"])
```

