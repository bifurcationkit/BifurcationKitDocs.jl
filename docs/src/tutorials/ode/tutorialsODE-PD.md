# Period doubling in Lur'e problem (PD aBS)

```@contents
Pages = ["tutorialsODE-PD.md"]
Depth = 3
```

The following model is an adaptive control system of Lur’e type. It is an example from the MatCont library.

$$\left\{\begin{array}{l}
\dot{x}=y \\
\dot{y}=z \\
\dot{z}=-\alpha z-\beta y-x+x^{2}
\end{array}\right.$$


The model is interesting because there is a period doubling bifurcation and we want to show the branch switching capabilities of `BifurcationKit.jl` in this case. We provide 3 different ways to compute this periodic orbits and highlight their pro / cons.

It is easy to encode the ODE as follows

```@example TUTLURE
using Revise, Parameters, Plots, LinearAlgebra
using BifurcationKit
const BK = BifurcationKit

norminf(x) = norm(x, Inf)
recordFromSolution(x, p) = (u1 = x[1], u2 = x[2])

function lur!(dz, z, p, t)
	@unpack α, β = p
	x, y, z = z
	dz[1] = y
	dz[2] =	z
	dz[3] = -α * z - β * y - x + x^2
	dz
end

lur(z, p) = lur!(similar(z), z, p, 0)

# parameters
par_lur = (α = -1.0, β = 1.)

# initial guess
z0 = zeros(3)

# bifurcation problem
prob = BifurcationProblem(lur, z0, par_lur, (@lens _.α);
    recordFromSolution = recordFromSolution)
nothing #hide
```

We first compute the branch of equilibria

```@example TUTLURE
# continuation options
opts_br = ContinuationPar(pMin = -1.4, pMax = 1.8, ds = 0.01, dsmax = 0.01, plotEveryStep = 20, maxSteps = 1000)

# computation of the branch
br = continuation(prob, PALC(), opts_br)

scene = plot(br)
```

With detailed information:

```@example TUTLURE
br
```

We note the Hopf bifurcation point which we shall investigate now.

## Branch of periodic orbits with finite differences

We compute the branch of periodic orbits from the Hopf bifurcation point.
We first define a plotting function and a record function which are used for all cases below:

```@example TUTLURE
# plotting function
function plotPO(x, p; k...)
	xtt = BK.getPeriodicOrbit(p.prob, x, p.p)
	plot!(xtt.t, xtt[1,:]; markersize = 2, k...)
	plot!(xtt.t, xtt[2,:]; k...)
	plot!(xtt.t, xtt[3,:]; legend = false, k...)
end

# record function
function recordPO(x, p)
	xtt = BK.getPeriodicOrbit(p.prob, x, p.p)
	period = BK.getPeriod(p.prob, x, p.p)
	return (max = maximum(xtt[1,:]), min = minimum(xtt[1,:]), period = period)
end
```

We use finite differences to discretize the problem for finding periodic orbits. We appeal to automatic branch switching as follows

```@example TUTLURE
# newton parameters
optn_po = NewtonPar(tol = 1e-8,  maxIter = 25)

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.02, ds= 0.01, dsmin = 1e-4, pMax = 1.8, pMin=-1., maxSteps = 80, newtonOptions =  optn_po, tolStability = 1e-4)

Mt = 120 # number of time sections
br_po = continuation(
	br, 1, opts_po_cont,
	PeriodicOrbitTrapProblem(M = Mt);
	δp = 0.01,
	recordFromSolution = recordPO,
	plotSolution = (x, p; k...) -> begin
		plotPO(x, p; k...)
		## plot previous branch
		plot!(br, subplot=1, putbifptlegend = false)
		end,
	normC = norminf)

scene = plot(br, br_po)
```

Two period doubling bifurcations were detected. We shall now compute the branch of periodic orbits from these PD points. We do not provide Automatic Branch Switching as we do not have the PD normal form computed in `BifurcationKit`. Hence, it takes some trial and error to find the `ampfactor` of the PD branch.

```@example TUTLURE
# aBS from PD
br_po_pd = continuation(br_po, 1, setproperties(br_po.contparams, maxSteps = 40);
	verbosity = 3, plot = true,
	ampfactor = .2, δp = -0.005,
	usedeflation = false,
	plotSolution = (x, p; k...) -> begin
		plotPO(x, p; k...)
		## add previous branch
		plot!(br_po; legend=false, subplot=1)
	end,
	recordFromSolution = recordPO,
	normC = norminf
	)
Scene = title!("")
```

```@example TUTLURE
plot(br, br_po, br_po_pd)
```

## Periodic orbits with Parallel Standard Shooting

We use a different method to compute periodic orbits: we rely on a fixed point of the flow. To compute the flow, we use `DifferentialEquations.jl`. This way of computing periodic orbits should be more precise than the previous one. We use a particular instance called multiple shooting which is computed in parallel. This is an additional advantage compared to the previous method. Finally, please note the close similarity to the code of the previous part. As before, we first rely on Hopf **aBS**.

```@example TUTLURE
using DifferentialEquations

# ODE problem for using DifferentialEquations
probsh = ODEProblem(lur!, copy(z0), (0., 1000.), par_lur; abstol = 1e-12, reltol = 1e-10)

# newton parameters
optn_po = NewtonPar(verbose = true, tol = 1e-12, maxIter = 25)

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.02, ds= -0.001, dsmin = 1e-4, maxSteps = 130, newtonOptions = optn_po, tolStability = 1e-5, detectBifurcation = 3, plotEveryStep = 10, nInversion = 6, nev = 2)

br_po = continuation(
	br, 1, opts_po_cont,
	# parallel shooting functional with 15 sections
	ShootingProblem(15, probsh, Rodas5(); parallel = true);
	# first parameter value on the branch
	δp = 0.005,
	verbosity = 3, plot = true,
	recordFromSolution = recordPO,
	plotSolution = plotPO,
	# limit the residual, useful to help DifferentialEquations
	callbackN = BK.cbMaxNorm(10),
	normC = norminf)

scene = title!("")
```

We do not provide Automatic Branch Switching as we do not have the PD normal form computed in `BifurcationKit`. Hence, it takes some trial and error to find the `ampfactor` of the PD branch.

```@example TUTLURE
# aBS from PD
br_po_pd = continuation(br_po, 1, setproperties(br_po.contparams, maxSteps = 40, dsmax = 0.03, plotEveryStep = 10, ds = 0.01);
	verbosity = 3, plot = true,
	ampfactor = .1, δp = -0.005,
	plotSolution = (x, p; k...) -> begin
		plotPO(x, p; k...)
		## add previous branch
		plot!(br_po; subplot = 1)
	end,
	recordFromSolution = recordPO,
	normC = norminf,
	callbackN = BK.cbMaxNorm(10),
	)

scene = plot(br_po, br_po_pd)
```

## Periodic orbits with orthogonal collocation

We now rely on a the state of the art method for computing periodic orbits of ODE: orthogonal collocation.

```@example TUTLURE
# newton parameters
optn_po = NewtonPar(tol = 1e-10,  maxIter = 25)

# continuation parameters
opts_po_cont = ContinuationPar(opts_br, dsmax = 0.03, ds= 0.01, dsmin = 1e-4, maxSteps = 80, newtonOptions = optn_po, tolStability = 1e-4, plotEveryStep = 20, nInversion = 6)

br_po = continuation(
	br, 1, opts_po_cont,
	PeriodicOrbitOCollProblem(20, 4);
	ampfactor = 1., δp = 0.01,
	verbosity = 2,	plot = true,
	recordFromSolution = recordPO,
	plotSolution = (x, p; k...) -> begin
		plotPO(x, p; k...)
		## plot previous branch
		plot!(br, subplot=1, putbifptlegend = false)
		end,
	normC = norminf)

scene = plot(br, br_po)
```
We do not provide Automatic Branch Switching as we do not have the PD normal form computed in `BifurcationKit`. Hence, it takes some trial and error to find the `ampfactor` of the PD branch.

```@example TUTLURE
# aBS from PD
br_po_pd = continuation(br_po, 1, setproperties(br_po.contparams, maxSteps = 50, ds = 0.01, plotEveryStep = 10);
	verbosity = 3, plot = true,
	ampfactor = .3, δp = -0.005,
	plotSolution = (x, p; k...) -> begin
		plotPO(x, p; k...)
		## add previous branch
		plot!(br_po; subplot = 1)
	end,
	recordFromSolution = recordPO,
	normC = norminf,
	callbackN = BK.cbMaxNorm(10),
	)

scene = plot(br_po, br_po_pd)
```
