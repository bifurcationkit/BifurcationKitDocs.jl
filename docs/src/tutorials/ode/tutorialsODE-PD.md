# [🟡 Period doubling in Lur'e problem (PD aBS)](@id pdlure)

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

The model is interesting because there is a period doubling bifurcation and we want to show the branch switching capabilities of `BifurcationKit.jl` in this case. We provide 3 different ways to compute the periodic orbits and highlight their pro / cons.

It is easy to encode the ODE as follows

```@example TUTLURE
using Revise, Plots
using BifurcationKit
const BK = BifurcationKit

recordFromSolution(x, p; k...) = (u1 = x[1], u2 = x[2])

function lur!(dz, u, p, t = 0)
	(;α, β) = p
	x, y, z = u
	dz[1] = y
	dz[2] =	z
	dz[3] = -α * z - β * y - x + x^2
	dz
end

# bifurcation problem
prob_bif = ODEBifProblem(lur!, zeros(3), (α = -1.0, β = 1.0), (@optic _.α);
    record_from_solution = recordFromSolution)
nothing #hide
```

We first compute the branch of equilibria

```@example TUTLURE
# continuation options
opts_br = ContinuationPar(p_min = -1.4, p_max = 1.8, nev = 3)

# computation of the branch
br = continuation(prob_bif, PALC(), opts_br)

scene = plot(br)
```

With detailed information:

```@example TUTLURE
br
```

We note the Hopf bifurcation point which we shall investigate now.

## Periodic orbits with orthogonal collocation

We compute the branch of periodic orbits from the Hopf bifurcation point.
We rely on a the state of the art method for computing periodic orbits of ODE: orthogonal collocation.

We first define a plotting function and a record function which are used for all cases below:

```@example TUTLURE
# plotting function
function plotPO(x, p; k...)
	xtt = get_periodic_orbit(p.prob, x, p.p)
	plot!(xtt.t, xtt[1,:]; markersize = 2, k...)
	plot!(xtt.t, xtt[2,:]; k...)
	plot!(xtt.t, xtt[3,:]; legend = false, k...)
end

# record function
function recordPO(x, p; k...)
	xtt = get_periodic_orbit(p.prob, x, p.p)
	period = getperiod(p.prob, x, p.p)
	mn, mx = extrema(xtt[1,:])
	return (;max = mx, min = mn, period)
end
```

Continuation of periodic orbits from the Hopf point:

```@example TUTLURE
# continuation parameters
opts_po_cont = ContinuationPar(opts_br, dsmax = 0.03, dsmin = 1e-4, max_steps = 80, tol_stability = 1e-4, plot_every_step = 20)

br_po = continuation(
	br, 1, opts_po_cont,
	PeriodicOrbitOCollProblem(40, 4);
	plot = true,
	record_from_solution = recordPO,
	plot_solution = (x, p; k...) -> begin
		plotPO(x, p; k...)
		# plot previous branch
		plot!(br, subplot = 1, putbifptlegend = false)
		end,
	normC = norminf)

scene = plot(br, br_po)
```

Note that you can compute the PD normal form

```@example TUTLURE
get_normal_form(br_po, 1)
```

We provide Automatic Branch Switching from the PD point and computing the bifurcated branch is as simple as:

```@example TUTLURE
# aBS from PD
br_po_pd = continuation(deepcopy(br_po), 1, ContinuationPar(br_po.contparams, max_steps = 100, dsmax = 0.02, plot_every_step = 10, ds = 0.005);
	plot_solution = (x, p; k...) -> begin
		plotPO(x, p; k...)
		## add previous branch
		plot!(br_po; subplot = 1)
	end,
	record_from_solution = recordPO,
	normC = norminf,
	)

scene = plot(br_po, br_po_pd, title = "Collocation based")
```

## Periodic orbits with Parallel Standard Shooting

We use a different method to compute periodic orbits: we rely on a fixed point of the flow. To compute the flow, we use `OrdinaryDiffEq.jl`. This way of computing periodic orbits should be less precise than the previous one. We rely on parallel multiple shooting. Finally, please note the close similarity to the code of the previous part. As before, we first rely on Hopf **aBS**.

```@example TUTLURE
import OrdinaryDiffEq as ODE

# ODE problem for using DifferentialEquations
prob_ode = ODE.ODEProblem(lur!, prob_bif.u0, (0, 1), prob_bif.params; abstol = 1e-12, reltol = 1e-10)

# continuation parameters
opts_po_cont = ContinuationPar(opts_br, dsmax = 0.03, newton_options = NewtonPar(tol = 1e-10), tol_stability = 1e-5, n_inversion = 8, max_steps = 100)

br_po = continuation(
	br, 1, opts_po_cont,
	# parallel shooting functional with 10 sections
	ShootingProblem(15, prob_ode, ODE.Vern9(); parallel = true);
	# plot = true,
	record_from_solution = recordPO,
	plot_solution = plotPO,
	# limit the residual, useful to help DifferentialEquations
	callback_newton = BK.cbMaxNorm(10),
	normC = norminf)

plot(br_po)
```

Note that you can compute the PD normal form

```@example TUTLURE
get_normal_form(br_po, 1; autodiff = true)
```

We provide Automatic Branch Switching from the PD point and computing the bifurcated branch is as simple as:

```@example TUTLURE
# aBS from PD
br_po_pd = continuation(deepcopy(br_po), 1, 
	ContinuationPar(br_po.contparams, max_steps = 10, ds = -0.005);
	plot = true, verbosity = 2,
	plot_solution = (x, p; k...) -> begin
		plotPO(x, p; k...)
		# add previous branch
		plot!(br_po; subplot = 1)
	end,
	record_from_solution = recordPO,
	autodiff_nf = true,
	normC = norminf,
	callback_newton = BK.cbMaxNorm(10),
	)

scene = plot(br, br_po, br_po_pd, title = "Shooting based")
```

## Branch of periodic orbits with finite differences

We use finite differences to discretize the problem for finding periodic orbits. We appeal to automatic branch switching from the Hopf point as follows

```@example TUTLURE
# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.02, dsmin = 1e-4, p_max = 1.1, max_steps = 80, tol_stability = 1e-4)

br_po = continuation(
	br, 1, opts_po_cont,
	PeriodicOrbitTrapProblem(M = 120);
	record_from_solution = recordPO,
	plot_solution = (x, p; k...) -> begin
		plotPO(x, p; k...)
		## plot previous branch
		plot!(br, subplot=1, putbifptlegend = false)
		end,
	normC = norminf)

scene = plot(br, br_po)
```

Two period doubling bifurcations were detected. We shall now compute the branch of periodic orbits from these PD points. We do not provide Automatic Branch Switching for Trapezoid method as we do not have yet the PD normal form computed for `PeriodicOrbitTrapProblem`. Hence, it takes some trial and error to find the `ampfactor` of the PD branch.

> This is like in MatCont and Auto-07p here...

```@example TUTLURE
# aBS from PD
br_po_pd = continuation(deepcopy(br_po), 1, ContinuationPar(br_po.contparams, max_steps = 70);
	# plot = true,
	ampfactor = .2, δp = -0.005,
	plot_solution = (x, p; k...) -> begin
		plotPO(x, p; k...)
		# add previous branch
		plot!(br_po; legend=false, subplot=1)
	end,
	record_from_solution = recordPO,
	normC = norminf
	)
Scene = title!("")
```

```@example TUTLURE
plot(br, br_po, br_po_pd, title = "Trapezoid based")
```