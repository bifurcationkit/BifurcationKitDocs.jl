# [🟠 Steinmetz-Larter model](@id steinmetz)

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

This tutorial is also useful in that we show how to start periodic orbits continuation from solutions obtained from solving ODEs. Being a relatively advanced tutorial, we do not give too many details.

We start by coding the bifurcation problem.

```@example STEINMETZ
using Revise
using ForwardDiff, Parameters, Plots
using BifurcationKit
const BK = BifurcationKit

function SL!(du, u, p, t = 0)
	@unpack k1, k2, k3, k4, k5, k6, k7, k₋₇, k8 = p
	A,B,X,Y = u

	du[1] = -k1*A*B*X - k3*A*B*Y + k7 - k₋₇*A
	du[2] = -k1*A*B*X - k3*A*B*Y + k8
	du[3] =  k1*A*B*X - 2*k2*X^2 + 2*k3*A*B*Y - k4*X+k6
	du[4] = -k3*A*B*Y + 2*k2*X^2 - k5*Y
	du
end

z0 = rand(4)
par_sl = (k1=0.1631021, k2=1250., k3=0.046875, k4=20., k5=1.104, k6=0.001, k₋₇=0.1175, k7=1.5, k8=0.75)
prob = BifurcationProblem(SL!, z0, par_sl, (@lens _.k8);)

# record variables for plotting
function recordFromSolution(x, p) 
	xtt = BK.get_periodic_orbit(p.prob, x, p.p)
	return (max = maximum(xtt[1,:]),
			min = minimum(xtt[1,:]),
			period = getperiod(p.prob, x, p.p))
end

# plotting function
function plotSolution(x, p; k...)
	xtt = BK.get_periodic_orbit(p.prob, x, p.p)
	plot!(xtt.t, xtt[:,:]'; label = "", k...)
end

# group parameters
argspo = (record_from_solution = recordFromSolution,
	plot_solution = plotSolution)

nothing #hide
```

We obtain some trajectories as seeds for computing periodic orbits.

```@example STEINMETZ
using DifferentialEquations
alg_ode = Rodas5()
prob_de = ODEProblem(SL!, z0, (0,136.), par_sl)
sol = solve(prob_de, alg_ode)
prob_de = ODEProblem(SL!, sol.u[end], (0,30.), sol.prob.p, reltol = 1e-11, abstol = 1e-13)
sol = solve(prob_de, alg_ode)
plot(sol)
```

## Computation with Shooting
We generate a shooting problem from the computed trajectories and continue the periodic orbits as function of $k_8$

```@example STEINMETZ
probsh, cish = generate_ci_problem( ShootingProblem(M=4), prob, prob_de, sol, 16.; reltol = 1e-10, abstol = 1e-12, parallel = true)

opts_po_cont = ContinuationPar(p_min = 0., p_max = 20.0, ds = 0.002, dsmax = 0.05, n_inversion = 8, detect_bifurcation = 3, max_bisection_steps = 25, nev = 4, max_steps = 60, save_eigenvectors = true, tol_stability = 1e-3)
@set! opts_po_cont.newton_options.verbose = false
@set! opts_po_cont.newton_options.max_iterations = 10
br_sh = continuation(deepcopy(probsh), cish, PALC(tangent = Bordered()), opts_po_cont;
	# verbosity = 3, plot = true,
	callback_newton = BK.cbMaxNorm(10),
	argspo...)
scene = plot(br_sh)
```

### Curve of Fold points of periodic orbits

```@example STEINMETZ
opts_posh_fold = ContinuationPar(br_sh.contparams, detect_bifurcation = 2, max_steps = 35, p_max = 1.9, plot_every_step = 10, dsmax = 4e-2, ds = 1e-2)
@set! opts_posh_fold.newton_options.tol = 1e-12
# @set! opts_posh_fold.newton_options.verbose = true
fold_po_sh = @time continuation(deepcopy(br_sh), 2, (@lens _.k7), opts_posh_fold;
		# verbosity = 2, plot = true,
		detect_codim2_bifurcation = 0,
		update_minaug_every_step = 1,
		start_with_eigen = true,
		usehessian = false,
		jacobian_ma = :minaug,
		normC = norminf,
		callback_newton = BK.cbMaxNorm(1e1),
		# bdlinsolver = BorderingBLS(solver = DefaultLS(), check_precision = false),
		)
plot(fold_po_sh)
```

### Curve of NS points of periodic orbits
```@example STEINMETZ
opts_posh_ns = ContinuationPar(br_sh.contparams, detect_bifurcation = 0, max_steps = 35, p_max = 1.9, plot_every_step = 10, dsmax = 4e-2, ds = 1e-2)
@set! opts_posh_ns.newton_options.tol = 1e-12
# @set! opts_posh_ns.newton_options.verbose = true
ns_po_sh = continuation(deepcopy(br_sh), 1, (@lens _.k7), opts_posh_ns;
		# verbosity = 2, plot = true,
		detect_codim2_bifurcation = 2,
		update_minaug_every_step = 1,
		start_with_eigen = false,
		jacobian_ma = :minaug,
		normC = norminf,
		callback_newton = BK.cbMaxNorm(1e1),
		)
```

```@example STEINMETZ
plot(ns_po_sh, fold_po_sh, branchlabel = ["NS","Fold"])
```

## Computation with collocation

```@example STEINMETZ
probcoll, cicoll = generate_ci_problem( PeriodicOrbitOCollProblem(50, 4), prob, sol, 16.)

opts_po_cont = ContinuationPar(p_min = 0., p_max = 2.0, 
	ds = 0.002, dsmax = 0.05, 
	# n_inversion = 6,
	nev = 4,
	max_steps = 50, 
	tol_stability = 1e-5)
br_coll = continuation(probcoll, cicoll, PALC(tangent = Bordered()), opts_po_cont;
    # verbosity = 3, plot = true,
    callback_newton = BK.cbMaxNorm(10),
    argspo...)
```

### Curve of Fold points of periodic orbits

```@example STEINMETZ
opts_pocl_fold = ContinuationPar(br_coll.contparams, detect_bifurcation = 1, plot_every_step = 10, dsmax = 4e-2)
fold_po_cl = @time continuation(br_coll, 2, (@lens _.k7), opts_pocl_fold;
        # verbosity = 3, plot = true,
        detect_codim2_bifurcation = 2,
        update_minaug_every_step = 1,
        start_with_eigen = false,
        usehessian = true,
        jacobian_ma = :minaug,
        normC = norminf,
        callback_newton = BK.cbMaxNorm(1e1),
        )
```

### Curve of NS points of periodic orbits

```@example STEINMETZ
opts_pocl_ns = ContinuationPar(br_coll.contparams, detect_bifurcation = 1, plot_every_step = 10, dsmax = 4e-2)
ns_po_cl = continuation(br_coll, 1, (@lens _.k7), opts_pocl_ns;
        # verbosity = 3, plot = true,
        detect_codim2_bifurcation = 2,
        update_minaug_every_step = 1,
        start_with_eigen = false,
        jacobian_ma = :minaug,
        normC = norminf,
        callback_newton = BK.cbMaxNorm(1e1),
        )
```

```@example STEINMETZ
plot(ns_po_cl, fold_po_cl, branchlabel = ["NS","Fold"])
```