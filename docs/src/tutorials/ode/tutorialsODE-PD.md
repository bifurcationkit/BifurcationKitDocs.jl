# [🟡 Branch points and period doubling of periodic orbits (aBS)](@id pdlure)

```@contents
Pages = ["tutorialsODE-PD.md"]
Depth = 3
```

We consider the following three-dimensional vector field

$$\left\{\begin{aligned}
r\dot{x} &= {-(\beta+\nu)x + \beta y - a_3 x^3 + b_3(y-x)^3},\\
\dot{y} &= \beta x - (\beta+\gamma)y - z - b_3(y-x)^3,\\
\dot{z} &= y,
\end{aligned}\right.$$

with parameters $\nu,\beta,\gamma,r,a_3,b_3$. We use $\nu$ as the continuation parameter. This model presents a Hopf bifurcation from the trivial equilibrium as well as branch points (fold / pitchfork) of periodic orbits. It is a nice test case to show **automatic branch switching (aBS)** for periodic orbits, both from a branch point and from a period-doubling point.

It is easy to encode the vector field as follows

```@example TUTFREIRE
using Revise, CairoMakie
import BifurcationKit as BK
BK.set_plot_backend!(BK.BK_Makie()) # hide

function tor!(du, u, p, t = 0)
	(;nu, beta, gamma, r, a3, b3) = p
	x, y, z = u
	du[1] = ( -(beta + nu) * x + beta * y - a3 * x^3 + b3 * (y - x)^3 ) / r
	du[2] =   beta * x - (beta + gamma) * y - z - b3 * (y - x)^3
	du[3] = y
	return du
end

par_tor = (nu = -0.9, beta = 0.5, gamma = -0.6, r = 0.6, a3 = 0.328578, b3 = 0.933578)

prob = BK.ODEBifProblem(tor!, zeros(3), par_tor, (BK.@optic _.nu);
	record_from_solution = (x, p; k...) -> (x = x[1], y = x[2], z = x[3]))

nothing #hide
```

## Branch of equilibria

We first compute the branch of equilibria w.r.t. $\nu$.

```@example TUTFREIRE
opts_br = BK.ContinuationPar(p_min = -1.0, p_max = 0.0, n_inversion = 8, ds = 0.001)
br = BK.continuation(prob, BK.PALC(), opts_br; normC = BK.norminf,)
f, ax = BK.plot(br, dash_unstable_style = true)
f
```

We detect a Hopf bifurcation (index 1) from which we shall compute a branch of periodic orbits, and a branch point of equilibria. We next define the recording and plotting functions used below.

```@example TUTFREIRE
function recordPO(x, p; iter, state, k...)
	probpo = BK.getprob(iter)
	par = BK.getparams(iter, state)
	xtt = BK.get_periodic_orbit(probpo, x, par)
	return (max = maximum(xtt[1, :]),
		min = minimum(xtt[1, :]),
		period = BK.getperiod(probpo, x, par))
end

function plotPO(ax, x, p; iter, state, k...)
	probpo = BK.getprob(iter)
	par = BK.getparams(iter, state)
	xtt = BK.get_periodic_orbit(probpo, x, par)
	lines!(ax, xtt.t, xtt[1, :]; label = "x", k...)
	lines!(ax, xtt.t, xtt[2, :]; label = "y", k...)
	lines!(ax, xtt.t, xtt[3, :]; label = "z", k...)
end

argspo = (record_from_solution = recordPO,  plot_solution = plotPO)

nothing #hide
```

## Periodic orbits with orthogonal collocation

We compute the branch of periodic orbits from the Hopf bifurcation point using orthogonal collocation.

We now compute the branch of periodic orbits:

```@example TUTFREIRE
opts_po_cont = BK.ContinuationPar(opts_br, ds = 0.001, dsmax = 0.015,
	tol_stability = 1e-9, p_min = -0.7)

prob_po = BK.Collocation(30, 4; meshadapt = false, jacobian = BK.DenseAnalytical())

br_po = BK.continuation(br, 1, opts_po_cont, prob_po;
	alg = BK.PALC(tangent = BK.Bordered()),
	argspo...,
	normC = BK.norminf,
	linear_algo = BK.COPBLS())

f, ax = BK.plot(br_po, dash_unstable_style = true)
f
```

The branch of periodic orbits presents several branch points. The first one is a fold of periodic orbits, the second one is a pitchfork.

### Normal form of the branch points

!!! tip "Skip the normal form computation"
    You can skip this in a first read. This is not needed in what follows.

We compute the normal form of the two first branch points using the Iooss method (`prm = Val(false)`).

```@example TUTFREIRE
bp_fold = BK.get_normal_form(br_po, 1; prm = Val(false))
```

The second branch point is a pitchfork:

```@example TUTFREIRE
bp_pf = BK.get_normal_form(br_po, 2; prm = Val(false))
```

### Automatic branch switching from a branch point

Since the second branch point is a pitchfork, automatic branch switching is available. We simply continue from the second branch point:

```@example TUTFREIRE
br_po2 = BK.continuation(deepcopy(br_po), 2;
	δp = 0.001,
	record_from_solution = recordPO,
	normC = BK.norminf)

f, ax = BK.plot(br_po, br_po2, branchlabel = ["BP branch", "aBS from BP"], dash_unstable_style = true)
f
```

The bifurcated branch contains a Neimark-Sacker point and a period-doubling point.

### Automatic branch switching from a period-doubling point

We now branch from the period-doubling point detected on the previous branch. 

!!! tip "Skip the normal form computation"
    You can skip this in a first read. This is not needed in what follows.

We can compute its normal form just for fun...

```@example TUTFREIRE
ipd = findfirst(sp -> sp.type == :pd, br_po2.specialpoint)
bp_pd = BK.get_normal_form(br_po2, ipd)
```

Let's branch:

```@example TUTFREIRE
br_po_pd = BK.continuation(deepcopy(br_po2), ipd;
	δp = 0.0005,
	record_from_solution = recordPO,
	normC = BK.norminf)

f, ax = BK.plot(br_po2, br_po_pd, branchlabel = ["BP branch", "aBS from PD"], dash_unstable_style = true)
f
```

## Periodic orbits with standard shooting

We now compute the branch of periodic orbits with standard shooting. For shooting, the normal form of branch points and period-doubling points is computed from the Poincaré return map.

```@example TUTFREIRE
import OrdinaryDiffEq as ODE

prob_ode = ODE.ODEProblem(tor!, zeros(3), (0, 1), par_tor; abstol = 1e-12, reltol = 1e-10)
prob_sh = BK.Shooting(3, prob_ode, ODE.Vern9(); parallel = true)

opts_po_sh = BK.ContinuationPar(opts_br, ds = -0.001, dsmax = 0.01, tol_stability = 1e-4, p_min = -0.7)

br_po_sh = BK.continuation(br, 1, opts_po_sh, prob_sh;
	record_from_solution = recordPO,
	normC = BK.norminf,
	δp = 0.001)

f, ax = BK.plot(br_po_sh, dash_unstable_style = true)
f
```

!!! tip "Skip the normal form computation"
    You can skip this in a first read. This is not needed in what follows.

For information, the normal form of the second branch point (a pitchfork) is computed with the Poincaré return map method (`prm = Val(true)`, the default for shooting):

```@example TUTFREIRE
bp_pf_sh = BK.get_normal_form(br_po_sh, 2)
```

### Automatic branch switching from a branch point
We branch from this branch point:

```@example TUTFREIRE
br_po_sh_bp = BK.continuation(deepcopy(br_po_sh), 2;
	δp = 0.0005,
	record_from_solution = recordPO,
	normC = BK.norminf)

f, ax = BK.plot(br_po_sh, br_po_sh_bp, branchlabel = ["BP branch", "aBS from BP"], dash_unstable_style = true)
f
```

### Automatic branch switching from a period-doubling point
The bifurcated branch has a Neimark-Sacker point and a period-doubling point. We branch from the period-doubling point:

```@example TUTFREIRE
ipd_sh = findfirst(sp -> sp.type == :pd, br_po_sh_bp.specialpoint)

br_po_sh_pd = BK.continuation(deepcopy(br_po_sh_bp), ipd_sh;
	δp = 0.0005, ampfactor = 0.5,
	record_from_solution = recordPO,
	callback_newton = BK.cbMaxNorm(1.0),
	normC = BK.norminf)

f, ax = BK.plot(br_po_sh_bp, br_po_sh_pd, branchlabel = ["BP branch", "aBS from PD"], dash_unstable_style = true)
f
```

## Periodic orbits with the trapezoid method

Finally, we use the trapezoid method. For this discretisation, the normal form of the branch point is not fully implemented: only the eigenvectors are computed and branch switching requires a manual choice of `δp` and `ampfactor`.

```@example TUTFREIRE
prob_trap = BK.Trapeze(M = 100; jacobian = BK.Dense())
opts_po_trap = BK.ContinuationPar(opts_br, ds = -0.001, dsmax = 0.01, tol_stability = 1e-4, p_min = -0.7)

br_po_trap = BK.continuation(br, 1, opts_po_trap, prob_trap;
	record_from_solution = recordPO,
	normC = BK.norminf,
	δp = 0.001)

f, ax = BK.plot(br_po_trap, dash_unstable_style = true)
f
```

The normal form of the second branch point only returns the eigenvectors:

```@example TUTFREIRE
bp_trap = BK.get_normal_form(br_po_trap, 2; detailed = Val(true))
```

Since automatic branch switching is not available for `Trapeze`, we provide the amplitude and parameter increment manually:

```@example TUTFREIRE
br_po_trap_bp = BK.continuation(deepcopy(br_po_trap), 2;
	δp = -0.001,
	ampfactor = 0.01,
	normal_form_kwargs = (detailed = Val(false), ),
	record_from_solution = recordPO,
	normC = BK.norminf)

f, ax = BK.plot(br_po_trap, br_po_trap_bp; branchlabel = ["BP branch", "manual BP switching"], dash_unstable_style = true)
f
```

The bifurcated branch has a Neimark-Sacker point and a period-doubling point. The normal form of the period-doubling point is available but its coefficients are not yet implemented for `Trapeze`:

```@example TUTFREIRE
BK.get_normal_form(br_po_trap_bp, 2; detailed = Val(true))
```

```@example TUTFREIRE
br_po_trap_pd = BK.continuation(deepcopy(br_po_trap_bp), 2;
	δp = -0.001,
	ampfactor = 0.01,
	normal_form_kwargs = (detailed = Val(false), ),
	record_from_solution = recordPO,
	normC = BK.norminf)

f, ax = BK.plot(br_po_trap_pd, br_po_trap_bp; dash_unstable_style = true)
f
```

## Automatic bifurcation diagram of periodic orbits

Finally, we compute the automatic bifurcation diagram of periodic orbits starting from `br_po`. The function `bifurcationdiagram` recursively branches at branch points and period-doubling points. We use a filter to avoid branches with a too large period.

```@example TUTFREIRE
function filter_bp(_, _, pt)
	return BK.saved_solution(pt.x)[end] < 20 # limit period
end

bd_po = @time BK.bifurcationdiagram(deepcopy(br_po.prob), deepcopy(br_po), 4,
	BK.ContinuationPar(BK.getcontparams(br_po2);
		dsmax = 0.005, p_min = -0.7, max_steps = 100);
	normC = BK.norminf,
	bothside = Val(true),
	verbosediagram = true,
	δp = 0.000125,
	normal_form_kwargs = (prm = Val(false),),
	filter_bif_points = filter_bp,
	callback_newton = BK.cbMaxNorm(10),
	linear_algo = BK.COPBLS(),
	)

nothing #hide
```

```@example TUTFREIRE
f, ax = BK.plot(bd_po, dash_unstable_style = true)
f
```
