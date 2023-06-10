# Neural mass equation - MTK


```@contents
Pages = ["NME-MTK.md"]
Depth = 3
```

In this tutorial, we study the following model taken from [^Cortes]. It is essentially the same tutorial as in [Neural mass equation (Hopf aBS)](@ref) but treated with [ModelingToolkit.jl](https://mtk.sciml.ai/stable/).

$$\left\{\begin{array}{l}
\tau \dot{E}=-E+g\left(J u x E+E_{0}\right) \\
\dot{x}=\tau_{D}^{-1}(1-x)-u E x \\
\dot{u}=U E(1-u)-\tau_{F}^{-1}(u-U)
\end{array}\right.$$


The model is interesting because the branch of periodic solutions converges to an homoclinic orbit which is challenging to compute with our methods.

It is easy to encode the ODE as follows

```@example TUTNMEMTK
using Revise, ModelingToolkit, LinearAlgebra
using DifferentialEquations, Plots
using BifurcationKit
const BK = BifurcationKit

norminf(x) = norm(x, Inf)
indexof(sym, syms) = findfirst(isequal(sym),syms)

@variables t E(t) x(t) u(t) SS0(t) SS1(t) 	# independent and dependent variables
@parameters U0 τ J E0 τD U0 τF τS α    		# parameters
D = Differential(t) 				# define an operator for the differentiation w.r.t. time

# define the model
@named NMmodel = ODESystem([SS0 ~ J * u * x * E + E0,
	SS1 ~ α * log(1 + exp(SS0 / α)),
	D(E) ~ (-E + SS1) / τ,
	D(x) ~ (1.0 - x) / τD - u * x * E,
	D(u) ~ (U0 - u) / τF +  U0 * (1 - u) * E],
	defaults = Dict(E => 0.238616, x => 0.982747, u => 0.367876,
	α => 1.5, τ => 0.013, J => 3.07, E0 => -2.0, τD => 0.200, U0 => 0.3, τF => 1.5, τS => 0.007))

# get the vector field and jacobian
odeprob = ODEProblem(structural_simplify(NMmodel), [], (0.,10.), [], jac = true)
odefun = odeprob.f
F = (u,p) -> odefun(u,p,0)
J = (u,p) -> odefun.jac(u,p,0)

id_E0 = indexof(E0, parameters(NMmodel))
par_tm = odeprob.p

# we collect the differentials together in a problem
prob = BifurcationProblem(F, odeprob.u0, par_tm, (@lens _[id_E0]); J = J,
    recordFromSolution = (x, p) -> (E = x[1], x = x[2], u = x[3]))
nothing #hide
```

We first compute the branch of equilibria

```@example TUTNMEMTK
# continuation options
opts_br = ContinuationPar(pMin = -10.0, pMax = -0.9,
	# parameters to have a smooth result
	ds = 0.04, dsmax = 0.05,
	# this is to detect bifurcation points precisely with bisection
	detectBifurcation = 3,
	# Optional: bisection options for locating bifurcations
	nInversion = 8, maxBisectionSteps = 25, nev = 3)

# continuation of equilibria
br = continuation(prob, PALC(tangent = Bordered()), opts_br; normC = norminf)

scene = plot(br, plotfold=false, markersize=3, legend=:topleft)
```

With detailed information:

```@example TUTNMEMTK
br
```

## Branch of periodic orbits with Collocation method

We then compute the branch of periodic orbits from the last Hopf bifurcation point (on the right). We use finite differences to discretize the problem of finding periodic orbits. Obviously, this will be problematic when the period of the limit cycle grows unbounded close to the homoclinic orbit.

```@example TUTNMEMTK
# newton parameters
optn_po = NewtonPar(verbose = true, tol = 1e-8,  maxIter = 10)

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.1, ds= -0.0001, dsmin = 1e-4, pMax = 0., pMin=-5.,
	maxSteps = 110, newtonOptions = (@set optn_po.tol = 1e-7),
	nev = 3, plotEveryStep = 10, detectBifurcation = 0)

# arguments for periodic orbits
# this is mainly for printing purposes
args_po = (	recordFromSolution = (x, p) -> begin
		xtt = BK.getPeriodicOrbit(p.prob, x, p.p)
		return (max = maximum(xtt[1,:]),
				min = minimum(xtt[1,:]),
				period = getPeriod(p.prob, x, p.p))
	end,
	plotSolution = (x, p; k...) -> begin
		xtt = BK.getPeriodicOrbit(p.prob, x, p.p)
		plot!(xtt.t, xtt[1,:]; label = "E", k...)
		plot!(xtt.t, xtt[2,:]; label = "x", k...)
		plot!(xtt.t, xtt[3,:]; label = "u", k...)
		plot!(br; subplot = 1, putspecialptlegend = false)
		end,
	normC = norminf)


Mt = 30 # number of time sections
	br_pocoll = @time continuation(
	# we want to branch form the 4th bif. point
	br, 4, opts_po_cont,
	# we want to use the Collocation method to locate PO, with polynomial degree 5
	PeriodicOrbitOCollProblem(Mt, 5);
	# regular continuation options
	args_po..., callbackN = BK.cbMaxNorm(10))

scene = plot(br, br_pocoll, markersize = 3)
plot!(scene, br_pocoll.param, br_pocoll.min, label = "")
```

We plot the maximum (resp. minimum) of the limit cycle. We can see that the min converges to the smallest equilibrium indicating a homoclinic orbit.

### Plot of some of the periodic orbits as function of $E_0$

We can plot some of the previously computed periodic orbits in the plane $(E,x)$ as function of $E_0$:

```@example TUTNMEMTK
plot()
# fetch the saved solutions
for sol in br_pocoll.sol[1:2:40]
	# periodic orbit
	po = sol.x
	# get the mesh and trajectory
	traj = BK.getPeriodicOrbit(br_pocoll.prob, po, @set par_tm[id_E0] = sol.p)
	plot!(traj[1,:], traj[2,:], xlabel = "E", ylabel = "x", label = "")
end
title!("")
```

## References

[^Cortes]:> Cortes, Jesus M., Mathieu Desroches, Serafim Rodrigues, Romain Veltz, Miguel A. Muñoz, and Terrence J. Sejnowski. **Short-Term Synaptic Plasticity in the Deterministic Tsodyks–Markram Model Leads to Unpredictable Network Dynamics.**” Proceedings of the National Academy of Sciences 110, no. 41 (October 8, 2013): 16610–15. https://doi.org/10.1073/pnas.1316071110.
