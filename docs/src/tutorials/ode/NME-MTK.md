# 🟢 Neural mass equation - MTK


```@contents
Pages = ["NME-MTK.md"]
Depth = 3
```

In this tutorial, we study the following model taken from [^Cortes]. It is essentially the same tutorial as in [Neural mass equation](@ref nmepo) but treated with [ModelingToolkit.jl](https://mtk.sciml.ai/stable/).

$$\left\{\begin{array}{l}
\tau \dot{E}=-E+g\left(J u x E+E_{0}\right) \\
\dot{x}=\tau_{D}^{-1}(1-x)-u E x \\
\dot{u}=U E(1-u)-\tau_{F}^{-1}(u-U)
\end{array}\right.$$


The model is interesting because the branch of periodic solutions converges to an homoclinic orbit which is challenging to compute with our methods.

It is easy to encode the ODE as follows

```@example TUTNMEMTK
using Revise
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using DifferentialEquations, Plots
using BifurcationKit
const BK = BifurcationKit

@variables E(t) x(t) u(t) SS0(t) SS1(t) 	# independent and dependent variables
@parameters U0 τ J E0 τD U0 τF α    		# parameters

# define the model
eqs = [ D(E) ~ (-E + α * log(1 + exp((J * u * x * E + E0) / α))) / τ,
            D(x) ~ (1.0 - x) / τD - u * x * E,
            D(u) ~ (U0 - u) / τF +  U0 * (1.0 - u) * E]
@mtkbuild nsys = ODESystem(eqs, t)
u0_guess = [E => 0.238616, x => 0.982747, u => 0.367876]
p_start = Dict(α => 1.4, τ => 0.013, J => 3.07, E0 => -2.0, τD => 0.20, U0 => 0.3, τF => 1.5)
bif_par = E0
plot_var = E

prob = BifurcationProblem(nsys,
    u0_guess,
    p_start,
    bif_par;
    plot_var = plot_var,
    jac = false)
```

We first compute the branch of equilibria

```@example TUTNMEMTK
# continuation options
opts_br = ContinuationPar(p_min = -10.0, p_max = -0.9,
	# parameters to have a smooth result
	ds = 0.04, dsmax = 0.05,
	nev = 3)

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
optn_po = NewtonPar(tol = 1e-8,  max_iterations = 10)

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.15, ds= -0.0001, dsmin = 1e-4, p_max = 0., p_min=-5.,
	max_steps = 150, newton_options = optn_po,
	nev = 3, plot_every_step = 10, detect_bifurcation = 0)

# arguments for periodic orbits
# this is mainly for printing purposes
args_po = (	record_from_solution = (x, p; k...) -> begin
		xtt = get_periodic_orbit(p.prob, x, p.p)
		return (max = maximum(xtt[1,:]),
				min = minimum(xtt[1,:]),
				period = getperiod(p.prob, x, p.p))
	end,
	plot_solution = (x, p; k...) -> begin
		xtt = get_periodic_orbit(p.prob, x, p.p)
		plot!(xtt.t, xtt[1,:]; label = "E", k...)
		plot!(xtt.t, xtt[2,:]; label = "x", k...)
		plot!(xtt.t, xtt[3,:]; label = "u", k...)
		plot!(br; subplot = 1, putspecialptlegend = false)
		end,
	normC = norminf)

br_pocoll = @time continuation(
	# we want to branch form the 4th bif. point
	br, 4, opts_po_cont,
	# we want to use the Collocation method to locate PO, with polynomial degree 5
	PeriodicOrbitOCollProblem(30, 5; meshadapt = true);
	# regular continuation options
	args_po..., callback_newton = BK.cbMaxNorm(10))

scene = plot(br, br_pocoll, markersize = 3)
plot!(scene, br_pocoll.param, br_pocoll.min, label = "")
```

We plot the maximum (resp. minimum) of the limit cycle. We can see that the min converges to the smallest equilibrium indicating a homoclinic orbit.

### Plot of some of the periodic orbits as function of $E_0$

We can plot some of the previously computed periodic orbits in the plane $(E,x)$ as function of $E_0$:

```@example TUTNMEMTK
plot()
# fetch the saved solutions
for i_sol in 1:2:40
	# periodic orbit
	traj = BK.get_periodic_orbit(br_pocoll, i_sol)
	plot!(traj[1,:], traj[2,:], xlabel = "E", ylabel = "x", label = "")
end
title!("")
```

## References

[^Cortes]:> Cortes, Jesus M., Mathieu Desroches, Serafim Rodrigues, Romain Veltz, Miguel A. Muñoz, and Terrence J. Sejnowski. **Short-Term Synaptic Plasticity in the Deterministic Tsodyks–Markram Model Leads to Unpredictable Network Dynamics.**” Proceedings of the National Academy of Sciences 110, no. 41 (October 8, 2013): 16610–15. https://doi.org/10.1073/pnas.1316071110.
