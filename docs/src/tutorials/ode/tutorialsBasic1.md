# 🟢 Neural mass equation

The following model is taken from [^Cortes]:

$$\left\{\begin{array}{l}
\tau \dot{E}=-E+g\left(J u x E+E_{0}\right) \\
\dot{x}=\tau_{D}^{-1}(1-x)-u E x \\
\dot{u}=U_0 E(1-u)-\tau_{F}^{-1}(u-U_0)
\end{array}\right.$$

with

$$g(y):=\alpha\log(1+exp(y/\alpha)).$$

We use this model as a mean to introduce the basics of `BifurcationKit.jl`, namely the continuation of equilibria.

It is easy to encode the ODE as follows

```@example TUTODE0
using Revise, Plots
import BifurcationKit as BK
import BifurcationKit: @optic

# vector field
function TMvf(z, p)
	(;J, α, E0, τ, τD, τF, U0) = p
	E, x, u = z
	SS0 = J * u * x * E + E0
	SS1 = α * log(1 + exp(SS0 / α))
	[
		(-E + SS1) / τ,
		(1.0 - x) / τD - u * x * E,
		(U0 - u) / τF +  U0 * (1.0 - u) * E
	]
end

# parameter values
par_tm = (α = 1.5, τ = 0.013, J = 3.07, E0 = -2.0, τD = 0.200, U0 = 0.3, τF = 1.5, τS = 0.007)

# initial condition
z0 = [0.238616, 0.982747, 0.367876]

# Bifurcation Problem
prob = BK.ODEBifProblem(TMvf, z0, par_tm, (@optic _.E0);
	record_from_solution = (x, p; k...) -> (E = x[1], x = x[2], u = x[3]),)

nothing #hide
```

We first compute the branch of equilibria 

```@example TUTODE0
# continuation options, we limit the parameter range for E0
opts_br = BK.ContinuationPar(p_min = -4.0, p_max = -0.9)

# continuation of equilibria
br = BK.continuation(prob, BK.PALC(), opts_br;
	# we want to compute both sides of the branch of the initial
	# value of E0 = -2
	bothside = true)

scene = plot(br, legend=:topleft)
```

With detailed information:

```@example TUTODE0
br
```

If you only  want to compute the branch without the bifurcations (more information is provided [Detection of bifurcation points of Equilibria](@ref) ), change the continuation options to

```@example TUTODE0
opts_br = BK.ContinuationPar(p_min = -4.0, p_max = -0.9,
	detect_bifurcation = 0)
	
# continuation of equilibria
br = BK.continuation(prob, BK.PALC(), opts_br;
	# we want to compute both sides of the branch of the initial
	# value of E0 = -2
	bothside = true)

scene = plot(br, plotfold=true)
```

## References

[^Cortes]:> Cortes, Jesus M., Mathieu Desroches, Serafim Rodrigues, Romain Veltz, Miguel A. Muñoz, and Terrence J. Sejnowski. **Short-Term Synaptic Plasticity in the Deterministic Tsodyks–Markram Model Leads to Unpredictable Network Dynamics.**” Proceedings of the National Academy of Sciences 110, no. 41 (October 8, 2013): 16610–15. https://doi.org/10.1073/pnas.1316071110.
