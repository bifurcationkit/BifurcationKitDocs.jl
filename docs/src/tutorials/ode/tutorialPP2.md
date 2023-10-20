# 🟢 pp2 example from AUTO07p (aBD + Hopf aBS)

```@contents
Pages = ["tutorialsPP2.md"]
Depth = 3
```

The goal of this example is to show how to use automatic bifurcation diagram computation for a simple ODE.

The following equations are a model of type predator-prey. The example is taken from Auto07p:

$$\begin{array}{l}
u_{1}^{\prime}=3 u_{1}\left(1-u_{1}\right)-u_{1} u_{2}-p_1\left(1-e^{-5 u_{1}}\right) \\
u_{2}^{\prime}=-u_{2}+3 u_{1} u_{2}
\end{array}$$

It is easy to encode the ODE as follows

```@example TUTPP2
using Revise, Parameters, Plots
using BifurcationKit

# function to record information from a solution
recordFromSolution(x, p) = (u1 = x[1], u2 = x[2])

function pp2!(dz, z, p, t = 0)
	@unpack p1, p2, p3, p4 = p
	u1, u2 = z
	dz[1] = p2 * u1 * (1 - u1) - u1 * u2 - p1 * (1 - exp(-p3 * u1))
	dz[2] =	-u2 + p4 * u1 * u2
	dz
end

# parameters of the model
par_pp2 = (p1 = 1., p2 = 3., p3 = 5., p4 = 3.)

# initial condition
z0 = zeros(2)

# bifurcation problem
prob = BifurcationProblem(pp2!, z0, par_pp2,
	# specify the continuation parameter
	(@lens _.p1), record_from_solution = recordFromSolution)

nothing #hide
```

## Automatic bifurcation diagram computation

We set up the options or the continuation

```@example TUTPP2
# continuation options
opts_br = ContinuationPar(p_min = 0.1, p_max = 1.0, dsmax = 0.01,
	# options to detect bifurcations
	detect_bifurcation = 3, n_inversion = 8, max_bisection_steps = 25,
	# number of eigenvalues
	nev = 2,
	# maximum number of continuation steps
	max_steps = 1000,)

nothing #hide
```

We are now ready to compute the diagram

```@example TUTPP2
diagram = bifurcationdiagram(prob, PALC(),
	# very important parameter. It specifies the maximum amount of recursion
	# when computing the bifurcation diagram. It means we allow computing branches of branches of branches
	# at most in the present case.
	3,
	(args...) -> setproperties(opts_br; ds = -0.001, dsmax = 0.01, n_inversion = 8, detect_bifurcation = 3);
	# δp = -0.01,
	verbosity = 0, plot = false)

scene = plot(diagram; code = (), title="$(size(diagram)) branches", legend = false)
```


## Branch of periodic orbits with finite differences

As you can see on the diagram, there is a Hopf bifurcation indicated by a red dot.  Let us compute the periodic orbit branching from the Hopf point.

We first find the branch

```@example TUTPP2
# branch of the diagram with Hopf point
brH = get_branch(diagram, (2,1)).γ

# newton parameters
optn_po = NewtonPar(tol = 1e-8,  max_iterations = 25)

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.1, ds= -0.001, dsmin = 1e-4,
 newton_options = (@set optn_po.tol = 1e-8), tol_stability = 1e-2,
 detect_bifurcation = 1)

Mt = 101 # number of time sections
	br_po = continuation(
	brH, 1, opts_po_cont,
	PeriodicOrbitTrapProblem(M = Mt;
	    # specific linear solver for ODEs
    	jacobian = :Dense);
	record_from_solution = (x, p) -> (xtt=reshape(x[1:end-1],2,Mt);
		return (max = maximum(xtt[1,:]),
			min = minimum(xtt[1,:]),
			period = x[end])),
	finalise_solution = (z, tau, step, contResult; prob = nothing, kwargs...) -> begin
		# limit the period
		z.u[end] < 100
		end,
	normC = norminf)


plot(diagram); plot!(br_po, label = "Periodic orbits", legend = :bottomright)
```

Let us now plot an orbit

```@example TUTPP2
# extract the different components
orbit = get_periodic_orbit(br_po, 10)
plot(orbit.t, orbit.u[1,:]; label = "u1", markersize = 2)
plot!(orbit.t, orbit.u[2,:]; label = "u2", xlabel = "time", title = "period = $(orbit.t[end])")
```
