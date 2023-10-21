# 🟠 Lorenz-84 model, take 2.


```@contents
Pages = ["lorenz84-PO.md"]
Depth = 3
```

In this tutorial, we study the extended Lorenz-84 model which is also treated in MatCont [^Kuznetsov]. We use this model to showcase the automatic branch switching procedure to
- Fold of periodic orbits from Bautin bifurcation point
- NS of periodic orbits from ZH bifurcation point
- NS of periodic orbits from HH bifurcation point.

> As this model has been studied in this [tutorial](@ref lorenz), we do not give much details and refer to the corresponding tutorial to get the 2 parameters curves of Hopf / Fold bifurcations. 

The model is as follows

$$\left\{\begin{array}{l}
\dot{X}=-Y^{2}-Z^{2}-\alpha X+\alpha F-\gamma U^{2} \\
\dot{Y}=X Y-\beta X Z-Y+G \\
\dot{Z}=\beta X Y+X Z-Z \\
\dot{U}=-\delta U+\gamma U X+T
\end{array}\right.\tag{E}$$

We recall the problem setting:

```@example LORENZ84V2
using Revise, ForwardDiff, Parameters, Plots, LinearAlgebra
using BifurcationKit
const BK = BifurcationKit

# vector field
function Lor(u, p, t = 0)
	@unpack α,β,γ,δ,G,F,T = p
	X,Y,Z,U = u
	[
		-Y^2 - Z^2 - α*X + α*F - γ*U^2,
		X*Y - β*X*Z - Y + G,
		β*X*Y + X*Z - Z,
		-δ*U + γ*U*X + T
	]
end

parlor = (α = 1//4, β = 1., G = .25, δ = 1.04, γ = 0.987, F = 1.7620532879639, T = .0001265)

z0 =  [2.9787004394953343, -0.03868302503393752,  0.058232737694740085, -0.02105288273117459]

recordFromSolutionLor(x, p) = (u = BK.getVec(x);(X = u[1], Y = u[2], Z = u[3], U = u[4]))
prob = BK.BifurcationProblem(Lor, z0, parlor, (@lens _.F);
	record_from_solution = (x, p) -> (X = x[1], Y = x[2], Z = x[3], U = x[4]),)

opts_br = ContinuationPar(p_min = -1.5, p_max = 3.0, ds = 0.002, dsmax = 0.05, n_inversion = 6, detect_bifurcation = 3, max_bisection_steps = 25, nev = 4, max_steps = 200, plot_every_step = 30)
	@set! opts_br.newton_options.verbose = false
	@set! opts_br.newton_options.tol = 1e-12
	br = @time continuation(re_make(prob, params = setproperties(parlor;T=0.04,F=3.)),
	 	PALC(), opts_br;
		normC = norminf, bothside = true)

scene = plot(br, plotfold=false, markersize=4, legend=:topleft)
```

## Two parameters curves of Fold / Hopf bifurcation

We follow the Fold points in the parameter plane $(T,F)$. We tell the solver to consider `br.specialpoint[5]` and continue it.

```@example LORENZ84V2
sn_codim2 = continuation(br, 5, (@lens _.T), ContinuationPar(opts_br, p_max = 3.2, p_min = -0.1, detect_bifurcation = 1, dsmin=1e-5, ds = -0.001, dsmax = 0.005, n_inversion = 10, save_sol_every_step = 1, max_steps = 130, max_bisection_steps = 55) ; plot = true,
	verbosity = 0,
	normC = norminf,
	detect_codim2_bifurcation = 2,
	update_minaug_every_step = 1,
	start_with_eigen = false,
	bothside = false,
	)

hp_codim2_1 = continuation(br, 3, (@lens _.T), ContinuationPar(opts_br, ds = -0.001, dsmax = 0.02, dsmin = 1e-4, n_inversion = 8, save_sol_every_step = 1, detect_bifurcation = 1) ; plot = false, verbosity = 0,
	normC = norminf,
	# tangentAlgo = BorderedPred(),
	detect_codim2_bifurcation = 2,
	update_minaug_every_step = 1,
	start_with_eigen = true,
	bothside = true,
	)

plot(sn_codim2, vars=(:F, :T), branchlabel = "SN")
plot!(hp_codim2_1, vars=(:F, :T), branchlabel = "Hopf1", xlims = (1,2.7), ylims = (-0.06,0.06))
```

## Fold bifurcations of periodic orbits from Bautin bifurcation

We compute the branch of Fold of periodic orbits from the Bautin bifurcation (labelled `:gh`) in the previous figure. In this tutorial, we focus on orthogonal collocation but standard shooting would do too.

```@example LORENZ84V2
opts_fold_po = ContinuationPar(hp_codim2_1.contparams, dsmax = 0.01, detect_bifurcation = 0, max_steps = 30, detect_event = 0, ds = 0.001, plot_every_step = 10, a = 0.8)
@set! opts_fold_po.newton_options.verbose = false
@set! opts_fold_po.newton_options.tol = 1e-8
fold_po = continuation(hp_codim2_1, 3, opts_fold_po, 
		PeriodicOrbitOCollProblem(20, 3, meshadapt = false);
		normC = norminf,
		δp = 0.02,
		update_minaug_every_step = 0,
		jacobian_ma = :minaug,
		verbosity = 0, plot = false,
	)
plot!(fold_po, vars=(:F, :T), branchlabel = "Fold-PO")
```

## NS bifurcations of periodic orbits from Hopf-Hopf bifurcation

When we computed the curve of Hopf points, we detected a Hopf-Hopf bifurcation. We can branch from it to get the curve of NS points. This is done as follows:

```@example LORENZ84V2
opts_ns_po = ContinuationPar(hp_codim2_1.contparams, dsmax = 0.02, detect_bifurcation = 1, max_steps = 20, ds = -0.001, detect_event = 0)
@set! opts_ns_po.newton_options.verbose = false
@set! opts_ns_po.newton_options.tol = 1e-9
@set! opts_ns_po.newton_options.max_iterations = 10
ns_po1 = continuation(hp_codim2_1, 4, opts_ns_po, 
		PeriodicOrbitOCollProblem(20, 3, update_section_every_step = 1);
		detect_codim2_bifurcation = 0,
		normC = norminf,
		δp = 0.02,
		update_minaug_every_step = 1,
		# which of the 2 NS curves should we compute?
		whichns = 1,
		jacobian_ma = :minaug,
		verbosity = 3,
		)
plot!(ns_po1, vars=(:F, :T), branchlabel = "NS1")
```

```@example LORENZ84V2
ns_po2 = continuation(hp_codim2_1, 4, opts_ns_po, 
		PeriodicOrbitOCollProblem(30, 3, update_section_every_step = 1);
		detect_codim2_bifurcation = 0,
		normC = norminf,
		δp = 0.02,
		update_minaug_every_step = 1,
		# which of the 2 NS curves should we compute?
		whichns = 2,
		jacobian_ma = :minaug,
		)
plot!(ns_po2, vars=(:F, :T), branchlabel = "NS2")
```

## References

[^Kuznetsov]:> Kuznetsov, Yu A., H. G. E. Meijer, W. Govaerts, and B. Sautois. “Switching to Nonhyperbolic Cycles from Codim 2 Bifurcations of Equilibria in ODEs.” Physica D: Nonlinear Phenomena 237, no. 23 (December 2008): 3061–68.
