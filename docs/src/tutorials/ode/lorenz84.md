# [Extended Lorenz-84 model (codim 2 + BT/ZH aBS)](@id lorenz)


```@contents
Pages = ["lorenz84.md"]
Depth = 3
```

In this tutorial, we study the extended Lorenz-84 model which is also treated in MatCont [^Kuznetsov]. This model is interesting because it features all codim 2 bifurcations. It is thus convenient to test our algorithms.

After this tutorial, you will be able to
- detect codim 1 bifurcation Fold / Hopf / Branch point
- follow Fold / Hopf points and detect codim 2 bifurcation points
- branch from the detected codim 2 points to curves of Fold / Hopf points (This part is still "work in progress")

The model is as follows

$$\left\{\begin{array}{l}
\dot{X}=-Y^{2}-Z^{2}-\alpha X+\alpha F-\gamma U^{2} \\
\dot{Y}=X Y-\beta X Z-Y+G \\
\dot{Z}=\beta X Y+X Z-Z \\
\dot{U}=-\delta U+\gamma U X+T
\end{array}\right.\tag{E}$$

We start with some imports that are useful in the following.

```@example LORENZ84
using Revise, ForwardDiff, Parameters, Setfield, Plots, LinearAlgebra
using BifurcationKit
const BK = BifurcationKit

# sup norm
norminf(x) = norm(x, Inf)
nothing #hide
```

## Problem setting
We can now encode the vector field (E) in a function and use automatic differentiation to compute its various derivatives.

```@example LORENZ84
# vector field
function Lor(u, p)
    @unpack α,β,γ,δ,G,F,T = p
	X,Y,Z,U = u
	[
		-Y^2 - Z^2 - α*X + α*F - γ*U^2,
		X*Y - β*X*Z - Y + G,
		β*X*Y + X*Z - Z,
		-δ*U + γ*U*X + T
	]
end

# parameter values
parlor = (α = 1//4, β = 1, G = .25, δ = 1.04, γ = 0.987, F = 1.7620532879639, T = .0001265)

# initial condition
z0 =  [2.9787004394953343, -0.03868302503393752,  0.058232737694740085, -0.02105288273117459]
nothing #hide
```

## Continuation and codim 1 bifurcations

Once the problem is set up, we can continue the state w.r.t. $F$ to and detect codim 1 bifurcations. This is achieved as follows:

```@example LORENZ84
# bifurcation problem
recordFromSolutionLor(x, p) = (X = x[1], Y = x[2], Z = x[3], U = x[4])
prob = BifurcationProblem(Lor, z0, setproperties(parlor; T=0.04, F=3.), (@lens _.F);
    recordFromSolution = recordFromSolutionLor)

# continuation options
opts_br = ContinuationPar(pMin = -1.5, pMax = 3.0, ds = 0.002, dsmax = 0.15,
	# Optional: bisection options for locating bifurcations
	nInversion = 6, maxBisectionSteps = 25,
	# number of eigenvalues
	nev = 4, maxSteps = 200)

# compute the branch of solutions
br = @time continuation(prob, PALC(), opts_br;
	normC = norminf,
	bothside = true)

scene = plot(br, plotfold=false, markersize=4, legend=:topleft)
```

With detailed information:

```@example LORENZ84
br
```

## Continuation of Fold points

We follow the Fold points in the parameter plane $(T,F)$. We tell the solver to consider `br.specialpoint[5]` and continue it.

```@example LORENZ84
# function to record the current state
sn_codim2 = continuation(br, 5, (@lens _.T), ContinuationPar(opts_br, pMax = 3.2, pMin = -0.1, detectBifurcation = 1, dsmin=1e-5, ds = -0.001, dsmax = 0.005, nInversion = 10, maxSteps = 130, maxBisectionSteps = 55) ; normC = norminf,
	# detection of codim 2 bifurcations with bisection
	detectCodim2Bifurcation = 2,
	# we update the Fold problem at every continuation step
	updateMinAugEveryStep = 1,
	startWithEigen = false,
	# we save the different components for plotting
	recordFromSolution = recordFromSolutionLor,
	)

scene = plot(sn_codim2, vars=(:X, :U), branchlabel = "Folds", ylims=(-0.5, 0.5))
```

with detailed information

```@example LORENZ84
sn_codim2
```

For example, we can compute the following normal form

```@example LORENZ84
getNormalForm(sn_codim2, 1; nev = 4)
```

## Continuation of Hopf points

We follow the Hopf points in the parameter plane $(T,F)$. We tell the solver to consider `br.specialpoint[3]` and continue it.

```@example LORENZ84
hp_codim2_1 = continuation((@set br.alg.tangent = Bordered()), 3, (@lens _.T), ContinuationPar(opts_br, ds = -0.001, dsmax = 0.02, dsmin = 1e-4, nInversion = 6, detectBifurcation = 1) ; normC = norminf,
	# detection of codim 2 bifurcations with bisection
	detectCodim2Bifurcation = 2,
	# we update the Fold problem at every continuation step
	updateMinAugEveryStep = 1,
	# we save the different components for plotting
	recordFromSolution = recordFromSolutionLor,
	# compute both sides of the initial condition
	bothside = true,
	)

plot(sn_codim2, vars=(:X, :U), branchlabel = "Folds")
	plot!(hp_codim2_1, vars=(:X, :U), branchlabel = "Hopfs")
	ylims!(-0.7,0.7);xlims!(1,1.3)
```

```@example LORENZ84
hp_codim2_1
```

For example, we can compute the following normal form

```@example LORENZ84
getNormalForm(hp_codim2_1, 3; nev = 4)
```

## Continuation of Hopf points from the Bogdanov-Takens point

When we computed the curve of Fold points, we detected a Bogdanov-Takens bifurcation. We can branch from it to get the curve of Hopf points. This is done as follows:

```@example LORENZ84
hp_from_bt = continuation((@set sn_codim2.alg.tangent = Bordered()), 4, ContinuationPar(opts_br, ds = -0.001, dsmax = 0.02, dsmin = 1e-4,
	nInversion = 6, detectBifurcation = 1) ; normC = norminf,
	# detection of codim 2 bifurcations with bisection
	detectCodim2Bifurcation = 2,
	# we update the Fold problem at every continuation step
	updateMinAugEveryStep = 1,
	# we save the different components for plotting
	recordFromSolution = recordFromSolutionLor,
	)

plot(sn_codim2, vars=(:X, :U), branchlabel = "SN")
	plot!(hp_codim2_1, vars=(:X, :U), branchlabel = "Hopf1")
	plot!(hp_from_bt, vars=(:X, :U), branchlabel = "Hopf2")
	ylims!(-0.7,0.75);xlims!(0.95,1.3)
```

with detailed information

```@example LORENZ84
hp_from_bt
```

## Continuation of Hopf points from the Zero-Hopf point

When we computed the curve of Fold points, we detected a Zero-Hopf bifurcation. We can branch from it to get the curve of Hopf points. This is done as follows:

```@example LORENZ84
hp_from_zh = continuation((@set sn_codim2.alg.tangent = Bordered()), 2, ContinuationPar(opts_br, ds = 0.001, dsmax = 0.02, dsmin = 1e-4, nInversion = 6, detectBifurcation = 1, maxSteps = 150) ;
	normC = norminf,
	detectCodim2Bifurcation = 2,
	updateMinAugEveryStep = 1,
	startWithEigen = true,
	recordFromSolution = recordFromSolutionLor,
	bothside = false,
	bdlinsolver = MatrixBLS(),
	)

plot(sn_codim2,vars=(:X, :U),)
	plot!(hp_codim2_1, vars=(:X, :U), branchlabel = "Hopf")
	plot!(hp_from_bt, vars=(:X, :U),  branchlabel = "Hopf2")
	plot!( hp_from_zh, vars=(:X, :U), branchlabel = "Hopf", plotspecialpoints = false, legend = :topleft)
	ylims!(-0.7,0.75);xlims!(0.95,1.3)
```

with detailed information

```@example LORENZ84
hp_from_zh
```

## References

[^Kuznetsov]:> Kuznetsov, Yu A., H. G. E. Meijer, W. Govaerts, and B. Sautois. “Switching to Nonhyperbolic Cycles from Codim 2 Bifurcations of Equilibria in ODEs.” Physica D: Nonlinear Phenomena 237, no. 23 (December 2008): 3061–68.
