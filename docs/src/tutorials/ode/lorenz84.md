# Extended Lorenz-84 model (codim 2 + BT/ZH aBS)


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

# we group the differentials together
jet = BK.getJet(Lor;matrixfree=false)

# parameter values
parlor = (α = 1//4, β = 1, G = .25, δ = 1.04, γ = 0.987, F = 1.7620532879639, T = .0001265)

# initial condition
z0 =  [2.9787004394953343, -0.03868302503393752,  0.058232737694740085, -0.02105288273117459]
nothing #hide
```

## Continuation and codim 1 bifurcations

Once the problem is set up, we can continue the state w.r.t. $F$ to and detect codim 1 bifurcations. This is achieved as follows:

```@example LORENZ84
# continuation options
opts_br = ContinuationPar(pMin = -1.5, pMax = 3.0, ds = 0.002, dsmax = 0.15,
	# options to detect codim 1 bifurcations using bisection
	detectBifurcation = 3,
	# Optional: bisection options for locating bifurcations
	nInversion = 6, maxBisectionSteps = 25,
	# number of eigenvalues
	nev = 4, maxSteps = 200)

# compute the branch of solutions
br, = @time continuation(jet[1], jet[2], z0, setproperties(parlor; T=0.04, F=3.), (@lens _.F), opts_br;
	recordFromSolution = (x, p) -> (X = x[1], Y = x[2], Z = x[3], U = x[4]),
	normC = norminf,
	bothside = true)

scene = plot(br, plotfold=false, markersize=4, legend=:topleft)
```

With detailed information:

```@example LORENZ84
br
```

## Continuation of Fold points

We follow the Fold points in the parameter plane $(T,F)$. We tell the solver to consider `br.specialpoint[4]` and continue it.

```@example LORENZ84
# function to record the current state
recordFromSolutionLor(x, p) = ((X = x.u[1], Y = x.u[2], Z = x.u[3], U = x.u[4]))

sn_codim2, = continuation(jet[1:2]..., br, 4, (@lens _.T), ContinuationPar(opts_br, pMax = 3.2, pMin = -0.1, detectBifurcation = 1, dsmin=1e-5, ds = -0.001, dsmax = 0.005, nInversion = 10, saveSolEveryStep = 1, maxSteps = 130, maxBisectionSteps = 55) ; normC = norminf,
	# detection of codim 2 bifurcations with bisection
	detectCodim2Bifurcation = 2,
	# we update the Fold problem at every continuation step
	updateMinAugEveryStep = 1,
	startWithEigen = false,
	# we save the different components for plotting
	recordFromSolution = recordFromSolutionLor,
	# give analytic higher differentials, useful for normal forms
	d2F = jet[3], d3F = jet[4],
	)

scene = plot(sn_codim2, vars=(:X, :U), branchlabel = "Folds", ylims=(-0.5, 0.5))
```

with detailed information

```@example LORENZ84
sn_codim2
```

For example, we can compute the following normal form

```@example LORENZ84
computeNormalForm(jet..., sn_codim2, 1; nev = 4)
```

## Continuation of Hopf points

We follow the Hopf points in the parameter plane $(T,F)$. We tell the solver to consider `br.specialpoint[2]` and continue it.

```@example LORENZ84
hp_codim2_1, = continuation(jet[1:2]..., br, 2, (@lens _.T), ContinuationPar(opts_br, ds = -0.001, dsmax = 0.02, dsmin = 1e-4, nInversion = 6, saveSolEveryStep = 1, detectBifurcation = 1) ; normC = norminf,
	# tangent algorithm
	tangentAlgo = BorderedPred(),
	# detection of codim 2 bifurcations with bisection
	detectCodim2Bifurcation = 2,
	# we update the Fold problem at every continuation step
	updateMinAugEveryStep = 1,
	# we save the different components for plotting
	recordFromSolution = recordFromSolutionLor,
	# give analytic higher differentials, useful for normal forms
	d2F = jet[3], d3F = jet[4],
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
computeNormalForm(jet..., hp_codim2_1, 2; nev = 4)
```

## Continuation of Hopf points from the Bogdanov-Takens bifurcation

When we computed the curve of Fold points, we detected a Bogdanov-Takens bifurcation. We can branch from it to get the curve of Hopf points. This is done as follows:

```@example LORENZ84
hp_from_bt, = continuation(jet..., sn_codim2, 4, ContinuationPar(opts_br, ds = -0.001, dsmax = 0.02, dsmin = 1e-4,
	nInversion = 6, detectBifurcation = 1) ; normC = norminf,
	# tangent algorithm
	tangentAlgo = BorderedPred(),
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

## Continuation of Hopf points from the Zero-Hopf bifurcation 

When we computed the curve of Fold points, we detected a Zero-Hopf bifurcation. We can branch from it to get the curve of Hopf points. This is done as follows:

```@example LORENZ84
hp_from_zh, = continuation(jet..., sn_codim2, 2, ContinuationPar(opts_br, ds = 0.001, dsmax = 0.02, dsmin = 1e-4, nInversion = 6, detectBifurcation = 1, maxSteps = 150) ;
	normC = norminf,
	tangentAlgo = BorderedPred(),
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
