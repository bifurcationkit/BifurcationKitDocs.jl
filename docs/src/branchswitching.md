# [Branch switching](@id Branch-switching-page)

The precise definition of the methods are given in [Branch switching (branch point)](@ref) and [Branch switching (Hopf point)](@ref).

```@contents
Pages = ["branchswitching.md"]
Depth = 2
```

## Summary of branching procedures

We collect in the following table the list of automatic branch switching (aBS) functions. Their detailed explanation follows in this page.

| function | ind-th bif. point | Type `T` | description |
|---|---|---|---|
|  `continuation(br::ContResult{T}, ind::Int; kw...)` | `:bp`, `:nd`| `EquilibriumCont`  |  aBS from equilibria to equilibria  |
|  `continuation(br::ContResult{T}, ind::Int, lens2::Lens; kw...)` | `:bp`, `:hopf`| `EquilibriumCont` | Fold/Hopf continuation w.r.t. parameters `getLens(br)` and `lens2`  |
|  `continuation(br::ContResult{T}, ind::Int; kw...)` | `:bt,:zh,:hh`| ` FoldCont,HopfCont` | switch to Fold/Hopf continuation from Hopf/Fold w.r.t. parameters of codim 2 `br`  |
| `continuation(br::ContResult{T}, ind_hopf::Int, ::ContinuationPar, prob::AbstractPeriodicOrbitProblem)`   | `:hopf` |  `EquilibriumCont` | Branch switching from Hopf point to periodic orbits |
| `continuation(br::ContResult{T}, ind::Int, kw...)`   | `:bp,:pd` |  `PeriodicOrbitCont` | Branch switching from Branch / Period-doubling point of periodic orbits to curve of periodic orbits |
| `continuation(br::ContResult{T}, ind::Int, kw...)`   | `:gh,:zh,:hh` |  `TwoParamCont` | Branch switching from Bautin / Zero-Hopf/ Hopf-Hopf point to curve of Fold/NS of periodic orbits |

## From simple branch point to equilibria

You can perform automatic branch switching by calling `continuation` with the following options:

```julia
continuation(br::ContResult, ind_bif::Int, optionsCont::ContinuationPar; kwargs...)
```

where `br` is a branch computed after a call to `continuation` with detection of bifurcation points enabled. This call computes the branch bifurcating from the `ind_bif `th bifurcation point in `br`. An example of use is provided in [2d generalized Bratu–Gelfand problem](@ref).

> See [Branch switching (branch point)](@ref) precise method definition

### Simple example

```@example TUT1
using BifurcationKit, Setfield, Plots

# vector field of transcritical bifurcation
F(x, p) = [x[1] * (p.μ - x[1])]

# parameters of the vector field
par = (μ = -0.2, )

# problem (automatic differentiation)
prob = BifurcationProblem(F, [0.], par, (@lens _.μ); recordFromSolution = (x, p) -> x[1])

# compute branch of trivial equilibria (=0) and detect a bifurcation point
opts_br = ContinuationPar(detectBifurcation = 3)
br = continuation(prob, PALC(), opts_br)
	
# perform branch switching on one side of the bifurcation point
br1Top = continuation(br, 1, setproperties(opts_br; maxSteps = 14) )

# on the other side
br1Bottom = continuation(br, 1, setproperties(opts_br; ds = -opts_br.ds, maxSteps = 14))

scene = plot(br, br1Top, br1Bottom; branchlabel = ["br", "br1Top", "br1Bottom"], legend = :topleft)
```

## Algorithms
- for the pitchfork bifurcation, the normal form is computed and non-trivial zeros are used to produce guesses for points on the bifurcated branch.


## [From non simple branch point to equilibria](@id abs-simple-eq)

We provide an automatic branch switching method in this case. The method is to first compute the reduced equation (see [Non-simple branch point](@ref)) and use it to compute the nearby solutions. These solutions are seeded as initial guess for [`continuation`](@ref). Hence, you can perform automatic branch switching by calling `continuation` with the following options:

```julia
continuation(br::ContResult, ind_bif::Int, optionsCont::ContinuationPar;
	kwargs...)
```

An example of use is provided in [2d generalized Bratu–Gelfand problem](@ref).

> See [Branch switching (branch point)](@ref) for the precise method definition

## From Hopf point to periodic orbits

In order to compute the bifurcated branch of periodic solutions at a Hopf bifurcation point, you need to choose a method to compute periodic orbits among:

- [Periodic orbits based on Trapezoidal rule](@ref)
- [Periodic orbits based on orthogonal collocation](@ref)
- [Periodic orbits based on the shooting method](@ref)

Once you have decided which method to use, you use the following call:

```julia
continuation(br::ContResult, ind_HOPF::Int, _contParams::ContinuationPar,
	prob::AbstractPeriodicOrbitProblem ;
	δp = nothing, ampfactor = 1, kwargs...)
```

We refer to [`continuation`](@ref) for more information about the arguments. Here, we just say a few words about how we can specify `prob::AbstractPeriodicOrbitProblem`.

- For [Periodic orbits based on Trapezoidal rule](@ref), you can pass `PeriodicOrbitTrapProblem(M = 51)` where `M` is the number of times slices in the periodic orbit.

- For [Periodic orbits based on orthogonal collocation](@ref), you can pass `PeriodicOrbitOCollProblem(M, m)` where `M` is the number of times slices in the periodic orbit and `m` is the degree of the collocation polynomials.

- For [Periodic orbits based on the shooting method](@ref), you need more parameters. For example, you can pass `ShootingProblem(M, odeprob, Euler())` or `PoincareShootingProblem(M, odeprob, Euler())` where `odeprob::ODEProblem` (see [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/types/ode_types/)) is an ODE problem to specify the Cauchy problem amd `M` is the number of sections.

Several examples are provided in [1d Brusselator (automatic)](@ref) or [2d Ginzburg-Landau equation (finite differences, codim 2, Hopf aBS)](@ref cgl).

> See [Branch switching (Hopf point)](@ref) for the precise method definition

## From Branch / Period-doubling point of curve of periodic orbits

We do not provide (for now) the automatic branching procedure for these bifurcations of periodic orbits. As a consequence, the user is asked to provide the amplitude of the bifurcated solution.

We provide the branching method for the following methods to compute periodic orbits: [`PeriodicOrbitTrapProblem`](@ref),[`ShootingProblem`](@ref) and [`PoincareShootingProblem`](@ref). The call is as follows. Please note that a deflation is included in this method to simplify branch switching.

An example of use is provided in [Period doubling in Lur'e problem (PD aBS)](@ref).

```julia
continuation(br::AbstractBranchResult, ind_PD::Int, contParams::ContinuationPar;
	δp = 0.1, ampfactor = 1, usedeflation = false, kwargs...)
```

## From Bogdanov-Takens (BT) point to Fold / Hopf curve

We provide an automatic branch switching method in this case (see for example [Extended Lorenz-84 model (codim 2 + BT/ZH aBS)](@ref lorenz) or [2d Ginzburg-Landau equation (finite differences, codim 2, Hopf aBS)](@ref cgl)). Hence, you can perform automatic branch switching by calling `continuation` with the following options:

```julia
continuation(br::ContResult, ind_BT::Int,
	options_cont::ContinuationPar = br.contparams;
	nev = options_cont.nev,
	detectCodim2Bifurcation::Int = 0,
	startWithEigen = false,
	autodiff = false,
	Teigvec = getvectortype(br),
	scaleζ = norm,
	kwargs...)
```

where `ind_BT` is the index of the BT point in `br`. Note that the BT has been detected during Fold or Hopf continuation. Calling the above method thus switches from Fold continuation to Hopf continuation (and vice-versa) automatically with the same parameter axis.

> Check the docs of [Fold / Hopf Continuation](@ref) and particularly [Setting the jacobian](@ref jac-fold) for improving the speed of computation for large scale systems.

## From Zero-Hopf (ZH) point to Fold / Hopf curve

We provide an automatic branch switching method in this case (see for example [Extended Lorenz-84 model (codim 2 + BT/ZH aBS)](@ref lorenz) or [2d Ginzburg-Landau equation (finite differences, codim 2, Hopf aBS)](@ref cgl)). Hence, you can perform automatic branch switching by calling `continuation` with the following options:

```julia
continuation(br::ContResult, ind_ZH::Int,
	options_cont::ContinuationPar = br.contparams;
	nev = options_cont.nev,
	detectCodim2Bifurcation::Int = 0,
	startWithEigen = false,
	autodiff = false,
	Teigvec = getvectortype(br),
	scaleζ = norm,
	kwargs...)
```

where `ind_ZH` is the index of the ZH point in `br`. Note that the ZH has been detected during Fold or Hopf continuation. Calling the above method thus switches from Fold continuation to Hopf continuation (and vice-versa) automatically with the same parameter axis.

> Check the docs of [Fold / Hopf Continuation](@ref) and particularly [Setting the jacobian](@ref jac-fold) for improving the speed of computation for large scale systems.

## From Hopf-Hopf (HH) point to Fold / Hopf curve

We provide an automatic branch switching method in this case (see for example [Extended Lorenz-84 model (codim 2 + BT/ZH aBS)](@ref lorenz) or [2d Ginzburg-Landau equation (finite differences, codim 2, Hopf aBS)](@ref cgl)). Hence, you can perform automatic branch switching by calling `continuation` with the following options:

```julia
continuation(br::ContResult, ind_HH::Int,
	options_cont::ContinuationPar = br.contparams;
	δp = nothing, ampfactor::Real = 1,
	nev = options_cont.nev,
	detectCodim2Bifurcation::Int = 0,
	startWithEigen = false,
	autodiff = false,
	Teigvec = getvectortype(br),
	scaleζ = norm,
	kwargs...)
```

where `ind_HH` is the index of the HH point in `br`. Note that the HH has been detected during Hopf continuation. Calling the above method thus switches from Hopf continuation to another Hopf branch automatically with the same parameter axis.

> Check the docs of [Fold / Hopf Continuation](@ref) and particularly [Setting the jacobian](@ref jac-fold) for improving the speed of computation for large scale systems.

## From Bautin point to curve Folds of period orbits

From the [Bautin normal form](http://scholarpedia.org/article/Bautin_bifurcation), we know that there is a curve of Fold of periodic orbits which branches from it.

We provide an automatic branch switching method in this case which reads as follows:

```julia
continuation(br::HopfCont, ind_BAUTIN::Int, 
	_contParams::ContinuationPar,
    prob::AbstractPeriodicOrbitProblem ;
    δp = nothing, ampfactor = 1, kwargs...)
``` 

where `prob` is a method to compute periodic orbits (see [From Hopf point to periodic orbits](@ref) for more information).

Note that the two parameters in `br` will be used for the continuation of Fold points of periodic orbits.

## From Zero-Hopf (ZH) point to curve NS of period orbits

From the [Zero-Hopf normal form](http://scholarpedia.org/article/Zero-Hopf_bifurcation), we know that there is a curve of Neimark-Sacker (NS) bifurcations of periodic orbits which branches from it.

We provide an automatic branch switching method in this case which reads as follows:

```julia
continuation(br::TwoParamCont, ind_ZH::Int, 
	_contParams::ContinuationPar,
    prob::AbstractPeriodicOrbitProblem ;
    δp = nothing, ampfactor = 1, kwargs...)
``` 

where `prob` is a method to compute periodic orbits (see [From Hopf point to periodic orbits](@ref) for more information).

Note that the two parameters in `br` will be used for the continuation of NS points of periodic orbits.

## From Hopf-Hopf (HH) point to curve NS of period orbits

From the [Hopf-Hopf normal form](http://scholarpedia.org/article/Hopf-Hopf_bifurcation), we know that there are two curves of Neimark-Sacker (NS) bifurcations of periodic orbits which branches from it.

We provide an automatic branch switching method in this case which reads as follows:

```julia
continuation(br::TwoParamCont, ind_HH::Int, 
	_contParams::ContinuationPar,
    prob::AbstractPeriodicOrbitProblem ;
    δp = nothing, ampfactor = 1, 
    whichns = 1,
    kwargs...)
``` 

where `prob` is a method to compute periodic orbits (see [From Hopf point to periodic orbits](@ref) for more information).
The option `whichns` which belongs to {1,2} controls which NS curve you want to compute. 

Note that the two parameters in `br` will be used for the continuation of NS points of periodic orbits.
	