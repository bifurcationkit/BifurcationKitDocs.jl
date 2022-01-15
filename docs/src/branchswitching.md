# [Branch switching](@id Branch-switching-page)

The precise definition of the methods are given in [Branch switching (branch point)](@ref) and [Branch switching (Hopf point)](@ref).

```@contents
Pages = ["branchswitching.md"]
Depth = 3
```

## Branch switching from simple branch point to equilibria

You can perform automatic branch switching by calling `continuation` with the following options:

```julia
continuation(F, dF, d2F, d3F, br::ContResult, ind_bif::Int, optionsCont::ContinuationPar;
	kwargs...)
```

where `br` is a branch computed after a call to `continuation` with detection of bifurcation points enabled. This call computes the branch bifurcating from the `ind_bif `th bifurcation point in `br`. An example of use is provided in [2d generalized Bratu–Gelfand problem](@ref).

> See [Branch switching (branch point)](@ref) precise method definition

### Simple example

```@example TUT1
using BifurcationKit, Setfield, Plots

# vector field of transcritical bifurcation
F(x, p) = [x[1] * (p.μ - x[1])]

# vector of differentials (automatic differentiation)
jet = BifurcationKit.getJet(F; matrixfree = false)

# parameters of the vector field
par = (μ = -0.2, )

# compute branch of trivial equilibria (=0) and detect a bifurcation point
opts_br = ContinuationPar(dsmax = 0.05, ds = 0.01, detectBifurcation = 3, nev = 2)
br, = continuation(jet[1], jet[2], [0.1], par, (@lens _.μ), opts_br;
	recordFromSolution = (x, p) -> x[1])
	
# perform branch switching on one side of the bifurcation point
br1Top, = continuation(jet..., br, 1, setproperties(opts_br; maxSteps = 14); recordFromSolution = (x, p) -> x[1])

# on the other side
br1Bottom, = continuation(jet..., br, 1, setproperties(opts_br; ds = -opts_br.ds, maxSteps = 14); recordFromSolution = (x, p) -> x[1])

scene = plot(br, br1Top, br1Bottom; branchlabel = ["br", "br1Top", "br1Bottom"], legend = :topleft)
```


## Branch switching from non simple branch point to equilibria

We provide an automatic branch switching method in this case. The method is to first compute the reduced equation (see [Non-simple branch point](@ref)) and use it to compute the nearby solutions. These solutions are seeded as initial guess for [`continuation`](@ref). Hence, you can perform automatic branch switching by calling `continuation` with the following options:

```julia
continuation(F, dF, d2F, d3F, br::ContResult, ind_bif::Int, optionsCont::ContinuationPar;
	kwargs...)
```

An example of use is provided in [2d generalized Bratu–Gelfand problem](@ref).

> See [Branch switching (branch point)](@ref) for the precise method definition

## Branch switching from Hopf point to periodic orbits

In order to compute the bifurcated branch of periodic solutions at a Hopf bifurcation point, you need to choose a method to compute periodic orbits among:

- [Periodic orbits based on Trapezoidal rule](@ref)
- [Periodic orbits based on orthogonal collocation](@ref)
- [Periodic orbits based on the shooting method](@ref)

Once you have decided which method to use, you can call the following method:

```julia
continuation(F, dF, d2F, d3F, br::ContResult, ind_bif::Int, _contParams::ContinuationPar,
	prob::AbstractPeriodicOrbitProblem ;
	δp = nothing, ampfactor = 1, kwargs...)
```

We refer to [`continuation`](@ref) for more information about the arguments. Here, we just say a few words about how we can specify `prob::AbstractPeriodicOrbitProblem`.

- For [Periodic orbits based on Trapezoidal rule](@ref), you can pass `prob = PeriodicOrbitTrapProblem(M = 51)` where `M` is the number of times slices in the periodic orbit.

- For [Periodic orbits based on orthogonal collocation](@ref), you can pass `prob = PeriodicOrbitOCollProblem(M, m)` where `M` is the number of times slices in the periodic orbit and `m` is the degree of the collocation polynomials.

- For [Periodic orbits based on the shooting method](@ref), you need more parameters. For example, you can pass `prob = ShootingProblem(2, prob, Euler())` or `prob = PoincareShootingProblem(2, prob, Euler())` where `prob::ODEProblem` is an ODE problem to specify the Cauchy problem.

Several examples are provided like [1d Brusselator (automatic)](@ref) or [2d Ginzburg-Landau equation (finite differences, codim 2, Hopf aBS)](@ref).

> See [Branch switching (Hopf point)](@ref) for the precise method definition

!!! tip "Precise options"
    Although very convenient, the automatic branch switching does not allow the very fine tuning of parameters. It must be used as a first attempt before recurring to manual branch switching

## Branch switching from Branch / Period-doubling point of curve of periodic orbits

We do not provide (for now) the associated normal forms to these bifurcations of periodic orbits. As a consequence, the user is asked to provide the amplitude of the bifurcated solution.

We provide the branching method for the following methods to compute periodic orbits: [`PeriodicOrbitTrapProblem`](@ref),[`ShootingProblem`](@ref) and [`PoincareShootingProblem`](@ref). The call is as follows. Please note that a deflation is included in this method to simplify branch switching.

An example of use is provided in [Period doubling in Lur'e problem (PD aBS)](@ref).

```julia
continuation(br::AbstractBranchResult, ind_bif::Int, contParams::ContinuationPar;
	δp = 0.1, ampfactor = 1, usedeflation = false, kwargs...)
```

## Branch switching from Bogdanov-Takens (BT) point to Fold / Hopf curve

We provide an automatic branch switching method in this case (see for example [Extended Lorenz-84 model (codim 2 + BT/ZH aBS)](@ref) or [2d Ginzburg-Landau equation (finite differences, codim 2, Hopf aBS)](@ref)). Hence, you can perform automatic branch switching by calling `continuation` with the following options:

```julia
continuation(F, dF, d2F, d3F, br::ContResult, ind_bif::Int, options_cont::ContinuationPar = br.contparams; Jᵗ = nothing, δ::Real = 1e-8, δp = nothing, ampfactor::Real = 1,
			nev = options_cont.nev,
			issymmetric = false,
			detectCodim2Bifurcation::Int = 0,
			startWithEigen = false,
			autodiff = false,
			Teigvec = getvectortype(br),
			scaleζ = norm,
			kwargs...) where {Ta, Teigvals, Teigvecbr, Biftype, Ts, Tfunc <: AbstractProblemMinimallyAugmented, Tpar, Tl <: Lens}
```

where `ind_bif` is the index of the BT point in `br`. Note that the BT has been detected during Fold or Hopf continuation. Calling the above method thus switches from Fold continuation to Hopf continuation (and vice-versa) automatically with the same parameter axis.

!!! warning "WIP"
    This is still work in progress. Please report any error.
