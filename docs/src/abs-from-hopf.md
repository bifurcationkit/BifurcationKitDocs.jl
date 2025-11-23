# From Hopf / PD / Branch point to periodic orbits

```@contents
Pages = ["abs-from-hopf.md"]
Depth = 3
```

## From Hopf point to periodic orbits
In order to compute the bifurcated branch of periodic solutions at a Hopf bifurcation point, you need to choose a method to compute periodic orbits among:

- [Periodic orbits based on Trapezoidal rule](@ref)
- [Periodic orbits based on orthogonal collocation](@ref)
- [Periodic orbits based on the shooting method](@ref)

Once you have decided which method to use, you can call the following:

```julia
continuation(br::ContResult, ind_HOPF::Int, _contParams::ContinuationPar,
	prob::AbstractPeriodicOrbitProblem ;
	δp = nothing, ampfactor = 1, kwargs...)
```

We refer to [`continuation`](@ref) for more information about the arguments. Here, we just say a few words about how we can specify `prob::AbstractPeriodicOrbitProblem`.

- For [Periodic orbits based on Trapezoidal rule](@ref), you can pass `PeriodicOrbitTrapProblem(M = 51)` where `M` is the number of times slices in the periodic orbit.

- For [Periodic orbits based on orthogonal collocation](@ref), you can pass `PeriodicOrbitOCollProblem(M, m)` where `M` is the number of times slices in the periodic orbit and `m` is the degree of the collocation polynomials.

- For [Periodic orbits based on the shooting method](@ref), you need more parameters. For example, you can pass `ShootingProblem(M, odeprob, Euler())` or `PoincareShootingProblem(M, odeprob, Euler())` where `odeprob::ODEProblem` (see [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/types/ode_types/)) is an ODE problem to specify the Cauchy problem amd `M` is the number of sections.

> See [Branch switching (Hopf point)](@ref) for the precise method definition

### Algorithm

The algorithm proceeds as follows. The normal form of the Hopf bifurcation is first computed. Then a predictor for the bifurcated branch of periodic orbits is generated from the normal form. Finally, this predictor is used as a guess for the computation of periodic orbits.

### Example

The simplest example is from the [getting-started section](@ref gt-hopf) which we repeat partially below.
Several examples are provided in [example ODE](@ref nmepo). In the case of PDE, you can have a look at [Brusselator](@ref brusauto) or [2d Ginzburg-Landau equation](@ref cgl).

We compute a branch with a Hopf bifurcation:

```@example hopf_abs
using BifurcationKit, Plots

function Fsl(X, p)
    (;r, μ, ν, c3) = p
    u, v = X
    ua = u^2 + v^2
    [
        r * u - ν * v - ua * (c3 * u - μ * v)
        r * v + ν * u - ua * (c3 * v + μ * u)
    ]
end

par_sl = (r = 0.1, μ = 0., ν = 1.0, c3 = 1.0)
u0 = zeros(2)
prob = BifurcationProblem(Fsl, u0, par_sl, (@optic _.r))
opts = ContinuationPar()
br = continuation(prob, PALC(), opts, bothside = true)
```

We then compute the branch of periodic solutions using orthogonal collocation (for example):

```@example hopf_abs
br_po = continuation(br, 2, opts,
        PeriodicOrbitOCollProblem(20, 5)
        )
plot(br, br_po)
```
## From Period-doubling point to curve of periodic orbits

For all cases, an example of use is provided in [Lur'e problem](@ref pdlure).

### Case of Shooting and Collocation

We provide an automatic branching procedure in this case. In essence, all you have to do is to call

```julia
continuation(br::ContResult, ind_PD::Int, _contParams::ContinuationPar;
    prm = Val(true), detailed = Val(true),
    kwargs...)
```

The option `prm = Val(true)` enforces that the period-doubling normal form is computed using the Poincaré return map ; this is only necessary in case of use of the collocation method. Indeed, in the case of the collocation method, an automatic procedure based on the Iooss normal form has yet to be implemented.

### Case of Trapezoid method

We do not provide (for now) the automatic branching procedure for these bifurcations of periodic orbits. As a consequence, the user is asked to provide the amplitude of the bifurcated solution.

The call is as follows. Please note that a deflation is included in this method to simplify branch switching.

```julia
continuation(br::AbstractBranchResult, ind_PD::Int, contParams::ContinuationPar;
	δp = 0.1, ampfactor = 1, 
	usedeflation = false,
	kwargs...)
```

### Algorithm

The algorithm proceeds as follows. The normal form of the Period-doubling bifurcation is first computed. Then a predictor for the bifurcated branch of periodic orbits is generated from the normal form. Finally, this predictor is used as a guess for the computation of periodic orbits.

## From Branch point to curve of periodic orbits

We do not provide (for now) the automatic branching procedure for these bifurcations of periodic orbits. As a consequence, the user is asked to provide the amplitude of the bifurcated solution.

We provide the branching method for the following methods to compute periodic orbits: [`PeriodicOrbitTrapProblem`](@ref), [`PeriodicOrbitOCollProblem`](@ref), [`ShootingProblem`](@ref) and [`PoincareShootingProblem`](@ref). The call is as follows. Please note that a deflation is included in this method to simplify branch switching.

An example of use is provided in [Lur'e problem](@ref pdlure).

```julia
continuation(br::AbstractBranchResult, ind_PD::Int, contParams::ContinuationPar;
	δp = 0.1, ampfactor = 1, usedeflation = false, kwargs...)
```
