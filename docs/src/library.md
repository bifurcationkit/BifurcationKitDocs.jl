# Library

```@contents
Pages = ["library.md"]
Depth = 3
```

## Parameters

```@docs
NewtonPar
```

```@docs
ContinuationPar
```

## Results


```@docs
NonLinearSolution
```

```@docs
ContResult
```

## Problems

```@docs
BifFunction
```

```@docs
BifurcationProblem
```

```@docs
DeflationOperator
```

```@docs
DeflatedProblem
```

### Periodic orbits

```@docs
PeriodicOrbitTrapProblem
```

```@docs
PeriodicOrbitOCollProblem
```

```@docs
ShootingProblem
```

```@docs
PoincareShootingProblem
```

### Waves

```@docs
BifurcationKit.TWProblem
```

## Newton

```@docs
solve
```

## [Continuation](@id Library-Continuation)

```@docs
BifurcationKit.DotTheta
```

```@docs
continuation
```

## Continuation algorithms

```@docs
Natural
```

```@docs
Secant
```


```@docs
Bordered
```

```@docs
Polynomial
```

```@docs
Multiple
```

```@docs
BifurcationKit.PALC
```

```@docs
BifurcationKit.AutoSwitch
```

```@docs
MoorePenrose
```

```@docs
BifurcationKit.DefCont
```

```@docs
AsymptoticNumericalMethod.ANM
```

## Events

```@docs
BifurcationKit.DiscreteEvent
```

```@docs
BifurcationKit.ContinuousEvent
```

```@docs
BifurcationKit.SetOfEvents
```

```@docs
BifurcationKit.PairOfEvents
```

## Branch switching (branch point)

```@docs
continuation(br::BifurcationKit.AbstractResult{BifurcationKit.EquilibriumCont}, ind_bif::Int, optionsCont::ContinuationPar ; kwargs...)
```

## Branch switching (Hopf point)
```@docs
continuation(br::BifurcationKit.AbstractBranchResult, ind_bif::Int, _contParams::ContinuationPar, prob::BifurcationKit.AbstractPeriodicOrbitProblem ; kwargs...)
```

## Bifurcation diagram

```@docs
bifurcationdiagram
```

```@docs
bifurcationdiagram!
```

```@docs
get_branch
```

```@docs
get_branches_from_BP
```

```@docs
BifurcationKit.SpecialPoint
```

## Utils for periodic orbits

```@docs
getperiod
```

```@docs
getamplitude
```

```@docs
getmaximum
```

```@docs
SectionSS
```

```@docs
SectionPS
```

## Misc.

```@docs
PrecPartialSchurKrylovKit
```

```@docs
PrecPartialSchurArnoldiMethod
```

```@docs
Flow
```

```@docs
FloquetQaD
```

```@docs
guess_from_hopf(br, ind_hopf, eigsolver::AbstractEigenSolver, M, amplitude; phase = 0)
```

```@docs
get_normal_form
```
