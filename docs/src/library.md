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
newton
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
PALC
```

```@docs
Polynomial
```

```@docs
MoorePenrose
```

```@docs
Multiple
```

```@docs
BifurcationKit.DefCont
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
continuation(br::ContResult, ind_bif::Int, optionsCont::ContinuationPar ; kwargs...)
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
getBranch
```

```@docs
getBranchesFromBP
```

```@docs
BifurcationKit.SpecialPoint
```

## Utils for periodic orbits

```@docs
getPeriod
```

```@docs
getAmplitude
```

```@docs
getMaximum
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
guessFromHopf(br, ind_hopf, eigsolver::AbstractEigenSolver, M, amplitude; phase = 0)
```

```@docs
getNormalForm
```