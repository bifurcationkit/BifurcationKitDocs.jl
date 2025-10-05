# Library

```@meta
CollapsedDocStrings = true
```

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

```@docs
BifurcationKit.cbMaxNorm
```

```@docs
BifurcationKit.cbMaxNormAndΔp
```

## Results


```@docs
NonLinearSolution
```

```@docs
ContResult
```

```@docs
BifurcationKit.Branch
```


```@docs
BifurcationKit.SpecialPoint
```

```@docs
BifurcationKit.BifDiagNode
```


## Problems

```@docs
BifFunction
```

```@docs
BifurcationProblem
```

```@docs
ODEBifProblem
```

```@docs
DAEBifProblem
```

```@docs
PDEBifProblem
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

## [Linear solvers](@id Library-LS)

```@docs
BifurcationKit.DefaultLS
```

```@docs
BifurcationKit.DefaultPILS
```

```@docs
BifurcationKit.GMRESIterativeSolvers
```

```@docs
BifurcationKit.GMRESKrylovKit
```

```@docs
BifurcationKit.KrylovLS
```

```@docs
BifurcationKit.KrylovLSInplace
```

## [Eigen solvers](@id Library-EIG)

```@docs
BifurcationKit.DefaultEig
```

```@docs
BifurcationKit.EigArpack
```

```@docs
BifurcationKit.EigKrylovKit
```

```@docs
BifurcationKit.EigArnoldiMethod
```

## [Floquet solvers](@id Library-FLOQUET)

```@docs
FloquetQaD
```

### Collocation

```@docs
BifurcationKit.FloquetGEV
```

```@docs
BifurcationKit.FloquetColl
```

## [Bordered linear solvers](@id Library-BLS)

```@docs
BifurcationKit.MatrixBLS
```

```@docs
BifurcationKit.BorderingBLS
```

```@docs
BifurcationKit.MatrixFreeBLS
```

```@docs
BifurcationKit.MatrixFreeBLSmap
```

```@docs
BifurcationKit.LSFromBLS
```
## Nonlinear solver

```@docs
BifurcationKit.solve
```

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
PALC
```

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

## Utils for periodic orbits

```@docs
getperiod
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
guess_from_hopf
```

```@docs
get_normal_form
```
