# [Branch switching](@id Branch-switching-page)

```@contents
Pages = ["intro-abs.md"]
Depth = 3
```

The precise definition of the methods is given in [Branch switching (branch point)](@ref) and [Branch switching (Hopf point)](@ref).

## Summary of branching procedures

We collect in the following table the list of automatic branch switching (aBS) functions. Their detailed explanation follows in this page.

| function | ind-th bif. point | Type `T` | description |
|---|---|---|---|
|  `continuation(br::ContResult{T}, ind::Int; kw...)` | `:bp`, `:nd`| `EquilibriumCont`  |  aBS from equilibria to equilibria  |
|  `continuation(br::ContResult{T}, ind::Int, lens2::Lens; kw...)` | `:bp`, `:hopf`| `EquilibriumCont` | Fold/Hopf continuation w.r.t. parameters `getlens(br)` and `lens2`  |
|  `continuation(br::ContResult{T}, ind::Int; kw...)` | `:bt,:zh,:hh`| ` FoldCont,HopfCont` | switch to Fold/Hopf continuation from Hopf/Fold w.r.t. parameters of codim 2 `br`  |
| `continuation(br::ContResult{T}, ind_hopf::Int, ::ContinuationPar, prob::AbstractPeriodicOrbitProblem)`   | `:hopf` |  `EquilibriumCont` | Branch switching from Hopf point to periodic orbits |
| `continuation(br::ContResult{T}, ind::Int, kw...)`   | `:bp,:pd` |  `PeriodicOrbitCont` | Branch switching from Branch / Period-doubling point of periodic orbits to curve of periodic orbits |
| `continuation(br::ContResult{T}, ind::Int, kw...)`   | `:gh,:zh,:hh` |  `TwoParamCont` | Branch switching from Bautin / Zero-Hopf/ Hopf-Hopf point to curve of Fold/NS of periodic orbits |
