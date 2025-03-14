# Migration from previous versions

```@contents
Pages = ["migration.md"]
Depth = 2
```

We only highlight changes that are potentially breaking for the user.

## Version 0.4.0

We rely on `Accessors.jl` instead of `Setfield.jl`. This basically amounts to changing `@lens` by `@optic` in your code.

## Version 0.3.4
- correct selection of default linear solver for MoorePenrose
- pass iterator for plotting
- `update_section_every_step` becomes a UInt
- add fields in `PeriodDoublingProblemMinimallyAugmented` and `NeimarkSackerProblemMinimallyAugmented` for holding Resonance test values
- add more abstract types `<: AbstractWrapperPOProblem`
- introduce function `get_lenses`
- introduce new struct `FinalisePO` to wrap finalizers for periodic orbits
- `AbstractProblemMinimallyAugmented` becomes parametric
- rewrite `get_bif_point_codim2`
- add callback `cbMaxNormAndΔp`
- `FloquetWrapper` becomes mutable

## Version 0.3.3

## Migration from v0.2.x to v0.3.x

A new version v0.3 has been tagged in which the function names, keyword arguments,... follow the Julia convention. There are a lot of breaking changes. For example, `callbackN` has been changed to `callback_newton`.

## Migration from v0.1.x to v0.2.x

New version of the package with modified interface. You are now required to define a `BifurcationProblem` to perform continuation or bifurcation analysis. You also need to pass your plot/record functions. 

The previous interface is available under the tag 0.1.12 which can be installed by doing

`] add BifurcationKit@0.1.12`

The new version provides many bugs fix though.
(Please note that the docs are up to date).

### Don't use AD yourself

There is nothing wrong with doing so but this is done in the constructor of `BifurcationPoblem`, so if `myJacAD` is the jacobian computed using `ForwardDiff`, the declaration

```
prob = BifurcationProblem(F, x, p, lens ; J = myJacAD) 
```

should be 

```
prob = BifurcationProblem(F, x, p, lens) 
```

> There is nothing wrong in passing your own jacobian though

### Error: no method matching iterate(::BifurcationKit.ContResult

This is because you use the old syntax 

```julia
br, = continuation(...)
```

instead of (no comma)

```julia
br = continuation(...)
```

### Arguments to `continuation`

`recordFromSolution` and `plotFromSolution` should be passed to `BifurcationProblem` instead of `continuation`.