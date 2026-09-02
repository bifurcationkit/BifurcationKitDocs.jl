# Migration from previous versions

```@contents
Pages = ["migration.md"]
Depth = 2
```

We only highlight changes that are potentially breaking for the user.

## Version 0.8.x

Breaking changes are detailed per minor release below.

### Version 0.8.4

- `get_adjoint_basis` is replaced by `_get_kernel_basis_1d_from_eigensolver` / `_get_kernel_basis_1d_from_bls` (resp. the Nd helpers), which now compute the right and left (adjoint) kernel vectors together

### Version 0.8.3

- `detect_codim2_parameters` is renamed `modify_contparams_for_codim2`
- the `_compute_bordered_vectors` method of the Fold/Hopf/PD/NS MA formulations is split into `__compute_bordered_vectors_fold` / `__compute_bordered_vectors_hopf`

### Version 0.8.2

- no breaking change

### Version 0.8.1

- `update!(prob, x)` becomes `restore_problem!(prob, x, pars)`
- `getlinsolver` / `getbls` become `get_bordered_linsolver`
- `OneParamCont`, `TwoParamCont`, `TwoParamPeriodicOrbitCont` become `AbstractOneParamCont`, `AbstractTwoParamCont`, `AbstractTwoParamPeriodicOrbitCont`
- `_getsolution` becomes `saved_solution`
- the `tangent` field in `MoorePenrose` becomes `predictor`
- the `autodiff` keyword argument for computing normal forms is removed
- `FlowDE` construction is now keyword-based
- the BT fields `nfsupp` are merged into `nf`

## Version 0.8.0

- new BVP interface: `BVPModel`, `PeriodicOrbitModel`, `DiscretizedBVP`, `discretize`, `generate_solution`, the problem `BVPBifProblem` and the discretizers `Shooting`, `Trapeze`, `Collocation` (with mesh adaptation)
- `SolPeriodicOrbit` becomes `BVPSolution`, `POSolution` becomes `POInterpolation` and `POSolutionAndState` becomes `POSavedSolutionAndState`
- `Trap` becomes `Trapeze` in the BVP interface
- `jacobian = :auto` becomes `jacobian = AutoDiffDense()` for BVP Shooting
- add `normal_form` alias for `get_normal_form`

## Version 0.7.x

- remove `AbstractPeriodicOrbitDiscretization`: the subtypes now inherit directly from `AbstractBoundaryValueDiscretization`
- `AbstractPODifferentialDiscretization` becomes `AbstractDifferentialDiscretization`
- `AbstractPOFiniteDifferencesDiscretization` becomes `AbstractFiniteDifferencesDiscretization`
- `AbstractPOShootingDiscretization` becomes `AbstractShootingDiscretization`
- `*ShootingProblem` becomes `*Shooting`; `PeriodicOrbitOCollProblem` becomes `Collocation`; `PeriodicOrbitTrapProblem` becomes `Trapeze`

## Version 0.6.x

- add `AbstractBoundaryValueDiscretization` and the discretization hierarchy `AbstractPeriodicOrbitDiscretization`, `AbstractPODifferentialDiscretization`, `AbstractPOFiniteDifferencesDiscretization`
- `AbstractPoincareShootingProblem` becomes `AbstractPoincareShootingDiscretization`
- `AbstractShootingProblem` becomes `AbstractPOShootingDiscretization`
- `BTProblemMinimallyAugmented` becomes `BTMinimallyAugmentedFormulation`
- `TWProblem` becomes `TWModel`
- `Fold/Hopf/PeriodDoubling/NeimarkSackerProblemMinimallyAugmented` become `Fold/Hopf/PeriodDoubling/NeimarkSackerMinimallyAugmentedFormulation`
- `AbstractProblemMinimallyAugmented` becomes `AbstractMinimallyAugmentedFormulation`
- `WrapPOColl` becomes `PeriodicOrbitFunctionalColl`, `WrapPOSh` becomes `PeriodicOrbitFunctionalSh`, `WrapPOTrap` becomes `PeriodicOrbitFunctionalTrap`
- `correct_bifurcation` becomes `_correct_event_labels`
- add `PeriodicOrbit{Tdisc}`, `TravellingWave{Tdisc}` to bridge discretization and problem
- add accessors `get_discretization`, `get_formulation`, `get_solution`, `getparams`
- change signature of `continuation` for periodic orbits: `prob::AbstractPeriodicOrbitProblem` becomes `disc::AbstractPeriodicOrbitDiscretization`
- remove the `params` and `lens` fields from `PeriodicOrbitFunctionalTrap`, `PeriodicOrbitFunctionalSh`, `PeriodicOrbitFunctionalColl`, `WrapTW` and the `param` field from the `*MAProblem` structs; `getparams` now delegates to `getparams(get_formulation(prob))`
- `DefaultLS` uses `VI.Zero()` / `VI.One()` instead of literal `0` / `1`

## Version 0.5.x

- `jacobian` parameter in `PeriodicOrbitTrapProblem` changes from `Symbol` to custom type
- `FloquetCollGEV` becomes `FloquetGEV` and `extract_period` becomes `_extract_period`
- export `ODEBifProblem` and `jad` becomes `jacobian_adjoint`
- `BorderedArray` complies with `VectorInterface.jl` and most of the code as well
- `multicontinuation` returns a `Branch` instead of `Vector{Branch}`

## Version 0.4.0

We rely on `Accessors.jl` instead of `Setfield.jl`. This basically amounts to changing `@lens` by `@optic` in your code.

- `jacobian_ma` argument changes from `Symbol` to type. For example, `continuation(br, 1; jacobian_ma = :minaug)` becomes `continuation(br, 1; jacobian_ma = MinAug())`.
- add `update!` function to `BifurcationProblem`, which allows to adapt the problem during continuation.

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
- `getvectortype` becomes `_getvectortype` and `hasstability` becomes `_hasstability`

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
