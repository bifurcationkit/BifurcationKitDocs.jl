# From codim 2 to equilibria


```@contents
Pages = ["abs-from-codim2-eq.md"]
Depth = 3
```


## From Bogdanov-Takens (BT) point to Fold / Hopf curve

We provide an automatic branch switching method in this case (see for example [Extended Lorenz-84 model](@ref lorenz) or [2d Ginzburg-Landau equation (finite differences, codim 2, Hopf aBS)](@ref cgl)). Hence, you can perform automatic branch switching by calling `continuation` with the following options:

```julia
continuation(br::ContResult, ind_BT::Int,
	options_cont::ContinuationPar = br.contparams;
	nev = options_cont.nev,
	detect_codim2_bifurcation::Int = 0,
	start_with_eigen = false,
	autodiff = false,
	Teigvec = getvectortype(br),
	scaleζ = norm,
	kwargs...)
```

where `ind_BT` is the index of the BT point in `br`. Note that the BT has been detected during Fold or Hopf continuation. Calling the above method thus switches from Fold continuation to Hopf continuation (and vice-versa) automatically with the same parameter axis.

> Check the docs of [Fold / Hopf Continuation](@ref) and particularly [Setting the jacobian](@ref jac-fold) for improving the speed of computation for large scale systems.

## From Zero-Hopf (ZH) point to Fold / Hopf curve

We provide an automatic branch switching method in this case (see for example [Extended Lorenz-84 model](@ref lorenz) or [2d Ginzburg-Landau](@ref cgl)). Hence, you can perform automatic branch switching by calling `continuation` with the following options:

```julia
continuation(br::ContResult, ind_ZH::Int,
	options_cont::ContinuationPar = br.contparams;
	nev = options_cont.nev,
	detect_codim2_bifurcation::Int = 0,
	start_with_eigen = false,
	autodiff = false,
	Teigvec = getvectortype(br),
	scaleζ = norm,
	kwargs...)
```

where `ind_ZH` is the index of the ZH point in `br`. Note that the ZH has been detected during Fold or Hopf continuation. Calling the above method thus switches from Fold continuation to Hopf continuation (and vice-versa) automatically with the same parameter axis.

> Check the docs of [Fold / Hopf Continuation](@ref) and particularly [Setting the jacobian](@ref jac-fold) for improving the speed of computation for large scale systems.

## From Hopf-Hopf (HH) point to Fold / Hopf curve

We provide an automatic branch switching method in this case (see for example [Extended Lorenz-84 model](@ref lorenz) or [2d Ginzburg-Landau equation](@ref cgl)). Hence, you can perform automatic branch switching by calling `continuation` with the following options:

```julia
continuation(br::ContResult, ind_HH::Int,
	options_cont::ContinuationPar = br.contparams;
	δp = nothing, ampfactor::Real = 1,
	nev = options_cont.nev,
	detect_codim2_bifurcation::Int = 0,
	start_with_eigen = false,
	autodiff = false,
	Teigvec = getvectortype(br),
	scaleζ = norm,
	kwargs...)
```

where `ind_HH` is the index of the HH point in `br`. Note that the HH has been detected during Hopf continuation. Calling the above method thus switches from Hopf continuation to another Hopf branch automatically with the same parameter axis.

> Check the docs of [Fold / Hopf Continuation](@ref) and particularly [Setting the jacobian](@ref jac-fold) for improving the speed of computation for large scale systems.
	