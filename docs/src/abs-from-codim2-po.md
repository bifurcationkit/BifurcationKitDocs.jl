# From codim 2 to periodic orbits

```@contents
Pages = ["abs-from-codim2-po.md"]
Depth = 3
```

## From Bautin point to curve Folds of periodic orbits

From the [Bautin normal form](http://scholarpedia.org/article/Bautin_bifurcation), we know that there is a curve of Fold of periodic orbits near the bifurcation point.

We provide an automatic branch switching method in this case which reads as follows:

```julia
continuation(br::HopfCont, ind_BAUTIN::Int, 
	_contParams::ContinuationPar,
    prob::AbstractPeriodicOrbitProblem ;
    δp = nothing, ampfactor = 1, kwargs...)
``` 

where `prob` is a method to compute periodic orbits (see [From Hopf point to periodic orbits](@ref) for more information).

Note that the two parameters in `br` will be used for the continuation of Fold points of periodic orbits.

See [ODE](@id lorenz98-take2) for an example of use.

## From Zero-Hopf (ZH) point to curve NS of periodic orbits

From the [Zero-Hopf normal form](http://scholarpedia.org/article/Zero-Hopf_bifurcation), we know that there is a curve of Neimark-Sacker (NS) bifurcations of periodic orbits near the bifurcation point.

We provide an automatic branch switching method in this case which reads as follows:

```julia
continuation(br::TwoParamCont, ind_ZH::Int, 
	_contParams::ContinuationPar,
    prob::AbstractPeriodicOrbitProblem ;
    δp = nothing, ampfactor = 1, kwargs...)
``` 

where `prob` is a method to compute periodic orbits (see [From Hopf point to periodic orbits](@ref) for more information).

Note that the two parameters in `br` will be used for the continuation of NS points of periodic orbits.

## From Hopf-Hopf (HH) point to curve NS of periodic orbits

From the [Hopf-Hopf normal form](http://scholarpedia.org/article/Hopf-Hopf_bifurcation), we know that there are two curves of Neimark-Sacker (NS) bifurcations of periodic orbits near the bifurcation point.

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

See [ODE](@id lorenz98-take2) for an example of use.
	