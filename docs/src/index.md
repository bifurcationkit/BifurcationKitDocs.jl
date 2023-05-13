# BifurcationKit.jl

This Julia package aims at performing **automatic bifurcation analysis** of possibly large dimensional equations F(u, λ)=0 where λ∈ℝ by taking advantage of iterative methods, dense / sparse formulation and specific hardwares (*e.g.* GPU).

It incorporates continuation algorithms (PALC, deflated continuation, ...) based on a Newton-Krylov method to correct the predictor step and a Matrix-Free/Dense/Sparse eigensolver is used to compute stability and bifurcation points.

> Despite initial focus on large scale problems, the package can easily handle low dimensional problems.

By leveraging on the above method, the package can also seek for periodic orbits of Cauchy problems by casting them into an equation $F(u,p)=0$ of high dimension. **It is by now, one of the only softwares which provides shooting methods AND methods based on finite differences or collocation to compute periodic orbits.**

The current package focuses on large scale nonlinear problems and multiple hardwares. Hence, the goal is to use Matrix Free methods on **GPU** (see [PDE example](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/tutorials2b/#d-Swift-Hohenberg-equation-(non-local)-on-the-GPU,-periodic-BC-(Advanced)) and [Periodic orbit example](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/tutorialsCGL/#Continuation-of-periodic-orbits-on-the-GPU-(Advanced))) or on a **cluster** to solve non linear PDE, nonlocal problems, compute sub-manifolds...

One design choice is that we try not to require `u` to be a subtype of an `AbstractArray` as this would forbid the use of spectral methods like the one from `ApproxFun.jl`. For now, our implementation does not allow this for Hopf continuation and computation of periodic orbits.

## Installation

This package requires Julia >= v1.3.0 because of the use of methods added to abstract types (see [#31916](https://github.com/JuliaLang/julia/pull/31916)).

To install it, please run

`] add BifurcationKit`

To install the bleeding edge version, please run

`] add BifurcationKit#master`

## Citing this work
If you use this package for your work, we ask that you **cite** the following paper!! Open source development strongly depends on this. It is referenced on HAL-Inria as follows:

```
@misc{veltz:hal-02902346,
  TITLE = {{BifurcationKit.jl}},
  AUTHOR = {Veltz, Romain},
  URL = {https://hal.archives-ouvertes.fr/hal-02902346},
  INSTITUTION = {{Inria Sophia-Antipolis}},
  YEAR = {2020},
  MONTH = Jul,
  KEYWORDS = {pseudo-arclength-continuation ; periodic-orbits ; floquet ; gpu ; bifurcation-diagram ; deflation ; newton-krylov},
  PDF = {https://hal.archives-ouvertes.fr/hal-02902346/file/354c9fb0d148262405609eed2cb7927818706f1f.tar.gz},
  HAL_ID = {hal-02902346},
  HAL_VERSION = {v1},
}
```

## Other softwares

There are many good softwares already available.

- For continuation in small dimension, most softwares are listed on [DSWeb](https://dsweb.siam.org/Software). One can mention the widely used AUTO-07p, or also, [XPPAUT](http://www.math.pitt.edu/~bard/xpp/xpp.html), [MATCONT](https://sourceforge.net/projects/matcont/), [PyDSTool](https://github.com/robclewley/pydstool) and [COCO](https://sourceforge.net/projects/cocotools/). All these are very reliable and some address high codimensional bifurcations.

- For large scale problems, there is the versatile and feature full [pde2path](http://www.staff.uni-oldenburg.de/hannes.uecker/pde2path/) but also [Trilinos-LOCA](https://trilinos.github.io/nox_and_loca.html), [CL_MATCONTL](https://github.com/careljonkhout/cl_matcontL), [COCO](https://sourceforge.net/projects/cocotools/), [GetFEM](https://getfem.org/userdoc/model_continuation.html) and the python libraries [pyNCT](https://pypi.org/project/PyNCT/) and [pacopy](https://github.com/nschloe/pacopy). There are also some nice [lectures](https://zenodo.org/record/3821169#.Y-zsAy8w08Q) by D. Avitabile *et al.* on matrix free secant continuation based on a concise Matlab implementation which are used in the influential paper [^Rankin]. At the time of initial release of `BifurcationKit.jl`, these [lectures](https://zenodo.org/record/3821169#.Y-zsAy8w08Q) provided one of the only libraries for matrix-free continuation, much easier to use than [Trilinos-LOCA](https://trilinos.github.io/nox_and_loca.html).
- For *deflated continuation*, there is [defcont](https://bitbucket.org/pefarrell/defcon/src/master/) (by the inventor of the algo. P. E. Farrell) and this [code](https://github.com/evstigneevnm/deflated_continuation) by N. M. Evstigneev.

In Julia, we also have [Bifurcations.jl](https://github.com/tkf/Bifurcations.jl) which is unmaintained.

## A word on performance

The examples which follow have not **all** been written with the goal of performance but rather simplicity (except maybe [2d Ginzburg-Landau equation (finite differences, codim 2, Hopf aBS)](@ref) and [1d Langmuir–Blodgett transfer model (advanced)](@ref)). One could surely turn them into more efficient codes. The intricacies of PDEs make the writing of efficient code highly problem dependent and one should take advantage of every particularity of the problem under study.

For example, in one of the simplest tutorials, [Temperature model (Simplest example)](@ref), one could use `BandedMatrices.jl` for the jacobian and an inplace modification when the jacobian is called ; using a composite type would be favored. Porting them to GPU could be another option.

## Requested methods for Custom State
Needless to say, if you use regular arrays, you don't need to worry about what follows.

We make the same requirements as `KrylovKit.jl`. Hence, we refer to its [docs](https://jutho.github.io/KrylovKit.jl/stable/#Package-features-and-alternatives-1) for more information. We additionally require the following methods to be available:

- `Base.length(x)`: it is used in the constraint equation of the pseudo arclength continuation method (see [`continuation`](@ref) for more details). If `length` is not available for your "vector", define it `length(x) = 1` and adjust tuning the parameter `theta` in `ContinuationPar`.
- `Base.copyto!(dest, in)` this is required to reduce the allocations by avoiding too many copies
- `Base.eltype` must be extended to your vector type. It is mainly used for branching.

## Citations
The papers citing this work are collected on [google scholar](https://scholar.google.fr/scholar?hl=fr&as_sdt=2005&cites=159498619004863176%2C8662907770106865595&scipsc=&as_ylo=&as_yhi=).

## References

[^Rankin]:> J. Rankin et al., "Continuation of localized coherent structures in nonlocal neural field equations", SIAM J. Scientific Computing 36, pp. B70–B93 (2014): https://epubs.siam.org/doi/10.1137/130918721

## Reproducibility
```@raw html
<details><summary>The documentation of BifurcationKit was built using these direct dependencies,</summary>
```
```@example
using Pkg # hide
Pkg.status() # hide
```
```@raw html
</details>
```
```@raw html
<details><summary>and using this machine and Julia version.</summary>
```
```@example
using InteractiveUtils # hide
versioninfo() # hide
```
```@raw html
</details>
```
```@raw html
<details><summary>A more complete overview of all dependencies and their versions is also provided.</summary>
```
```@example
using Pkg # hide
Pkg.status(;mode = PKGMODE_MANIFEST) # hide
```
```@raw html
</details>
```
```@raw html
You can also download the
<a href="
```
```@eval
using TOML
version = TOML.parse(read("../../Project.toml",String))["version"]
name = TOML.parse(read("../../Project.toml",String))["name"]
link = "https://github.com/bifurcationkit/"*name*".jl/tree/gh-pages/v"*version*"/assets/Manifest.toml"
```
```@raw html
">manifest</a> file and the
<a href="
```
```@eval
using TOML
version = TOML.parse(read("../../Project.toml",String))["version"]
name = TOML.parse(read("../../Project.toml",String))["name"]
link = "https://github.com/bifurcationkit/"*name*".jl/tree/gh-pages/v"*version*"/assets/Project.toml"
```
```@raw html
">project</a> file.
```
