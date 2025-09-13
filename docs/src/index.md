# BifurcationKit.jl

```@contents
Pages = ["index.md"]
Depth = 2
```

This Julia package aims at performing **automatic bifurcation analysis** of possibly large dimensional equations F(u, λ)=0 where λ is real by taking advantage of iterative methods, dense / sparse formulation and specific hardwares (*e.g.* GPU).

It incorporates continuation algorithms (PALC, deflated continuation, ...) based on a Newton-Krylov method to correct the predictor step and a Matrix-Free/Dense/Sparse eigensolver is used to compute stability and bifurcation points.

> Despite initial focus on large scale problems, the package can easily handle low dimensional problems.

The package can also seek for periodic orbits of Cauchy problems . **It is by now, one of the only softwares which provides shooting methods and methods based on finite differences or collocation to compute periodic orbits.**

Hence, it is possible to study large scale nonlinear problems on different hardwares. 

One design choice is that we try not to require `u` to be a subtype of an `AbstractArray` as this would forbid the use of spectral methods like the one from `ApproxFun.jl`. For now, our implementation does not allow this for all methods of the package.

## 📦 Installation

This package requires Julia >= v1.3.0 because of the use of methods added to abstract types (see [#31916](https://github.com/JuliaLang/julia/pull/31916)).

To install it, please run

`] add BifurcationKit`

To install the bleeding edge version, please run

`] add BifurcationKit#master`

## 📚 Citing this work
If you use this package for your work, we ask that you **cite** the following paper!! Open source development strongly depends on this. It is referenced on [HAL-Inria](https://hal.archives-ouvertes.fr/hal-02902346) with *bibtex* entry [CITATION.bib](https://github.com/bifurcationkit/BifurcationKit.jl/blob/master/CITATION.bib).

## 🧑‍💻 Other softwares

There are many good softwares already available.

- For continuation in small dimension, most softwares are listed on [DSWeb](https://dsweb.siam.org/Software). One can mention the widely used [AUTO-07p](https://github.com/auto-07p/auto-07p), or also, [XPPAUT](http://www.math.pitt.edu/~bard/xpp/xpp.html), [MATCONT](https://sourceforge.net/projects/matcont/), [PyDSTool](https://github.com/robclewley/pydstool) and [COCO](https://sourceforge.net/projects/cocotools/). All these are very reliable and some address high codimensional bifurcations.

- For large scale problems, there is the versatile and feature full [pde2path](http://www.staff.uni-oldenburg.de/hannes.uecker/pde2path/) but also [Trilinos-LOCA](https://trilinos.github.io/nox_and_loca.html), [CL_MATCONTL](https://github.com/careljonkhout/cl_matcontL), [COCO](https://sourceforge.net/projects/cocotools/), [GetFEM](https://getfem.org/userdoc/model_continuation.html) and the python libraries [pyNCT](https://pypi.org/project/PyNCT/) and [pacopy](https://github.com/nschloe/pacopy). There are also some nice [lectures](https://zenodo.org/record/3821169#.Y-zsAy8w08Q) by D. Avitabile *et al.* on matrix free secant continuation based on a concise Matlab implementation which are used in the influential paper [^Rankin]. At the time of initial release of `BifurcationKit.jl`, these [lectures](https://zenodo.org/record/3821169#.Y-zsAy8w08Q) provided one of the only libraries for matrix-free continuation, much easier to use than [Trilinos-LOCA](https://trilinos.github.io/nox_and_loca.html).
- For *deflated continuation*, there is the python package [defcont](https://bitbucket.org/pefarrell/defcon/src/master/) (by the inventor of the algo. P. E. Farrell) and this C++ [code](https://github.com/evstigneevnm/deflated_continuation) by N. M. Evstigneev.

In Julia, we have the following tools:

- [HomotopyContinuation.jl](https://www.juliahomotopycontinuation.org) which focuses on continuation of zeros of polynomials.
- [Bifurcations.jl](https://github.com/tkf/Bifurcations.jl) which is unmaintained.
- [NumericalContinuation.jl](https://github.com/dawbarton/NumericalContinuation.jl) by David Barton.

## 🚀 A word on performance

The tutorials have not **all** been written with the goal of performance but rather simplicity (except maybe [2d Ginzburg-Landau equation](@ref cgl) and [Langmuir–Blodgett model](@ref langmuir)). One could surely turn them into more efficient codes. The intricacies of PDEs make the writing of efficient code highly problem dependent and one should take advantage of every particularity of the problem under study.

## Requested methods for custom arrays
Needless to say, if you use "regular" arrays, you don't need to worry about what follows.

We make the same requirements as `KrylovKit.jl`. Hence, we refer to its [docs](https://jutho.github.io/KrylovKit.jl/stable/#Package-features-and-alternatives-1) and the [vector-like interface](https://github.com/Jutho/VectorInterface.jl?tab=readme-ov-file) for more information. We additionally require the following methods to be available:

- `Base.length(x)`: it is used in the constraint equation of the pseudo arclength continuation method (see [`continuation`](@ref) for more details). If `length` is not available for your "vector", define `length(x) = 1` and adjust the parameter `θ` in `PALC`.
- `Base.copyto!(dest, in)` this is required to reduce the allocations
- `Base.eltype` must be extended to your vector type.

## Citations
The papers citing this work are collected on [google scholar](https://scholar.google.com/scholar?q=bifurcationkit&hl=en&as_sdt=0,5) and [zotero](https://www.zotero.org/groups/6097154/citations_of_bifurcationkit/library).

See also this other [google scholar](https://scholar.google.com/scholar?oi=bibs&hl=en&cites=159498619004863176,12573642401780006854,8662907770106865595) but less papers are referenced.


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

## References

[^Rankin]:> J. Rankin et al., "Continuation of localized coherent structures in nonlocal neural field equations", SIAM J. Scientific Computing 36, pp. B70–B93 (2014): https://epubs.siam.org/doi/10.1137/130918721
