# BifurcationKit.jl

```@contents
Pages = ["index.md"]
Depth = 2
```

`BifurcationKit.jl` is a comprehensive Julia package for **automatic bifurcation analysis** of large-scale equations F(u, λ)=0 where λ is a real parameter. It leverages iterative methods, dense/sparse formulations, and specialized hardware (e.g., GPUs) for high-performance computations.

## Key features

- **Continuation algorithms**: PALC, deflated continuation, and more
- **Newton-Krylov methods**: Efficient predictor-corrector schemes
- **Flexible eigensolvers**: Matrix-free, dense, and sparse eigensolvers for stability analysis
- **Periodic orbit computation**: Shooting, finite differences and collocation methods
- **Hardware flexibility**: CPU and GPU support for massive computations
- **General vector interface**: Does not require `u` to be an `AbstractArray`

> Despite its focus on large-scale problems, the package handles low-dimensional problems with equal ease.

**Unique capabilities**: `BifurcationKit.jl` is one of the few software packages that provides shooting methods alongside finite difference and collocation methods for computing periodic orbits of dynamical systems.

**Design philosophy**: The package avoids requiring `u` to be a subtype of `AbstractArray`, enabling the use of `KrylovKit.jl` for linear solvers. Note that this flexibility is not yet available for all methods in the package.

## 📦 Installation

This package requires Julia >= v1.3.0 (due to methods added to abstract types, see [#31916](https://github.com/JuliaLang/julia/pull/31916)).

To install the stable version:

```julia
using Pkg
Pkg.add("BifurcationKit")
```

Or in the Julia REPL package mode:

```julia
] add BifurcationKit
```

To install the development version:

```julia
] add BifurcationKit#master
```


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

The tutorials prioritize **clarity and simplicity** over maximum performance (with exceptions like the [2d Ginzburg-Landau equation](@ref cgl) and [Langmuir–Blodgett model](@ref langmuir)). These examples can be optimized further for specific use cases. The intricacies of PDEs make efficient code highly problem-dependent, and users should leverage the specific characteristics of their problems.

!!! note "Large scale computations"
    This package has been used for massive computations (millions of unknowns), on the GPU. You can look at the citations below, or for example at [article1](
    https://doi.org/10.48550/arXiv.2508.19134), [article2](https://comptes-rendus.academie-sciences.fr/mathematique/item/CRMATH_2022__360_G1_59_0/)

!!! danger "Performance of tutorials"
    Do not get fooled by the small dimensions of the discretization of the PDE in the tutorials. The documentation has to be run on github CI so the machines are not powerful. You can increase the discretization on your computer. 

## Required methods for custom arrays

If you use standard arrays, you can skip this section.

We make the same requirements as `KrylovKit.jl`. Hence, we refer to its [docs](https://jutho.github.io/KrylovKit.jl/stable/#Package-features-and-alternatives-1) and the package [VectorInterface.jl](https://github.com/Jutho/VectorInterface.jl?tab=readme-ov-file) for more information. We additionally require the following methods to be available:

- `Base.length(x)`: it is used in the constraint equation of the pseudo arclength continuation method (see [`continuation`](@ref) for more details). If `length` is not available for your "vector" type, define `length(x) = 1` and adjust the parameter `θ` in `PALC`.
- `Base.copyto!(dest, in)` this is used to reduce the allocations.
- `real` this is used to compute normal forms.
- `conj` this is used to compute normal forms.

## Citations

Papers citing this work are collected on [Zotero](https://www.zotero.org/groups/6097154/citations_of_bifurcationkit/library).

These citations are aggregated from [Google Scholar (search)](https://scholar.google.com/scholar?q=bifurcationkit&hl=en&as_sdt=0,5) and [Google Scholar (citations)](https://scholar.google.com/scholar?oi=bibs&hl=en&cites=159498619004863176,12573642401780006854,8662907770106865595). Note that each link may reference different subsets of papers.


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
