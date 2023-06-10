# Overview of capabilities

```@contents
Pages = ["capabilities.md"]
Depth = 3
```

## Main features

- Newton-Krylov solver with generic linear / eigen *preconditioned* solver. Idem for the arc-length continuation.
- Newton-Krylov solver with nonlinear deflation and preconditioner. It can be used for branch switching for example. It is used for deflated continuation.
- Continuation written as an [iterator](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/iterator/)
- Monitoring user functions along curves computed by continuation, see [events](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/EventCallback/)
- Continuation methods: PALC, Moore Penrose, Multiple, Polynomial, Deflated continuation, ANM, ...
- Bifurcation points located with bisection

## Capabilities related to equilibria
- Detection of Branch, Fold, Hopf bifurcation points of stationary solutions and computation of their normal form.
- Automatic branch switching at branch points (**whatever the dimension of the kernel**) to equilibria
- **Automatic computation of bifurcation diagrams of equilibria**
- Fold / Hopf continuation based on Minimally Augmented formulation, with Matrix Free / Sparse / Dense Jacobian.
- Detection of all codim 2 bifurcations of equilibria and computation of the normal forms of Bogdanov-Takens, Bautin and Cusp
- Branching from Bogdanov-Takens / Zero-Hopf / Hopf-Hopf points to Fold / Hopf curve

**Note that you can combine most solvers, like use Deflation for Periodic orbit computation or Fold of periodic orbits family.**

> Custom state means, we can use something else than `AbstractArray`, for example your own `struct`.

|Features|Matrix Free|Custom state| [Tutorial](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/tutorials/) | GPU |
|---|---|---|---|---|
| (Deflated) Krylov-Newton| Yes| Yes| [link](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/deflatedproblem/)| Yes|
| Continuation PALC (Natural, Secant, Tangent, Polynomial) | Yes| Yes |All  | Yes |
| Deflated Continuation | Yes| Yes| [link](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/tutorialCarrier/#Deflated-Continuation-in-the-Carrier-Problem-1) |Yes  |
| Bifurcation / Fold / Hopf point detection | Yes| Yes| All / All / [link](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/tutorials/#Bifurcation-diagrams-with-periodic-orbits-1) | Yes |
| Fold Point continuation | Yes| Yes| [link](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/tutorials1/#Temperature-model-(simplest-example-for-equilibria)-1), [link](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/tutorialsCGL/#Complex-Ginzburg-Landau-2d-1), [ODE](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/lorenz84/#Extended-Lorenz-84-model-(codim-2-BT/ZH-aBS)) | Yes |
| Hopf Point continuation | Yes| `AbstractArray` | [ODE](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/lorenz84/#Extended-Lorenz-84-model-(codim-2-BT/ZH-aBS)) ||
| Branch point / Fold / Hopf normal form | Yes| Yes|  | Yes | [link](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/tutorials3/#Continuation-of-Hopf-points-1) | |
| Branch switching at Branch points | Yes| `AbstractArray` | [link](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/tutorials/#Bifurcation-diagrams-with-periodic-orbits-1) | Yes |
| **Automatic bifurcation diagram computation of equilibria** | Yes| `AbstractArray` |  [link](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/tutorialPP2/#pp2-example-from-AUTO07p-(aBD-Hopf-aBS)) | |
| Bogdanov-Takens / Bautin / Cusp / Zero-Hopf / Hopf-Hopf point detection | Yes| Yes | [ODE](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/lorenz84/#Extended-Lorenz-84-model-(codim-2-BT/ZH-aBS)) | |
| Bogdanov-Takens / Bautin / Cusp normal forms | Yes| `AbstractArray`| [ODE](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/lorenz84/#Extended-Lorenz-84-model-(codim-2-BT/ZH-aBS))| Yes |
| Branching from Bogdanov-Takens / Zero-Hopf / Hopf-Hopf to Fold / Hopf curve | Yes | `AbstractArray` | [ODE](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/lorenz84/#Extended-Lorenz-84-model-(codim-2-BT/ZH-aBS))|  |


## Capabilities related to Periodic orbits

- Periodic orbit computation and continuation using **parallel** (Standard or Poincaré) Shooting, Finite Differences or Orthogonal Collocation.
- Automatic branch switching at simple Hopf points to periodic orbits
- Detection of Branch, Fold, Neimark-Sacker, Period Doubling bifurcation points of periodic orbits.
- Continuation of Fold of periodic orbits

> Legend for the table: Standard shooting (SS), Poincaré shooting (PS), Orthogonal collocation (OC), trapezoid (T).

|Features|Method|Matrix Free|Custom state| [Tutorial](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/tutorials/) | GPU |
|---|---|---|---|---|---|
| Branch switching at Hopf points |SS/PS/OC/T| See each|  | [link](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/tutorialsODE-PD/#Period-doubling-in-Lur'e-problem-(PD-aBS)) | |
| Newton / continuation | T | Yes| `AbstractVector` | [link](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/tutorials3/#Brusselator-1d-(automatic)-1), [link](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/tutorialsCGL/#Complex-Ginzburg-Landau-2d-1) | Yes|
| Newton / continuation |OC| | `AbstractVector` | [link](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/tutorialsODE/#Neural-mass-equation-(Hopf-aBS)) | |
| Newton / continuation |SS| Yes| `AbstractArray` |  [link](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/tutorialsODE-PD/#Period-doubling-in-Lur'e-problem-(PD-aBS)) | Yes|
| Newton / continuation |PS| Yes| `AbstractArray` |  [link](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/tutorials3/#d-Brusselator-(automatic)) | Yes|
| Fold, Neimark-Sacker, Period doubling detection |SS/PS/OC/T| See each| `AbstractVector` | [link](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/tutorials3/#d-Brusselator-(automatic))  | |
| Branch switching at Branch point |SS/PS/OC/T| See each|  | [link](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/tutorialsODE-PD/#Period-doubling-in-Lur'e-problem-(PD-aBS)) | |
| Branch switching at PD point |SS/PS/OC/T| See each|  | [link](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/tutorialsODE-PD/#Period-doubling-in-Lur'e-problem-(PD-aBS)) | |
| Continuation of Fold points |SS/PS/OC/T| See each| `AbstractVector` | [link](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/tutorialsCGL/#Continuation-of-Fold-of-periodic-orbits) | Yes |
| Continuation of Period-doubling points |SS/OC| | `AbstractVector` |  [link](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/tutorialsCodim2PO/#Continuation-of-Fold/PD-of-periodic-orbits-with-Shooting) | |
| Continuation of Neimark-Sacker points |SS/OC| | `AbstractVector` | [link](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/steinmetz/#Curve-of-NS-points-of-periodic-orbits) | |
| detection of codim 2 bifurcations of periodic orbits |SS/OC| | `AbstractVector` | [link](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/steinmetz/#Curve-of-NS-points-of-periodic-orbits) | |
| Branch switching at Bautin point to curve of Fold of periodic orbits |SS/OC| | `AbstractVector` |  | |
| Branch switching at ZH/HH point to curve of NS of periodic orbits |SS/OC| | `AbstractVector` |  | |

## Capabilities related to Homoclinic orbits

This is available through the plugin [HclinicBifurcationKit.jl](https://github.com/bifurcationkit/HclinicBifurcationKit.jl). Please see the [specific docs](https://bifurcationkit.github.io/HclinicBifurcationKit.jl/dev/) for more information.

- compute Homoclinic to Hyperbolic Saddle Orbits (HomHS) using Orthogonal collocation or Standard shooting
- compute bifurcation of HomHS
- start HomHS from a direct simulation
- automatic branch switching to HomHS from Bogdanov-Takes bifurcation point