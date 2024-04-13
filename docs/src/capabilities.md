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
- Bifurcation points / events located with bisection
- Compatible with GPU

## Capabilities related to equilibria
- Detection of Branch, Fold, Hopf bifurcation points of stationary solutions and computation of their normal form.
- Automatic branch switching at branch points (**whatever the dimension of the kernel**) to equilibria
- **Automatic computation of bifurcation diagrams of equilibria**
- Fold / Hopf continuation based on Minimally Augmented formulation, with Matrix Free / Sparse / Dense Jacobian.
- Detection of all codim 2 bifurcations of equilibria and computation of the normal forms of Bogdanov-Takens, Bautin, Cusp, Zero-Hopf. (Hopf-Hopf normal form not implemented)
- Branching from Bogdanov-Takens / Zero-Hopf / Hopf-Hopf points to Fold / Hopf curve

## (limited) Capabilities related to maps
- continuation of fixed points of maps
- computation of normal form of Period-doubling, Neimark-Sacker, Branch point bifurcations.

**Note that you can combine most solvers, like use Deflation for Periodic orbit computation or Fold of periodic orbits family.**

> Custom state means, you can use something else than `AbstractArray`, for example your own `struct`.

|Features|Matrix Free|Custom state| [Tutorial](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/tutorials/) | GPU |
|---|---|---|---|---|
| (Deflated) Krylov-Newton| Yes| Yes| [link](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/deflatedproblem/)| Yes|
| Continuation PALC (Natural, Secant, Tangent, Polynomial) | Yes| Yes |All  | Yes |
| Deflated Continuation | Yes| Yes| [link](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/tutorialCarrier/#Deflated-Continuation-in-the-Carrier-Problem-1) |Yes  |
| Bifurcation / Fold / Hopf point detection | Yes| Yes| All / All / [link](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/tutorials/#Bifurcation-diagrams-with-periodic-orbits-1) | Yes |
| Fold Point continuation | Yes| Yes| [PDE](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/tutorials1/#Temperature-model-(simplest-example-for-equilibria)-1), [PDE](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/tutorialsCGL/#Complex-Ginzburg-Landau-2d-1), [ODE](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/lorenz84/#Extended-Lorenz-84-model-(codim-2-BT/ZH-aBS)) | Yes |
| Hopf Point continuation | Yes| `AbstractArray` | [ODE](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/lorenz84/#Extended-Lorenz-84-model-(codim-2-BT/ZH-aBS)) ||
| Branch point / Fold / Hopf normal form | Yes| Yes|  | Yes | [link](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/tutorials3/#Continuation-of-Hopf-points-1) | |
| Branch switching at Branch points | Yes| `AbstractArray` | [link](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/branchswitching/#From-simple-branch-point-to-equilibria) | Yes |
| **Automatic bifurcation diagram computation of equilibria** | Yes| `AbstractArray` |  [link](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/tutorialPP2/#pp2-example-from-AUTO07p-(aBD-Hopf-aBS)) | |
| Bogdanov-Takens / Bautin / Cusp / Zero-Hopf / Hopf-Hopf point detection | Yes| Yes | [ODE](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/lorenz84/#Extended-Lorenz-84-model-(codim-2-BT/ZH-aBS)) | |
| Bogdanov-Takens / Bautin / Cusp normal forms | Yes| `AbstractArray`| [ODE](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/lorenz84/#Extended-Lorenz-84-model-(codim-2-BT/ZH-aBS))| Yes |
| Branching from Bogdanov-Takens / Zero-Hopf / Hopf-Hopf to Fold / Hopf curve | Yes | `AbstractArray` | [ODE](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/lorenz84/#Extended-Lorenz-84-model-(codim-2-BT/ZH-aBS))|  |


## Capabilities related to Periodic orbits (PO)

- PO computation and continuation using **parallel** (Standard or Poincaré) Shooting, Finite Differences or Orthogonal Collocation (mesh adaptive).
- Automatic branch switching from simple Hopf points to PO
- Automatic branch switching from simple Period-Doubling points to PO
- Assisted branch switching from simple Branch points to PO
- Detection of Branch, Fold, Neimark-Sacker (NS), Period Doubling (PD) bifurcation points of PO.
- Fold / PD / NS continuation based on Minimally Augmented formulation (for shooting and collocation). Trapezoid method only allows continuing Fold of PO.
- Detection of all codim 2 bifurcations of PO (R1, R2, R3, R4, GPD, NS-NS, Chenciner, Fold-Flip, Fold-NS, PD-NS)
- Computation of the normal forms of PD, NS (for shooting and collocation) using the method based on Poincaré return map or the Iooss normal form (more precise).
- automatic branching from Bautin to curve of Fold of PO
- automatic branching from Zero-Hopf to curve of NS of PO
- automatic branching from Hopf-Hopf to curve of NS of PO

> Legend for the table: Standard shooting (SS), Poincaré shooting (PS), Orthogonal collocation (OC), trapezoid (T).

|Features|Method|Matrix Free|Custom state| [Tutorial](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/tutorials/) | GPU |
|---|---|---|---|---|---|
| Branch switching at Hopf points |SS/PS/OC/T| See each|  | [ODE](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/tutorialsODE-PD/#Period-doubling-in-Lur'e-problem-(PD-aBS)) | |
| Newton / continuation | T | Yes| `AbstractVector` | [PDE](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/tutorials3/#Brusselator-1d-(automatic)-1), [PDE](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/tutorialsCGL/#Complex-Ginzburg-Landau-2d-1) | Yes|
| Newton / continuation |OC| | `AbstractVector` | [ODE](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/tutorialsODE/#Neural-mass-equation-(Hopf-aBS)) | |
| Newton / continuation |SS| Yes| `AbstractArray` |  [ODE](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/tutorialsODE-PD/#Period-doubling-in-Lur'e-problem-(PD-aBS)) | Yes|
| Newton / continuation |PS| Yes| `AbstractArray` |  [PDE](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/tutorials3/#d-Brusselator-(automatic)) | Yes|
| Fold, Neimark-Sacker, Period doubling detection |SS/PS/OC/T| See each| `AbstractVector` | [link](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/tutorials3/#d-Brusselator-(automatic))  | |
| Branch switching at Branch point |SS/PS/OC/T| See each|  | [ODE](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/tutorialsODE-PD/#Period-doubling-in-Lur'e-problem-(PD-aBS)) | |
| Branch switching at PD point |SS/PS/OC/T| See each|  | [ODE](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/tutorialsODE-PD/#Period-doubling-in-Lur'e-problem-(PD-aBS)) | |
| Continuation of Fold points |SS/PS/OC/T| See each| `AbstractVector` |[ODE](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/tutorialsCodim2PO/#Periodic-predator-prey-model) [PDE](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/tutorialsCGL/#Continuation-of-Fold-of-periodic-orbits) | Yes |
| Continuation of Period-doubling points |SS/OC| | `AbstractVector` |  [ODE](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/tutorialsCodim2PO/#Continuation-of-Fold/PD-of-periodic-orbits-with-Shooting) | |
| Continuation of Neimark-Sacker points |SS/OC| | `AbstractVector` | [ODE](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/steinmetz/#Curve-of-NS-points-of-periodic-orbits) | |
| detection of codim 2 bifurcations of periodic orbits |SS/OC| | `AbstractVector` | [ODE](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/steinmetz/#Curve-of-NS-points-of-periodic-orbits) | |
| Branch switching at Bautin point to curve of Fold of periodic orbits |SS/OC| | `AbstractVector` | [ODE](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/lorenz84-PO/#Lorenz-84-model,-take-2.) | |
| Branch switching at ZH/HH point to curve of NS of periodic orbits |SS/OC| | `AbstractVector` | [ODE](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/lorenz84-PO/#Lorenz-84-model,-take-2.) | |

## Capabilities related to Homoclinic orbits

This is available through the plugin [HclinicBifurcationKit.jl](https://github.com/bifurcationkit/HclinicBifurcationKit.jl). Please see the [specific docs](https://bifurcationkit.github.io/HclinicBifurcationKit.jl/dev/) for more information.

- compute Homoclinic to Hyperbolic Saddle Orbits (HomHS) using Orthogonal collocation or Standard shooting
- compute bifurcation of HomHS
- start HomHS from a direct simulation
- automatic branch switching to HomHS from Bogdanov-Takes bifurcation point

## List of detected bifurcations

A **left-to-right** arrow in the following graph from $E_1$ to $E_2$ means that $E_2$ can be detected when continuing an object of type $E_1$.

A **right-to-left** arrow from $E_2$ to $E_1$ means that we can start the computation of object of type $E_1$ from $E_2$.

Each object of codim 0 (resp. 1) can be continued with 1 (resp. 2) parameters.

```mermaid
graph LR
    S[  ]
    C[ Equilibrium ]
    PO[ Periodic orbit ]
    BP[ Fold/simple branch point ]
    H[ Hopf \n :hopf]
    CP[Cusp]
    BT[ Bogdanov-Takens \n :bt ]
    ZH[Zero-Hopf \n :zh]
    GH[Bautin \n :gh]
    HH[Hopf-Hopf \n :hh]
    FPO[ Fold Periodic orbit ]
    NS[ Neimark-Sacker \n :ns]
    PD[ Period Doubling \n :pd ]
    BPC[BPC]
    CH[Chenciner \n :ch]
    GPD[Generalized period doubling \n :gpd]
    BPC[Branch point PO]
    LPPD[Fold-Flip]
    LPNS[Fold-NeimarkSacker]
    R1[1:1 resonance point\n :R1]
    R2[1:2 resonance point\n :R2]
    R3[1:3 resonance point\n :R3]    
    R4[1:4 resonance point\n :R4]

    S --> C
    S --> PO
    C --> nBP[ non simple\n branch point ]
    C --> BP
    C --> H

    BP --> CP
    BP <--> BT

    PO --> H
    PO --> FPO
    PO --> NS
    PO --> PD

    FPO <--> GH
    FPO <--> BPC
    FPO --> R1

    NS --> R1
    NS --> R3
    NS --> R4
    NS --> CH
    NS --> LPNS
    NS --> NSNS
    NS --> R2
    NS --> PDNS

    PD --> PDNS
    PD --> R2
    PD --> LPPD
    PD --> GPD

    H <--> BT
    H <--> ZH
    BP <--> ZH
    H <--> HH
    H <--> GH
    NS <--> ZH
    PO <--> BPC
    NS <--> HH
    FPO --> LPNS
    FPO --> LPPD

    _ --> Codim0 --> Codim1 --> Codim2
```