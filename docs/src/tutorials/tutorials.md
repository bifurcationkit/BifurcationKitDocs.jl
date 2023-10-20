# Tutorials

```@contents
Pages = ["tutorials.md"]
Depth = 2
```

There are three levels of tutorials:

1. fully **automatic bifurcation diagram** (**aBD**) computation (only for equilibria): one uses the function `bifurcationdiagram` and let it compute the diagram fully automatically. Another possibility is to use **deflated continuation**.
2. semi-automatic bifurcation diagram computation: one uses **automatic branch switching** (**aBS**) to compute branches at specified bifurcation points
3. manual bifurcation diagram computation: one does not use automatic branch switching. This has only educational purposes or for complex problems where aBS fails.

## ODE examples

These examples are specific to ODEs. 

### Computation of equilibria

```@contents
Pages = ["ode/tutorialsODE.md","ode/tutorialCO.md","ode/lorenz84.md", "ode/tutorialPP2.md",]
Depth = 1
```

### Periodic orbits
We provide some examples focused on the computation of periodic orbits.
Here is one where we present the different ways to compute periodic orbits. 
```@contents
Pages = ["ode/tutorialsODE-PD.md"]
Depth = 1
```

In the next tutorial, we show how to refine a periodic orbit guess obtained from numerical simulation. We also show how to perform **continuation of PD/NS** points using Shooting or Collocation. 

```@contents
Pages = ["ode/tutorialsCodim2PO.md"]
Depth = 1
```

In the next tutorial, we showcase the detection of **Chenciner** bifurcations. This is a relatively advanced tutorial, so we don't give much explanations. The reader should get first familiar with the above simpler examples.

```@contents
Pages = ["ode/steinmetz.md",]
Depth = 1
```

In the next tutorial, we showcase aBS from Bautin/HH to curve of Fold/NS of periodic orbits.

```@contents
Pages = ["ode/lorenz84-PO.md",]
Depth = 1
```

### Homoclinic orbits

Based on the package [HclinicBifurcationKit.jl](https://github.com/bifurcationkit/HclinicBifurcationKit.jl) and its [docs](https://bifurcationkit.github.io/HclinicBifurcationKit.jl/dev/).

- [Autonomous electronic circuit (aBS from BT)](https://bifurcationkit.github.io/HclinicBifurcationKit.jl/dev/tutorials/ode/tutorialsFreire/#Autonomous-electronic-circuit-(aBS-from-BT))
- [Nonlinear laser model](https://bifurcationkit.github.io/HclinicBifurcationKit.jl/dev/tutorials/ode/OPL/#Nonlinear-laser-model)

## DAE examples

```@contents
Pages = ["ode/Colpitts.md"]
Depth = 1
```

## DDE examples

See the [tutorials](https://bifurcationkit.github.io/DDEBifurcationKit.jl/dev/tutorials/tutorials/) of [DDEBifurcationKit.jl](https://github.com/bifurcationkit/DDEBifurcationKit.jl).

## Examples based on ModelingToolkit

```@contents
Pages = ["ode/NME-MTK.md"]
Depth = 1
```

## PDEs: bifurcations of equilibria
```@contents
Pages = ["tutorials1.md", "tutorials1b.md", "tutorials2.md", "mittelmann.md", "tutorials2b.md", "tutorialsSH3d.md"]
Depth = 1
```

## PDEs: automatic bifurcation diagram
```@contents
Pages = ["Swift-Hohenberg1d.md", "tutorialCarrier.md", "ks1d.md", "mittelmannAuto.md", "ks2d.md"]
Depth = 1
```

## PDEs: bifurcations of periodic orbits
```@contents
Pages = ["tutorials3.md","tutorials3b.md", "BrusselatorFF.md", "tutorialsPD.md", "tutorialsCGL.md", "tutorialsCGLShoot.md","Langmuir.md"]
Depth = 1
```

## PDEs based on FEM with [Gridap.jl](https://github.com/gridap/Gridap.jl)
```@contents
Pages = ["mittelmannGridap.md"]
Depth = 1
```

## Symmetries, freezing, waves, fronts

```@contents
Pages = ["autocatalyticAuto.md", "autocatalytic.md", "cgl1dwave.md", "detonationEngine.md"]
Depth = 1
```
