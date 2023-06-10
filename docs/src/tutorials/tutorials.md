# Tutorials

There are three levels of tutorials:

1. fully **automatic bifurcation diagram** (**aBD**) computation (only for equilibria): one uses the function `bifurcationdiagram` and let it compute the diagram fully automatically. Another possibility is to use **deflated continuation**.
2. semi-automatic bifurcation diagram computation: one uses **automatic branch switching** (**aBS**) to compute branches at specified bifurcation points
3. manual bifurcation diagram computation: one does not use automatic branch switching. This has only educational purposes or for complex problems where aBS fails.

## ODE examples

We present examples in the case of ODEs. Although `BifurcationKit.jl` is not geared towards them, we provide some specific methods which allow to study the bifurcations of ODE in a relatively efficient way.

```@contents
Pages = ["ode/tutorialsODE.md","ode/tutorialCO.md","ode/lorenz84.md", "ode/tutorialPP2.md",]
Depth = 1
```

Here are some examples more oriented towards the computation of periodic orbits. Here is one for aBS from **period-doubling** bifurcations of periodic orbits
```@contents
Pages = ["ode/tutorialsODE-PD.md"]
Depth = 1
```

In the next tutorial, we show how to refine a periodic orbit guess obtained from numerical simulation. We also show how to perform **continuation of PD/NS** points using Shooting or Collocation. 

```@contents
Pages = ["ode/tutorialsCodim2PO.md"]
Depth = 1
```

In the next tutorial, we showcase the detection of **Chenciner** bifurcations.

```@contents
Pages = ["ode/steinmetz.md",]
Depth = 1
```

In the next tutorial, we showcase aBS from Bautin/HH to curve of Fold/NS of periodic orbits.

```@contents
Pages = ["ode/lorenz84-PO.md",]
Depth = 1
```

## DAE examples

```@contents
Pages = ["ode/Colpitts.md"]
Depth = 1
```

## Examples based on ModelingToolkit

```@contents
Pages = ["ode/NME-MTK.md"]
Depth = 1
```

## Bifurcation of Equilibria
```@contents
Pages = ["tutorials1.md", "tutorials1b.md", "tutorials2.md", "mittelmann.md", "tutorials2b.md", "tutorialsSH3d.md"]
Depth = 1
```

### Automatic bifurcation diagram
```@contents
Pages = ["Swift-Hohenberg1d.md", "tutorialCarrier.md", "ks1d.md", "mittelmannAuto.md", "ks2d.md"]
Depth = 1
```

### Solving PDEs using Finite elements with [Gridap.jl](https://github.com/gridap/Gridap.jl)
```@contents
Pages = ["mittelmannGridap.md"]
Depth = 1
```

## Bifurcation diagrams with periodic orbits
```@contents
Pages = ["tutorials3.md","tutorials3b.md", "BrusselatorFF.md", "tutorialsPD.md", "tutorialsCGL.md", "tutorialsCGLShoot.md","Langmuir.md"]
Depth = 1
```

## Symmetries, freezing, waves, fronts

```@contents
Pages = ["autocatalyticAuto.md", "autocatalytic.md", "cgl1dwave.md", "detonationEngine.md"]
Depth = 1
```
