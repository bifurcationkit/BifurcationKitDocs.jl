# Tutorials

```@contents
Pages = ["tutorials.md"]
Depth = 3
```

The tutorials are rated by the following scale of difficulty

1. 🟢 basic knowledge of (numerical) bifurcation theory (following equilibria / periodic orbits)
2. 🟡 advanced knowledge of (numerical) bifurcation theory (codim 2 bifurcations of equilibria)
2. 🟠 high level of knowledge of (numerical) bifurcation theory (codim 2 bifurcations of periodic orbits, tweaking the methods)
2. 🟤 very advanced tutorial, research level

There are three levels of automatization of the computation in these tutorials:

1. fully **automatic bifurcation diagram** (**aBD**) computation (only for equilibria): one uses `bifurcationdiagram` and let it compute the diagram fully automatically. Another possibility is to use **deflated continuation**.
2. semi-automatic bifurcation diagram computation: one uses **automatic branch switching** (**aBS**) to compute branches at specified bifurcation points
3. manual bifurcation diagram computation: one does not use automatic branch switching. This has only educational purposes or for complex problems where aBS fails.


## ODE examples

These examples are specific to ODEs. 

### Computation of equilibria

```@contents
Pages = ["ode/tutorialsBasic1.md", "ode/tutorials1.md", "ode/tutorialPP2.md",]
Depth = 1
```

### Codimension 2 bifurcations of equilibria

```@contents
Pages = ["ode/tutorialCO.md","ode/lorenz84.md",]
Depth = 1
```

### Periodic orbits
We provide some examples focused on the computation of periodic orbits.
Here is one where we present the different ways to compute periodic orbits. 

```@contents
Pages = ["ode/tutorialsODE.md"]
Depth = 1
```

Here is one for aBS from **period-doubling** bifurcations of periodic orbits
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

- 🟡 [Autonomous electronic circuit (aBS from BT)](https://bifurcationkit.github.io/HclinicBifurcationKit.jl/dev/tutorials/ode/tutorialsFreire/#Autonomous-electronic-circuit-(aBS-from-BT))
- 🟡 [Nonlinear laser model](https://bifurcationkit.github.io/HclinicBifurcationKit.jl/dev/tutorials/ode/OPL/#Nonlinear-laser-model)


## DAE examples

```@contents
Pages = ["ode/Colpitts.md"]
Depth = 1
```

## DDE examples

See the [tutorials](https://bifurcationkit.github.io/DDEBifurcationKit.jl/dev/tutorials/tutorials/) of [DDEBifurcationKit.jl](https://github.com/bifurcationkit/DDEBifurcationKit.jl).

## Examples based on ModelingToolkit

`ModelingToolkit` provides a tailored interface to `BifurcationKit` and the user is encouraged to have a look at its specific [documentation](https://docs.sciml.ai/ModelingToolkit/stable/tutorials/bifurcation_diagram_computation/).

We also provide some alternative example(s): 

```@contents
Pages = ["ode/NME-MTK.md"]
Depth = 1
```
## Examples based on Catalyst

`Catalyst` provides a tailored interface to `BifurcationKit` and the user is encouraged to have a look at its specific [documentation](https://docs.sciml.ai/Catalyst/stable/steady_state_functionality/bifurcation_diagrams/)

## Multi-parameters


Based on the package [MultiParamContinuation.jl](https://github.com/bifurcationkit/MultiParamContinuation.jl) and its [docs](https://bifurcationkit.github.io/MultiParamContinuation.jl/dev/). This functionality allows to continue solutions in two parameters.

- 🟢 [sphere example](https://bifurcationkit.github.io/MultiParamContinuation.jl/dev/tutorials/ode/sphere/#Sphere)
- 🟢 [sphere example based in BifurcationKit](https://bifurcationkit.github.io/MultiParamContinuation.jl/dev/tutorials/ode/sphereBK/#Sphere-based-on-BifurcationKit.jl)

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
