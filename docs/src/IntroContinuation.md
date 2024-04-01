# Continuation methods: introduction

```@contents
Pages = ["IntroContinuation.md"]
Depth = 3
```

The goal of these methods[^Kuz],[^Govaerts],[^Rabinowitz],[^Mei],[^Keller] is to find solutions $x\in\mathbb R^n$ to nonlinear equations

$$\mathbb R^n\ni F(x,p) = 0 \quad\tag{E}$$

as function of a real parameter $p$. Given a known solution $(x_0,p_0)$, we can, under reasonable assumptions, continue it by computing a 1d curve of solutions $\gamma = (x(s),p(s))_{s\in I}$ passing through $(x_0,p_0)$.

For the sequel, it is convenient to use the following formalism [^Kuz]

1. prediction of the next point
2. correction
3. step size control.


## Natural continuation

> More information is available at [Predictors - Correctors](@ref)

We just use this simple continuation method to give a **trivial** example of the formalism.
Knowing $(x_0, p_0)$, we form the predictor $(x_0, p_0+ds)$ for some $ds$ and use it as a guess for a Newton corrector applied to $x\to F(x, p_0+ds)$. The corrector is thus the **Newton** algorithm.

This continuation method is `continuation(prob, Natural(), options)`.

!!! danger "Usage"
    You should almost never use this predictor for computations. It fails at turning points, is not adaptive, ... 

## Step size control

Each time the corrector phase failed, the step size ``ds`` is halved. This has the disadvantage of having **unnecessary** Newton iterations (which costs time) and imposing small steps (which can be slow as well). To prevent this, the step size can be controlled internally with the idea of having a constant number of Newton iterations per point. This is in part controlled by the aggressiveness factor `a` in `ContinuationPar`.


## References

[^Kuz]:> Kuznetsov, Elements of Applied Bifurcation Theory.

[^Govaerts]:> Govaerts, Numerical Methods for Bifurcations of Dynamical Equilibria; Allgower and Georg, Numerical Continuation Methods

[^Rabinowitz]:> Rabinowitz, Applications of Bifurcation Theory; Dankowicz and Schilder, Recipes for Continuation

[^Mei]:> Mei, Numerical Bifurcation Analysis for Reaction-Diffusion Equations

[^Keller]:> Keller, Lectures on Numerical Methods in Bifurcation Problems
