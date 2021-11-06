# Freezing problems, symmetries and waves

This section is dedicated to the study of an equation (in `x`) `F(x,p)=0` where one wishes to freeze a continuous symmetry. When the equation $F(x, p) = 0$ has a continuous symmetry described by a Lie group $G$ and action $g\cdot x$ for $g\in G$, one can reduce the symmetry of the problem by considering the constrained problem[^Beyn]:

$$\left\{
\begin{array}{l}
F(x, p) - s\cdot T\cdot x=0 \\
\langle T\cdot x_{ref},x-x_{ref}\rangle=0
\end{array}\right.$$

where $T$ is a generator of the Lie algebra associated to $G$, $x_{ref}$ is a reference solution and $s$ is the speed. This is known as the *freezing method*.

Similarly, one can reduce several symmetries by considering

$$\left\{
\begin{array}{l}
F(x, p) - \sum\limits_{i=1}^N\ s_i\cdot T_i\cdot x=0 \\
\langle T_i\cdot x_{ref},x-x_{ref}\rangle=0,\quad i=1,\cdots,N.
\end{array}\right.$$


## Encoding of the functional for the freezed problem

The freezing method is encoded in the composite type [`TWProblem`](@ref) which we loosely refer to ass a Travelling Wave (TW) problem.

## Computation with `newton`

We provide a simplified call to `newton` to locate the freezed solution (*e.g.* a travelling wave)

```
newton(prob::TWProblem, orbitguess, par0, options::NewtonPar; kwargs...)
```

## Continuation

We also provide a simplified call to `continuation` to continue the freezed solution as function of a parameter:

```
continuation(prob::TWProblem, orbitguess, par, lens::Lens, contParams::ContinuationPar; jacobian = :MatrixFree, kwargs...)
```

Note that in this case, the eigen solver passed in `contParams` is converted into an appropriate generalized eigenvector.

## References

[^Beyn]:> Beyn and Thümmler, **Phase Conditions, Symmetries and PDE Continuation.**
