# Freezing problems, symmetries and waves

This section is dedicated to the study of the equation `F(x,p)=0` where one wishes to freeze a continuous symmetry. When the equation $F(x, p) = 0$ has a continuous symmetry described by a Lie group $G$ and action $g\cdot x$ for $g\in G$, one can reduce the symmetry of `F` by considering the constrained problem[^Beyn]:

$$\left\{
\begin{array}{l}\tag{W}
F(x, p) - s\cdot T\cdot x=0 \\
\langle T\cdot x_{ref},x-x_{ref}\rangle=0
\end{array}\right.$$

where $T$ is a generator of the Lie algebra associated to $G$, $x_{ref}$ is a reference solution and $s$ is the speed. This is known as the *freezing method*.

Similarly, one can reduce several symmetries by considering

$$\left\{
\begin{array}{l}
F(x, p) - \sum\limits_{i=1}^{N_g}\ s_i\cdot T_i\cdot x=0 \\
\langle T_i\cdot x_{ref},x-x_{ref}\rangle=0,\quad i=1,\cdots,N_g.
\end{array}\right.$$

> We call a solution `(x,p,s)` of (W) a **wave**.

## Wave stability

There are several ways to compute the stability of a wave $(x^w,p,s^w)$. From [^Sandstede], this requires to compute the spectrum of

$$d_1F(x,p)- \sum\limits_{i=1}^{N_g}\ s_i\cdot T_i\tag{EV}.$$

However, there is (potentially) the zero eigenvalue associated to the eigenvectors $T_i\cdot x^w$. In practice, because the symmetry is discrete numerically, we find a small eigenvalue instead of 0.

Another way to compute the same spectrum is to proceed as follows. Using (W) as a definition of the functional $G((x,s),p)\in\mathbb R^{N+1}$, the eigenproblem for computing the stability of a wave $(x^w,s^w)$ is

$$A x = σ Bx\tag{GEV}$$

where $B = diag(1,\cdots,1,0)$ and $A:=dG$. An advantage of (GEV) over (EV) is that the trivial eigenvalues are removed but it comes at an increased cost. We can improved this situation as follows.

## Case $N_g=1$
Let us have a look at (GEV) more closely. We need to solve for the eigenvalues $\sigma$ and the eigenvectors $(x_1,c_1)$ solutions of

$$\left\{
\begin{array}{l}\tag{W}
J x_1+c_1A_{12} = \sigma x_1 \\
\langle A_{21},x_1\rangle + A_{22}c_1=0
\end{array}\right.$$

### Case $A_{22}\neq 0$

If $A_{22}\neq 0$, the eigen problem is equivalent to

$$Jx_1 - c_1\frac{\langle A_{21},x_1\rangle}{A_{22}} A_{12}= \sigma x_1$$

### Case $A_{22} = 0$

If $A_{22} = 0$, the eigen problem is equivalent to $x_1=α A_{21} + x_1^\bot$ with $\langle A_{21},x_1^\bot\rangle=0$. Hence, I find $\langle A_{21},Jx_1^\bot\rangle+c_1\langle A_{21},A_{12}\rangle=0$

$$Jx_1^\bot-\frac{\langle A_{21},Jx_1^\bot\rangle}{\langle A_{21},A_{12}\rangle}A_{21}=σ x_1^⊥$$

## Encoding of the functional for the freezed problem

The freezing method is encoded in the composite type [`TWProblem`](@ref) which we loosely refer to as a Travelling Wave (TW) problem. 

## Computation with `newton`

We provide a simplified call to `newton` to locate the freezed solution

```
newton(prob::TWProblem, orbitguess, options::NewtonPar; kwargs...)
```

## Continuation

We also provide a simplified call to `continuation` to continue the freezed solution as function of a parameter:

```
continuation(prob::TWProblem, orbitguess, lens::Lens, contParams::ContinuationPar; jacobian = :MatrixFree, kwargs...)
```

Note that in this case, the eigen solver passed in `contParams` is converted into an appropriate generalized eigensolver.

## References

[^Beyn]:> Beyn and Thümmler, **Phase Conditions, Symmetries and PDE Continuation.**

[^Sandstede]:> Sandstede, Björn. “Stability of Travelling Waves.” In Handbook of Dynamical Systems, 2:983–1055. Elsevier, 2002. https://doi.org/10.1016/S1874-575X(02)80039-X.
