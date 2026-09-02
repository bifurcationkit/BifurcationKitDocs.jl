# Freezing problems, symmetries and waves

```@contents
Pages = ["intro_wave.md"]
Depth = 3
```

This section is dedicated to the study of the equation `F(x,p)=0` where one wishes to freeze a continuous symmetry. When the equation $F(x, p) = 0$ has a continuous symmetry described by a Lie group $G$ and action $g\cdot x$ for $g\in G$, one can reduce the symmetry of `F` by considering the constrained problem[^Beyn]:

$$\left\{
\begin{array}{l}\tag{W}
F(x, p) - s\cdot T\cdot x=0 \\
\langle T\cdot x_{ref},x-x_{ref}\rangle=0
\end{array}\right.$$

where $T$ is a generator of the Lie algebra associated to $G$ (e.g. the spatial derivative $\partial_x$ for the translation group in a 1d problem), $x_{ref}$ is a reference solution and $s$ is the speed. This is known as the *freezing method*.

Similarly, one can reduce several symmetries by considering

$$\left\{
\begin{array}{l}
F(x, p) - \sum\limits_{i=1}^{N_g}\ s_i\cdot T_i\cdot x=0 \\
\langle T_i\cdot x_{ref},x-x_{ref}\rangle=0,\quad i=1,\cdots,N_g.
\end{array}\right.$$

> We call a solution `(x,p,s)` of (W) a **wave**.

## The functional and its unknowns

The constrained system (W) is solved for the couple $(x,s)$ (the parameter $p$ plays its usual role of continuation parameter). Given the generators $T_i$, the functional $G((x,s),p)$ of (W) is

$$G((x,s),p):=\begin{bmatrix}
F(x,p) - \sum_i s_i T_i x\\
\langle x - x_{ref}, T_1 x_{ref}\rangle\\
\vdots\\
\langle x - x_{ref}, T_{N_g} x_{ref}\rangle
\end{bmatrix}\in\mathbb R^{N+N_g}.$$

It is encoded in the composite type [`TWModel`](@ref) which we loosely refer to as a Travelling Wave (TW) problem. The tutorials (e.g. on the 1d Ginzburg-Landau equation or on reaction-diffusion front problems) show how to build such a problem and to continue traveling waves.

## Wave stability

There are several ways to compute the stability of a wave $(x^w,p,s^w)$. 

- From [^Sandstede], this requires to compute the spectrum of
$$d_1F(x^w,p) - \sum\limits_{i=1}^{N_g}\ s_i\cdot T_i \tag{EV}.$$
However, there is (potentially) the zero eigenvalue associated to the eigenvectors $T_i\cdot x^w$ (the direction of the symmetry). In practice, because the symmetry is discrete numerically, we find a small eigenvalue instead of 0, it can be removed by setting `ContinuationPar.tol_stability`. This is implemented in the solver [`BifurcationKit.EigenWave`](@ref).

- Another way to compute the same spectrum is to proceed as follows. Using (W) as a definition of the functional $G((x,s),p)\in\mathbb R^{N+N_g}$, the eigenproblem for computing the stability of a wave $(x^w,s^w)$ is the *generalized* eigenproblem
$$A x = \sigma Bx\tag{GEV}$$
where $A:=dG(x^w,s^w)$ and $B = diag(1,\cdots,1,0)$ in the case $N_g=1$ (in general, the last $N_g$ entries of $B$ vanish because the equations of the constraints do not contain time derivatives). An advantage of (GEV) over (EV) is that the trivial eigenvalues are removed but it comes at an increased cost. This is implemented in the solver [`BifurcationKit.GEigenWave`](@ref). We can improve this situation as follows.

## Case $N_g=1$
Let us have a look at (GEV) more closely. Writing the linearized operator of the functional of (W) in block form

$$dG(x^w,s^w)=\begin{bmatrix} J & A_{12}\\ A_{21}^{\top} & A_{22}\end{bmatrix},\qquad J:=d_1F(x^w,p) - s^w T,$$

we need to solve for the eigenvalues $\sigma$ and the eigenvectors $(x_1,c_1)$ solutions of

$$\left\{
\begin{array}{l}
J x_1+c_1A_{12} = \sigma x_1 \\
A_{21}^{\top} x_1 + A_{22}c_1=0
\end{array}\right.$$

### Case $A_{22}\neq 0$

If $A_{22}\neq 0$, the constraint gives $c_1 = -A_{21}^{\top}x_1/A_{22}$ and the eigenproblem reduces to the (smaller) standard one

$$\left(J - \frac{A_{12} A_{21}^{\top}}{A_{22}}\right)x_1 = \sigma x_1.$$

### Case $A_{22} = 0$

If $A_{22} = 0$, the constraint imposes $A_{21}^{\top}x_1=0$. Testing the first equation of the eigenproblem against $A_{21}$ (and using $A_{21}^{\top}x_1=0$) gives $c_1\langle A_{21},A_{12}\rangle=-A_{21}^{\top}J x_1$, and the eigenproblem reduces to

$$Jx_1-\frac{\langle A_{21},Jx_1\rangle}{\langle A_{21},A_{12}\rangle}A_{12}=σ x_1,\qquad A_{21}^{\top}x_1=0.$$

For the frozen problem (W) one has $A_{12}=-T x^w$ and $A_{21}\approx T x_{ref}$, so $A_{12}$ and $A_{21}$ are both aligned with the direction $T x^w$ of the symmetry: the reduced problem above is the (EV) problem restricted to the hyperplane orthogonal to that direction, the rank-one term being the corresponding Schur complement.

## Encoding of the functional for the frozen problem

The freezing method is encoded in the composite type [`TWModel`](@ref). In particular, the user must provide the vector field `F`, the generators `(T_1,\cdots,T_{N_g})` (matrices or linear operators) and a guess `u0` for the reference solution $x_{ref}$:

```julia
probTW = TWModel(prob_vf, (∂,), u0)
```

where `prob_vf` is an `AbstractBifurcationProblem`. The reference solution $x_{ref}$ (resp. the vectors $T_i x_{ref}$ entering the constraints) is automatically updated during continuation, see [`TWModel`](@ref) for details.

## Computation with `newton`

We provide a simplified call to `newton` to locate the frozen solution

```
newton(prob::TWModel, orbitguess, options::NewtonPar; kwargs...)
```

where `orbitguess` is a guess `(x, s)` for the wave, i.e. the state together with the speeds appended at the end.

## Continuation

We also provide a simplified call to `continuation` to continue the frozen solution as function of a parameter:

```
continuation(prob::TWModel, orbitguess, alg::BifurcationKit.AbstractContinuationAlgorithm,
             contParams::ContinuationPar;
             eigsolver = BifurcationKit.GEigenWave(), kwargs...)
```

Note that in this case, the eigen solver passed in `contParams` is converted into an appropriate generalized eigensolver (see [Wave stability](@ref)); the default `GEigenWave` implements the (GEV) point of view while passing `eigsolver = BifurcationKit.EigenWave()` switches to the (EV) point of view.

## References

[^Beyn]:> Beyn and Thümmler, **Phase Conditions, Symmetries and PDE Continuation.**

[^Sandstede]:> Sandstede, Björn. “Stability of Travelling Waves.” In Handbook of Dynamical Systems, 2:983–1055. Elsevier, 2002. https://doi.org/10.1016/S1874-575X(02)80039-X.
