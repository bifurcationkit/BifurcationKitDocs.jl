# Fold / Hopf Continuation

```@contents
Pages = ["codim2Continuation.md"]
Depth = 2
```

In this page, we explain how to perform continuation of Fold / Hopf points and detect the associated bifurcations.

For this to work best, it is advised to have an analytical expression for the jacobian. See the tutorial [Temperature model](@ref temperature) for more details although `BifurcationProblem` implements it with AD by default.

A quite complete example for detection of codim 2 bifurcations of equilibria is [Extended Lorenz-84 model (codim 2 + BT/ZH aBS)](@ref lorenz).

### List of detected codim 2 bifurcation points
|Bifurcation|symbol used|
|---|---|
| Bogdanov-Takens | bt |
| Bautin | gh |
| Cusp | cusp |
| Zero-Hopf | zh |
| Hopf-Hopf | hh |

In a nutshell, all you have to do (see below) is to call `continuation(br, ind_bif, lens2)` to continue the bifurcation point stored in `br.specialpoint[ind_bif]` and set proper options.

## Fold continuation (theory)

The continuation of Fold bifurcation points is based on a **Minimally Augmented**[^Govaerts] formulation which is an efficient way to detect singularities. The continuation of Fold points is based on the formulation

$$G(u,p) = (F(u,p), \sigma(u,p))\in\mathbb R^{n+1}\quad\quad (F_f)$$

where the test function $\sigma$ is solution of

$$\left[\begin{array}{cc}
dF(u,p) & a \\
b^{\top} & 0
\end{array}\right]\left[\begin{array}{c}
v \\
\sigma(u,p)
\end{array}\right]=\left[\begin{array}{c}0_{n} \\1\end{array}\right]\quad\quad (M_f)$$

where $a,b$ are chosen in order to have a non-singular matrix $(M_f)$. More precisely, $b$ (resp. $a$) should be close to a null vector of `dF(u,p)` (resp. `dF(u,p)'`). During continuation, the vectors $a,b$ are updated so that the matrix $(M_f)$ remains non-singular ; this is controlled with the argument `update_minaug_every_step` (see below).

> note that there are very simplified calls for this, see **Newton refinement** below. In particular, you don't need to set up the Fold Minimally Augmented problem yourself. This is done in the background.

!!! warning "Linear Method"
    You can pass the bordered linear solver to solve $(M_f)$ using the option `bdlinsolver ` (see below). Note that the choice `bdlinsolver = BorderingBLS()` can lead to singular systems. Indeed, in this case, $(M_f)$ is solved by inverting `dF(u,p)` which is singular at Fold points.

### Detection of codim 2 bifurcation points

You can detect the following codim 2 bifurcation points by using the option `detect_codim2_bifurcation` in the method `continuation`. Under the hood, the detection of these bifurcations is done by using Event detection as explained in [Event Handling](@ref).

- the detection of Cusp (Cusp) is done by the detection of Fold bifurcation points along the curve of Folds by monitoring the parameter component of the tangent.
- the detection of Bogdanov-Takens (BT) is performed using the test function[^Bindel] $\psi_{BT}(p) = \langle a(p),b(p)\rangle$
- the detection of Zero-Hopf (ZH) is performed by monitoring the number of eigenvalues $\lambda$ such that $\Re\lambda > \min\limits_{\nu\in\Sigma(dF)}|\Re\nu|$ and $\Im\lambda > \epsilon$ where $\epsilon$ is the Newton tolerance.

## Hopf continuation (theory)

The continuation of Hopf bifurcation points is based on a **Minimally Augmented** (see [^Govaerts] p. 87) formulation which is an efficient way to detect singularities for the Cauchy problem
$$M(u, p)\frac{du}{dt} = F(u,p)$$
where $M$ is a mass matrix. By default, it is $M=I_n$.
The continuation of Hopf points is based on the formulation

$$G(u,\omega,p) = (F(u,\omega,p), \Re\sigma(u,\omega,p), \Im\sigma(u,\omega,p))\in\mathbb R^{n+2}\quad\quad (F_h)$$

where the test function $\sigma$ is solution of

$$\left[\begin{array}{cc}
dF(u,p)-i\omega M(u,p) & M(u,p)\cdot a \\
(M(u,p)\cdot b)^{\top} & 0
\end{array}\right]\left[\begin{array}{c}
v \\
\sigma(u,\omega,p)
\end{array}\right]=\left[\begin{array}{c}
0_{n} \\
1
\end{array}\right]\quad\quad (M_h)$$

where $a,b$ are chosen in order to have a non-singular matrix $(M_h)$. More precisely, $a$ (resp. $b$) should be a left (resp. right) approximate null vector of $dF(u,p)-i\omega M(u,p)$. During continuation, the vectors $a,b$ are updated so that the matrix $(M_h)$ remains non-singular ; this is controlled with the argument `update_minaug_every_step ` (see below).

> note that there are very simplified calls to this, see **Newton refinement** below. In particular, you don't need to set up the Hopf Minimally Augmented problem yourself. This is done in the background.

!!! warning "Linear Method"
    You can pass the bordered linear solver to solve $(M_h)$ using the option `bdlinsolver ` (see below). Note that the choice `bdlinsolver = BorderingBLS()` can lead to singular systems. Indeed, in this case, $(M_h)$ is solved by inverting `dF(u,p)-iω M` which is singular at Hopf points.

!!! danger "Mass matrix $M$"
    For now, the package only deals with the case $M = I_n$.

### Detection of codim 2 bifurcation points

You can detect the following codim 2 bifurcation points by using the option `detect_codim2_bifurcation` in the method `continuation`. Under the hood, the detection of these bifurcations is done by using Event detection as explained in [Event Handling](@ref).

- the detection of Bogdanov-Takens (BT) is performed using the test function[^Bindel],[^Blank] $\psi_{BT}(p) = 	\langle a(p),b(p)\rangle$
- the detection of Bautin (GH) is based on the test function $\psi_{GH}(p) = \Re(l_1(p))$ where $l_1$ is the Lyapunov coefficient defined in [Simple Hopf point](@ref).
- the detection of Zero-Hopf (ZH) is performed by monitoring the eigenvalues.
- the detection of Hopf-Hopf (HH) is performed by monitoring the eigenvalues.

> The continuation of Hopf points is stopped at BT and when $\omega<100\epsilon$ where $\epsilon$ is the newton tolerance.

## [Setting the jacobian](@id jac-fold)

In order to apply the newton algorithm to $F_f$ or $F_h$, one needs to invert the jacobian. This is not completely trivial as one must compute this jacobian and then invert it. You can select the following jacobians for your computations (see below):

- `jacobian_ma = AutoDiff()` [Default]: automatic differentiation is applied to $F_f$ (or $F_h$) and the jacobian matrix is then inverted using the provided linear solver. In particular, the jacobian is formed. This is very well suited for small dimensions (say < 100) but **quite slow** in large dimensions because the jacobian matrix is stored and inverted (with a dense or sparse LU factorization by default).
- `jacobian_ma = MinAug()`: a specific procedure for evaluating the jacobian of $F_f$ (or $F_h$) and inverting it **without forming the jacobian matrix** is used. The bordered systems $(M_f)$ / $(M_h)$ are inverted using the linear solver of the underlying vector field (see below). This is the **recommended choice for large dimensions** and it works with matrix-free (*e.g.* iterative) solvers and preconditioners.
- `jacobian_ma = MinAugMatrixBased()`: the jacobian matrix is evaluated using an analytical formula, which allows for example to form a sparse matrix when the underlying problem has a sparse jacobian. It is faster than `AutoDiff()` and a good middle-ground when a sparse (or banded) jacobian is available.
- `jacobian_ma = FiniteDifferencesMF()`: the jacobian is evaluated in a matrix-free version using finite differences. **Mainly for debugging purposes**.
- `jacobian_ma = FiniteDifferences()`: the jacobian matrix is evaluated using finite differences. **Mainly for debugging purposes.**

> When `jacobian_ma = MinAug()` (or `MinAugMatrixBased()`) is used on a non-symmetric problem, the **adjoint** of the jacobian of the vector field is required to build the bordered systems. If you do not provide one (through the option `Jᵗ` of the bifurcation problem, or more generally through `jacobian_adjoint`), it is computed internally using `transpose(J)` which works for `AbstractArray`. For matrix-free jacobians you must provide the adjoint yourself (e.g. the CGL tutorial relies on this for the Hopf continuation, see below). See also the tips in [`newton_fold`](@ref) / [`newton_hopf`](@ref).

## [Linear solvers & preconditioners for large scale](@id ls-large)

Even with `jacobian_ma = MinAug()`, the linear systems associated to the vector field jacobian $J = dF(u,p)$ (or $J-i\omega I$) still have to be solved at each Newton step, see $(M_f)$ / $(M_h)$ above. In large dimensions one must therefore choose carefully

- the **linear solver** `options.newton_options.linsolver` used to invert these systems. The default `DefaultLS()` (based on LU/Cholesky) becomes prohibitive for large sparse systems. It should be replaced by a **preconditioned iterative solver**, e.g. `GMRESIterativeSolvers` or `GMRESKrylovKit` with an ILU or AMG preconditioner, see [Linear solvers (LS)](@ref). See the tutorial [2d Ginzburg-Landau equation (finite differences, codim 2, Hopf aBS)](@ref cgl) for a complete example based on ILU preconditioners.
- the **bordered linear solver** `bdlinsolver` used to solve $(M_f)$ / $(M_h)$. For large scale problems, `MatrixFreeBLS` or `BorderingBLS` (with a good preconditioned solver) should be preferred over `MatrixBLS` which forms the full bordered matrix, see [Bordered linear solvers (BLS)](@ref). Note that the bordering strategy can fail close to a Fold / Hopf point since it requires inverting the singular operator $J$ (resp. $J-i\omega I$), see the warning above.

!!! tip "Preconditioner update"
    It can be advantageous to recompute the preconditioner during the continuation, e.g. every few steps, using the `callback_newton` mechanism. See the example [2d Ginzburg-Landau equation (finite differences, codim 2, Hopf aBS)](@ref cgl).

!!! tip "Hopf and non-symmetric problems"
    For Hopf continuation the bordered system involves both $J-i\omega I$ and its adjoint. When the jacobian is not symmetric (or is matrix-free), you can pass dedicated solvers through the options `linsolve_adjoint` and `bdlinsolver_adjoint` of `continuation` (see [`continuation_hopf`](@ref)). In the tutorial [2d Ginzburg-Landau equation (finite differences, codim 2, Hopf aBS)](@ref cgl), the option `start_with_eigen = true` is used for the same reason: the left eigenvector of the jacobian is not the conjugate of the right one.

## [start_with_eigen, update_minaug_every_step & compute_eigen_elements](@id start-with-eigen)

Some options deserve special attention in large dimensions:

- `start_with_eigen = true`: initializes the vectors $a, b$ of the Minimally Augmented formulation from the eigen-elements (right and left eigenvectors) of the bifurcation point stored in `br`. This is **recommended**, especially for Hopf continuation where the left eigenvector of the (possibly non-symmetric) jacobian is not the conjugate of the right one. It is also useful for Fold continuation as it removes the need for an initial bordered solve with random vectors.
- `update_minaug_every_step`: controls how often the vectors $a, b$ are recomputed so that the matrix $(M_f)$ / $(M_h)$ remains well-conditioned. Keeping the default value `= 1` (update at every continuation step) is strongly recommended; when detecting codim 2 bifurcations, a warning is issued if you set it to `0` because the detection may then be unreliable.
- `compute_eigen_elements`: whether to compute eigenelements along the curve of Fold / Hopf points. It is required for the detection of Zero-Hopf (ZH) / Hopf-Hopf (HH) points, see [Event Handling](@ref).

## Recap for large dimensions

For the continuation of Fold / Hopf points in large dimensions, a typical call is

```julia
const BK = BifurcationKit

# preconditioned iterative solver for the linear systems
ls = GMRESIterativeSolvers(reltol = 1e-4, N = n, Pl = ilu(J0, τ = 0.005)) # or GMRESKrylovKit
opts = ContinuationPar(br.contparams; newton_options = NewtonPar(linsolver = ls))

brfold = continuation(br, ind_bif, lens2, opts;
    # matrix-free evaluation of the jacobian of the Fold / Hopf functional
    jacobian_ma = BK.MinAug(),
    # bordered linear solver for the MA formulation
    bdlinsolver = BorderingBLS(solver = ls, check_precision = false),
    # recommended, esp. for Hopf
    start_with_eigen = true,
    # keep the MA vectors up to date
    update_minaug_every_step = 1,
    # detect codim 2 bifurcations (BT, CP, ZH, ...)
    detect_codim2_bifurcation = 2,
    normC = norminf)
```

You can find working examples in the tutorials [Temperature model (codim 2)](@ref temperature), [Extended Lorenz-84 model (codim 2 + BT/ZH aBS)](@ref lorenz), [1d Brusselator (automatic)](@ref brusauto), [1d Langmuir–Blodgett transfer model](@ref langmuir) and, for the most advanced matrix-free settings with preconditioners, [2d Ginzburg-Landau equation (finite differences, codim 2, Hopf aBS)](@ref cgl).

## Newton refinement

Once a Fold / Hopf point has been detected after a call to `br = continuation(...)`, it can be refined using `newton` iterations. Let us say that `ind_bif` is the index in `br.specialpoint` of a Fold / Hopf point. This guess can be refined as follows:

```julia
outfold = newton(br::AbstractBranchResult, ind_bif::Int;  
	normN = norm, options = br.contparams.newton_options,
	start_with_eigen = false,
	lens2 = (@optic _), kwargs...)
```

For the options parameters, we refer to [Krylov-Newton algorithm](@ref). Note that you can pass a bordered linear solver through the option `bdlinsolver`, see [`newton_fold`](@ref) / [`newton_hopf`](@ref). In large dimensions, we recommend using `start_with_eigen = true` (see [above](@ref start-with-eigen)) together with the preconditioned iterative solver setup described in [Linear solvers & preconditioners for large scale](@ref ls-large).

It is important to note that for improved performances, a function implementing the expression of the **hessian** should be provided. This is by far the fastest. `BifurcationProblem` provides it by default using AD though. In the matrix-free case (`jacobian_ma = MinAug()`), the hessian is optional: if not provided, the derivative $\partial_x\sigma$ is evaluated with finite differences (see [Algorithmic details (Fold)](@ref)).

Reader interested in this advanced usage should look at the code `example/chan.jl` of the tutorial [Temperature model](@ref temperature).

## Fold / Hopf continuation

To compute the curve of Fold / Hopf points, one can call [`continuation`](@ref) with the following options

```@docs
 continuation(br::BifurcationKit.AbstractBranchResult, ind_bif,
				lens2::BifurcationKit.AllOpticTypes, options_cont::ContinuationPar = br.contparams ;
				kwargs...)
```

where the options are as above except that we have an additional parameter axis `lens2` which is used to locate the bifurcation points.


See [Temperature model](@ref temperature) for an example of use.

## Advanced use

Here, we expose the solvers that are used to perform newton refinement or codim 2 continuation in case the above methods fails. This is useful in case it is too involved to expose the linear solver options. An example of advanced use is the continuation of Folds of periodic orbits, see [Continuation of Fold of periodic orbits](@ref fold-po).

```@docs
BifurcationKit.newton_fold
```

```@docs
BifurcationKit.newton_hopf
```


```@docs
continuation_fold
```

```@docs
continuation_hopf
```

## Algorithmic details (Fold)

Here we detail the computation of the jacobian of the Fold functional $G=(F,\sigma)$ required by the Newton algorithm. During the differentiation, the bordering vectors $a,b$ are **kept fixed** (they are only updated between continuation steps, see [start_with_eigen, update_minaug_every_step & compute_eigen_elements](@ref start-with-eigen)). Let $J(u,p)=\partial_uF(u,p)$ and consider the bordered system

$$\left[\begin{array}{cc}J(u,p)&a\\ b^{\top}&0\end{array}\right]\left[\begin{array}{c}v\\ \sigma(u,p)\end{array}\right]=\left[\begin{array}{c}0_{n}\\1\end{array}\right],$$

of which $(v,\sigma)$ is the solution. Because the bordered system and its adjoint share the same right-hand side $(0_n,1)$, the second component of the adjoint solution equals the test function, i.e. $\tau=\sigma$: the **adjoint bordered system**

$$\left[\begin{array}{cc}J(u,p)^{\top}&b\\ a^{\top}&0\end{array}\right]\left[\begin{array}{c}w\\ \sigma\end{array}\right]=\left[\begin{array}{c}0_{n}\\1\end{array}\right]\quad\Longleftrightarrow\quad J^{\top}w+b\sigma=0,\quad a^{\top}w=1$$

only determines the vector $w$, close to a null-vector of $J^{\top}$, that is to a left null-vector of $J$. Differentiating the bordered system along a direction $\dot z=(\dot u,\dot p)$, multiplying the first $n$ equations by $w^{\top}$ and using $a^{\top}w=1$ together with $b^{\top}v=1$ (so that $b^{\top}\partial_z v=0$), one can show[^Govaerts] that the differential of $\sigma$ with respect to $z$ satisfies:

$$\partial_z \sigma\cdot\dot z = -w^{\top}\partial_z dF(u,p)[\dot z]\,v \quad\Longleftrightarrow\quad \partial_z \sigma + \langle w,\partial_z dF \cdot v\rangle = 0$$

This allows to compute the jacobian of the Fold functional to use for the Newton algorithm:

$$\left[\begin{array}{cc}
\partial_{u}F(u,p) & \partial_pF(u,p) \\
\partial_u\sigma(u,p) & \partial_p\sigma(u,p)
\end{array}\right],\qquad
(\partial_u\sigma)_i=-w^{\top}(\partial_{u_i}J)\,v,\quad
\partial_p\sigma=-w^{\top}(\partial_pJ)\,v ,$$

for $i=1,\dots,n$. The bottom row requires, on top of the bordered solve giving $v,\sigma$ (which is already needed to evaluate $G$), one **adjoint bordered solve** giving $w$ and the Hessian contractions $(\partial_u\sigma)_i=-\langle w,\partial^2F(u,p)[e_i,v]\rangle$, where $\partial^2F(u,p)[\cdot,\cdot]$ is the second derivative of $F$ with respect to $u$. If no Hessian is available, this row is evaluated with finite differences; in the matrix-free case this is done without ever forming the jacobian.

## Algorithmic details (Hopf)

We recall that the unknowns are $(u,p,\omega)$ and that the bordered system $(M_h)$ involves the complex matrix $A(u,p,\omega)=dF(u,p)-i\omega M(u,p)$, where the mass matrix $M=M(u,p)$ now depends on $u$ and $p$ (by default $M=I_n$). Here $(v,\sigma)$ denotes the solution of the bordered system $(M_h)$: its first $n$ components $v$ form an approximate right null-vector of $A$ (i.e. $Av\simeq0$ at the Hopf point) while $\sigma$ is the complex test function used in $(F_h)$. As a consequence, the borders $M(u,p)\,a$ and $(M(u,p)\,b)^{\top}$ of $(M_h)$ also depend on $(u,p)$: contrary to the Fold case, they can no longer be considered constant when differentiating $\sigma$.

As in the Fold case, the second component of the adjoint solution coincides with $\sigma$ (i.e. $\tau=\sigma$), so that the adjoint bordered system

$$\left[\begin{array}{cc}
A(u,p,\omega)^{\top} & M(u,p)\,b\\
(M(u,p)\,a)^{\top} & 0
\end{array}\right]\left[\begin{array}{c}w\\ \sigma\end{array}\right]=\left[\begin{array}{c}0_{n}\\ 1\end{array}\right],$$

i.e. $A^{\top}w+Mb\sigma=0$ and $(Ma)^{\top}w=1$, determines the vector $w$. Differentiating the bordered system along a direction $\dot z$ of the unknowns $(u,p,\omega)$ (the vectors $a,b$ being kept fixed, as in the Fold case) and proceeding as above yields

$$\partial_z\sigma\cdot\dot z = -w^{\top}\partial_z A(u,p,\omega)[\dot z]\,v
- \sigma\Big( w^{\top}\partial_z M(u,p)[\dot z]\,a + b^{\top}\partial_z M(u,p)[\dot z]^{\top}v\Big),$$

where the term in parentheses gathers the contributions of $\partial_zM$ through the $(u,p)$-dependent borders $Ma$ and $(Mb)^{\top}$. Since $\partial_\omega M=0$ and $\partial_\omega A(u,p,\omega)=-iM(u,p)$, this simplifies to

$$\partial_\omega\sigma = i\,w^{\top} M(u,p)\, v ,$$

while the $u$- and $p$-derivatives of $\sigma$ keep the terms involving $\partial_uM$ and $\partial_pM$:

$$\partial_{u_i}\sigma = -w^{\top}\partial_{u_i}A(u,p,\omega)\,v
- \sigma\Big( w^{\top}\partial_{u_i}M\,a + b^{\top}(\partial_{u_i}M)^{\top}v\Big) ,\qquad i=1,\dots,n,$$

$$\partial_p\sigma = -w^{\top}\partial_pA(u,p,\omega)\,v
- \sigma\Big( w^{\top}\partial_pM\,a + b^{\top}(\partial_pM)^{\top}v\Big) .$$

Since the functional $G=(F,\Re\sigma,\Im\sigma)$ uses only the real and imaginary parts of $\sigma$ and since $F$ does not depend on $\omega$, the jacobian of the Hopf functional to use for the Newton algorithm is

$$\left[\begin{array}{ccc}
\partial_{u}F & \partial_pF & 0 \\
\Re\,\partial_{u}\sigma & \Re\,\partial_{p}\sigma & \Re\,\partial_{\omega}\sigma\\
\Im\,\partial_{u}\sigma & \Im\,\partial_{p}\sigma & \Im\,\partial_{\omega}\sigma
\end{array}\right].$$

- **The case $M=I_n$.** When the mass matrix is the identity (the default), $M$ does not depend on $(u,p)$: the borders $Ma$ and $(Mb)^{\top}$ are constant again and the extra terms in $\partial_u\sigma,\partial_p\sigma$ drop out ($\partial_uM=\partial_pM=0$). The formulas then reduce to their Fold counterparts with $dF-i\omega I$ in place of $J$:
$$\partial_{u_i}\sigma=-w^{\top}\partial_{u_i}dF\,v,\qquad \partial_p\sigma=-w^{\top}\partial_pdF\,v,\qquad \partial_\omega\sigma=i\,w^{\top}v .$$

- **The case of a constant $M$.** The previous remark is a particular case of a constant mass matrix, i.e. $M$ independent of $(u,p)$ but not necessarily the identity. Then $\partial_uM=\partial_pM=0$: the borders $Ma$ and $(Mb)^{\top}$ are constant again, the extra terms drop out and the formulas reduce to their Fold counterparts with $dF-i\omega M$ in place of $J$:
$$\partial_{u_i}\sigma=-w^{\top}\partial_{u_i}dF\,v,\qquad \partial_p\sigma=-w^{\top}\partial_pdF\,v,\qquad \partial_\omega\sigma=i\,w^{\top}M\,v .$$

- **Scalar case ($n=1$).** All the quantities become scalars and the formulas above can be checked explicitly. Solving the bordered system $(M_h)$ and its adjoint gives $v=1/(bM)$, $w=1/(aM)$ and $\sigma=-A/(abM^{2})$, so that $\sigma+wAv=0$. Differentiating $\sigma$ directly with $a,b$ fixed yields (SCA)
which is exactly what the formula above predicts: the first term is $-w\partial_pA\,v$ while the border term reads $-\sigma(w\partial_pM\,a+b\partial_pM\,v)=-\sigma\big(\partial_pM/M+\partial_pM/M\big)=2A\,\partial_pM/(abM^{3})$. In particular this border term is $\propto\partial_pM/M$ and only vanishes when $\partial_pM=0$, which justifies keeping the extra $\sigma(\cdots)$ terms when the mass matrix depends on $(u,p)$.
$$\partial_p\sigma=-\frac{\partial_pA}{abM^{2}}+\frac{2A\,\partial_pM}{abM^{3}}\quad\quad (SCA).$$

## References

[^Govaerts]: > Govaerts, Willy J. F. Numerical Methods for Bifurcations of Dynamical Equilibria. Philadelphia, Pa: Society for Industrial and Applied Mathematics, 2000.

[^Blank]: > Blank, H. J. de, Yu. A. Kuznetsov, M. J. Pekkér, and D. W. M. Veldman. “Degenerate Bogdanov–Takens Bifurcations in a One-Dimensional Transport Model of a Fusion Plasma.” Physica D: Nonlinear Phenomena 331 (September 15, 2016): 13–26. https://doi.org/10.1016/j.physd.2016.05.008.

[^Bindel]: > Bindel, D., M. Friedman, W. Govaerts, J. Hughes, and Yu.A. Kuznetsov. “Numerical Computation of Bifurcations in Large Equilibrium Systems in Matlab.” Journal of Computational and Applied Mathematics 261 (May 2014): 232–48. https://doi.org/10.1016/j.cam.2013.10.034.
