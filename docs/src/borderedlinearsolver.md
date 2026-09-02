# Bordered linear solvers (BLS)

```@contents
Pages = ["borderedlinearsolver.md"]
Depth = 3
```

> The bordered linear solvers must be subtypes of `AbstractBorderedLinearSolver <: AbstractLinearSolver`.

The methods provided here solve **bordered linear systems**, i.e. linear systems whose matrix is an operator $A$ augmented with one (or several) rows and columns. Such systems appear at several places in `BifurcationKit`:

- in pseudo-arclength continuation (see [`PALC`](@ref)), where the Newton correction is sought in the enlarged space $\mathbb R^{n+1}$ for the couple `(state, parameter)`,
- in the continuation / refinement of codimension 2 bifurcation points (Fold, Hopf, ...) using minimally augmented systems (see [Fold / Hopf Continuation](@ref)),
- in the computation of the normal forms at these points.

The detailed information for each one of them is located in the [API](@ref Library-BLS).

We look for the solution $u$ of $J\cdot u = v$ where

$$\tag E J=\left(\begin{array}{ll}
{A} & {b} \\
{c^T} & {d}
\end{array}\right) \text { and } v=\left(\begin{array}{l}
{v_1} \\
{v_2}
\end{array}\right)$$

where $A$ is the (large) operator to invert, $b$, $c$ are vectors and $d$ is a scalar (or, more generally, a small matrix when several constraints are appended). In practice, the bordered solvers are callable objects whose arguments depend on the context in which they are used. For example, in the pseudo-arclength continuation the system solved is

$$\left(\begin{array}{cc}
J & dR \\
\xi_u\,\mathrm{dzu}^T & \xi_p\,\mathrm{dzp}
\end{array}\right)
\left(\begin{array}{c}
dX \\ dl
\end{array}\right)
=
\left(\begin{array}{c}
R \\ n
\end{array}\right)$$

where `dzu` and `dzp` encode the tangent to the branch while $\xi_u$, $\xi_p$ are scalar weights (for the pseudo-arclength constraint, $\xi_u = \theta/|du|$ and $\xi_p = 1-\theta$). The corresponding call is

```julia
dX, dl, success, itnumber = bls(J, dR, dzu, dzp, R, n, ξu, ξp; kwargs...)
```

The keyword arguments allow one to pass a `shift` (so that the top-left block reads $\mathrm{shift}\cdot I + J$), an inner product `dotp` (used for non standard scalar products on the state space) and `applyξu!`. This is the interface implemented by all the bordered solvers below. The simple two-by-two block form $(E)$ above is recovered by taking $A=J$, $b=dR$, $c^T=\xi_u\,\mathrm{dzu}^T$, $d=\xi_p\,\mathrm{dzp}$, $v_1=R$ and $v_2=n$.

!!! warning "Complex numbers"
    In the case where $c\in\mathbb C^N$, please note that the adjoint operator $c^T$ involves a conjugate.

## Full matrix `MatrixBLS`

The easiest way to solve $(E)$ is by forming the full matrix $J$ and inverting it with the backslash operator `\` (in case it is sparse, this is relatively efficient). You can create such bordered linear solver using `bls = MatrixBLS()`; the optional argument `bls = MatrixBLS(ls)` records the underlying linear solver `ls::AbstractLinearSolver` so that the bordered solver follows the linear solver of the problem. This is the **default** bordered solver used by the pseudo-arclength continuation (`PALC`) and by the minimally augmented algorithms. It is robust and well adapted to ODE problems, but it requires assembling the $(n+1)\times(n+1)$ matrix so it may be too memory hungry for large scale problems.

## Bordering method `BorderingBLS`

The general solution $u=\left(\begin{array}{l}
{u_1} \\
{u_2}
\end{array}\right)$ to $(E)$ when $A$ is non singular is

$$\begin{array}{l}
u_2 = \frac{1}{d - c\cdot x_2}(v_2 - c\cdot x_1) \\
u_1=x_1-u_2x_2
\end{array}$$

where $x_1=A^{-1}v_1, x_2=A^{-1}b$.

It is very efficient for large scale problems because it is entirely Matrix-Free and one can use preconditioners: the two linear solves with $A$ (for $x_1$ and $x_2$) can be performed with any (preconditioned, possibly Matrix-Free) linear solver `ls`. You can create such bordered linear solver using `bls = BorderingBLS(ls)` where `ls::AbstractLinearSolver` is a linear solver which defaults to `\`.

> 1. In the case where `ls = DefaultLS()`, the factorisation of `A` is cached so the second linear solve is very fast
> 2. This method requires $A$ to be invertible. When $A$ is singular but the full bordered matrix is not (which is precisely the situation encountered at a bifurcation point), use `MatrixFreeBLS` below.

There are more options to `BorderingBLS`. First, the residual can be checked using the option `check_precision = true`. If the residual is above a prescribed tolerance `tol`, an iterative method is used based on several bordering transformations, the number of recursions being set by the option `k`. This is the *BEC+k* algorithm in [^Govaerts].

## Full Matrix-Free `MatrixFreeBLS`

In cases where $A$ is singular but $J$ is not, the bordering method may fail. It can thus be advantageous to apply a generic (Matrix-Free) linear solver to the *full* bordered operator. You can create such bordered linear solver using `bls = MatrixFreeBLS(ls)` where `ls::AbstractLinearSolver` is a (Matrix Free) linear solver which is used to invert the full bordered matrix, e.g. `GMRESIterativeSolvers` or `GMRESKrylovKit`. By default the augmented unknown `(x, p)` is stored in a `BorderedArray`; since some linear solvers only handle `AbstractArray`s (like `GMRESIterativeSolvers`), the storage can be changed to a plain `Vector` `vcat(x, p)` by passing `use_bordered_array = false`: `bls = MatrixFreeBLS(ls, false)`.

> This bordered solver requires that the full bordered operator, not just $A$, can be handled by `ls`.


## References

[^Govaerts]:> Govaerts, W. “Stable Solvers and Block Elimination for Bordered Systems.” SIAM Journal on Matrix Analysis and Applications 12, no. 3 (July 1, 1991): 469–83. https://doi.org/10.1137/0612034.
