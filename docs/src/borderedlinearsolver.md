# Bordered linear solvers (BLS)

> The bordered linear solvers must be subtypes of `AbstractBorderedLinearSolver <: AbstractLinearSolver`.

The methods provided here solve bordered linear equations. More precisely, one is interested in the solution $u$ to $J\cdot u = v$ where

$$\tag E J=\left(\begin{array}{ll}
{A} & {b} \\
{c^T} & {d}
\end{array}\right) \text { and } v=\left(\begin{array}{l}
{v_1} \\
{v_2}
\end{array}\right)$$

Such linear solver `bdlsolve` will be called like `sol, success, itnumber = bdlsolve(A, b, c, d, v1, v2)` throughout the package.

!!! warning "Complex numbers"
    In the case where $c\in\mathbb C^N$, please note that the adjoint operator $c^T$ involves a conjugate.

## Full matrix `MatrixBLS`
This easiest way to solve $(E)$ is by forming the matrix $J$. In case it is sparse, it should be relatively efficient. You can create such bordered linear solver using `bls = MatrixBLS(ls)` where `ls::AbstractLinearSolver` is a linear solver (which defaults to `\`) used to solve the linear problem associated to $J$. This is the default method used in the package. 

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

It is very efficient for large scale problems because it is entirely Matrix-Free and one can use preconditioners. You can create such bordered linear solver using `bls = BorderingBLS(ls)` where `ls::AbstractLinearSolver` is a linear solver which defaults to `\`. The intermediate solutions $x_1=A^{-1}v_1, x_2=A^{-1}b$ are formed using `ls`.

> 1. In the case where `ls = DefaultLS()`, the factorisation of `A` is cached so the second linear solve is very fast

There are more options to `BorderingBLS`. First, the residual can be checked using the option `check_precision = true`. If the residual is above a prescribed tolerance, an iterative method is used based on several bordering transformations. This is the *BEC+k* algorithm in [^Govaerts].

## Full Matrix-Free `MatrixFreeBLS`

In cases where $A$ is singular but $J$ is not, the bordering method may fail. It can thus be advantageous to form the Matrix-Free version of $J$ and call a generic linear solver to find the solution to $(E)$. You can create such bordered linear solver using `bls = MatrixFreeBLS(ls)` where `ls::AbstractLinearSolver` is a (Matrix Free) linear solver which is used to invert `J`.

> For now, this linear solver only works with `AbstractArray`


## References

[^Govaerts]:> Govaerts, W. “Stable Solvers and Block Elimination for Bordered Systems.” SIAM Journal on Matrix Analysis and Applications 12, no. 3 (July 1, 1991): 469–83. https://doi.org/10.1137/0612034.

