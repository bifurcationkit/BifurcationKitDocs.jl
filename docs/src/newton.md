# Krylov-Newton algorithm

```@contents
Pages = ["newton.md"]
Depth = 3
```

`BifurcationKit` is built upon the newton algorithm for solving (large-dimensional) nonlinear equations

$$F(x)=0\in\mathbb R^n,\quad x\in\mathbb R^n.$$

Writing $J(x)\in\mathcal L(\mathbb R^n)$ the jacobian, the algorithm reads

$$x_{n+1} = x_n - J(x_n)^{-1}F(x_n)$$

with initial guess $x_0$.

The crux of the algorithm is to solve the linear system in $y$:

$$J(x_n)\cdot y = F(x_n).$$

To this end, we never form $J^{-1}$ like with `pinv(J)` but solve the linear system directly.

## Convergence criterion

The Newton iteration is stopped when the (scaled) residual becomes smaller than a user tolerance. More precisely, at iteration $k$ we compute $r_k:=N\big(F(x_k)\big)$ where $N$ is a norm that can be specified by the keyword `normN` of [`BifurcationKit.solve`](@ref) (by default `normN = norm`, the ``\ell^2`` norm). The iteration is declared **converged** as soon as

$$r_k < \mathtt{tol} \quad\text{and}\quad \text{callback}(x_k, p, F(x_k)) \text{ returns true},$$

where `tol` is the field `NewtonPar.tol` (default ``10^{-12}``) and `callback` is a keyword of `solve` which allows to reject (or stop the iteration on) steps based on extra information, e.g. the norm of the residual or the norm of the Newton step. `BifurcationKit` provides the helpers [`BifurcationKit.cbMaxNorm`](@ref) and [`BifurcationKit.cbMaxNormAndΔp`](@ref) to build such callbacks, and they can be passed to `solve`/`continuation` (as `callback`/`callback_newton` respectively). We refer to the output of `continuation` iterations for examples of use.

!!! tip "Tip"
    When the Newton iteration is used as the corrector of a continuation method, the same criterion is used but the norm is then the one given by the argument `normC` of `continuation` (see [Pseudo arclength  continuation](@ref)).

## Options of the Newton algorithm

The Newton iterations are controlled by the composite type [`NewtonPar`](@ref). The most useful fields are:

| field | default | meaning |
|:------|:--------|:--------|
| `tol` | ``10^{-12}`` | tolerance on the residual ``\|F(x)\|`` to declare convergence |
| `max_iterations` | `25` | maximal number of Newton iterations |
| `verbose` | `false` | print the iterations? |
| `linsolver` | `DefaultLS()` | linear solver used to invert the jacobian, must be `<: AbstractLinearSolver` (see [Linear solvers (LS)](@ref)) |
| `eigsolver` | `DefaultEig()` | eigen solver used to compute eigenvalues, must be `<: AbstractEigenSolver` (see [Eigen solvers (Eig)](@ref)) |
| `linesearch` | `false` | use a line search algorithm (i.e. Newton with Armijo's rule) |
| `α` | `1.0` | initial value of the damping factor used by the line search |
| `αmin` | `0.001` | minimal value of the damping factor |

!!! note "Line search"
    The fields `linesearch`, `α`, `αmin` are used by the *PALC corrector* `newton_palc` during continuation (see [Pseudo arclength  continuation](@ref)). The plain Newton algorithm described on this page does not use them.

## Space of solutions

For the algorithm to be defined, a certain number of operations on `x` need to be available. If you pass `x::AbstractArray`, you should not have any problem. Otherwise, your `x` must comply with the requirements listed in [Required methods for custom arrays](@ref Required-Methods).

## Different Jacobians

There are basically two ways to specify the jacobian:

1. Matrix based
2. Matrix-free.

In case you pass a matrix (in effect an `AbstractMatrix` like a sparse one,...), you can use the default linear solver from `LinearAlgebra` termed the backslash operator `\`. This is a **direct** method. This is the case 1 above.

Another possibility is to pass a function `J(dx)` and to use **iterative** linear solvers. In this case, this is termed a **Krylov-Newton** method. This is the case 2 above. In comparison to the Matrix-based case, there is no restriction to the number of unknowns $n$.

> The available linear solvers are explained in the section [Linear solvers (LS)](@ref).

One can find a full description of the Krylov-Newton method in the docstring of [`BifurcationKit.solve`](@ref).

## Simple example

Here is a quick example to show how the basics work. In particular, the problem generates a matrix based jacobian using automatic differentiation.

```@example NEWTON
using BifurcationKit
F(x, p) = x.^3 .- 1
x0 = rand(10)
prob = ODEBifProblem(F, x0, nothing)
sol = BifurcationKit.solve(prob, Newton(), NewtonPar(verbose = true))
```

> The function `solve` is not exported so we call it through `BifurcationKit.solve`. The type `ODEBifProblem` is a shortcut to define a problem of the form $\dot x = F(x,p)$ ; the generic (autonomous) problem type is `BifurcationProblem`. Both are equivalent for the Newton algorithm above.

The returned object `sol` is a `NonLinearSolution`; it contains the solution `u`, the residuals at each iteration, the number of iterations and whether the algorithm converged (use `BifurcationKit.converged(sol)`).

## Other flavours of `newton`

The Newton algorithm is also implemented for specific functionals, in which case a **simplified** entry point `newton` (exported) is provided which dispatches on the type of the problem:

- deflated problems: pass a [`DeflationOperator`](@ref) to `solve` in order to find *other* solutions than the known ones, see [Deflated problems](@ref);
- periodic orbits computed by finite differences, shooting or collocation, see the corresponding pages [Periodic orbits based on Trapezoidal rule](@ref), [Periodic orbits based on the shooting method](@ref), [Periodic orbits based on orthogonal collocation](@ref);
- travelling waves, see [Freezing problems, symmetries and waves](@ref);
- refining branch points of codimension 2, see for example `newton_hopf`, `newton_fold`, `newton_bt` and the pages [Fold / Hopf Continuation](@ref) and [Bogdanov-Takens refinement](@ref).

## Example

The tutorial [Temperature model](@ref temperature) presents the various jacobians (direct and iterative ones) in a concrete setting.
