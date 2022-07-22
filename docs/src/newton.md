# Krylov-Newton algorithm

`BifurcationKit` is built upon the newton algorithm for solving (large-dimensional) nonlinear equations 

$$F(x)=0\in\mathbb R^n,\quad x\in\mathbb R^n.$$

Writing $J(x)\in\mathcal L(\mathbb R^n)$ the jacobian, the algorithm reads

$$x_{n+1} = x_n - J(x_n)^{-1}F(x_n)$$

with initial guess $x_0$.

The crux of the algorithm is to solve the linear system in $y$:

$$J(x_n)\cdot y = F(x_n).$$

To this end, we never form $J^{-1}$ like with `pinv(J)` but solve the linear system directly. 

## Space of solutions

For the algorithm to be defined, a certain number of operations on `x` need to be available. If you pass `x::AbstractArray`, you should not have any problem. Otherwise, your `x` must comply with the requirements listed in [Requested methods for Custom State](@ref).

## Different Jacobians

There are basically two ways to specify the jacobian:

1. Matrix based
2. Matrix-free.

In case you pass a matrix (in effect an `AbstractMatrix` like a sparse one,...), you can use the default linear solver from `LinearAlgebra` termed the backslash operator `\`. This is a **direct** method. This is the case 1 above.

Another possibility is to pass a function `J(dx)` and to use **iterative** linear solvers. In this case, this is termed a **Krylov-Newton** method. This is the case 2 above. In comparison to the Matrix-based case, there is no restriction to the number of unknowns $n$.

> The available linear solvers are explained in the section [Linear solvers (LS)](@ref).

One can find a full description of the Krylov-Newton method in the [API](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/library/#Newton). 

## Simple example

Here is a quick example to show how the basics work. In particular, the problem generates a matrix based jacobian using automatic differentiation.

```@example NEWTON
using BifurcationKit
F(x, p) = x.^3 .- 1
x0 = rand(10)
prob = BifurcationProblem(F, x0, nothing)
sol = newton(prob, NewtonPar(verbose = true))
```

## Example

The (basic) tutorial [Temperature model (Simplest example)](@ref) presents all cases (direct and iterative ones).

