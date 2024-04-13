# [Deflated problems](@id Deflated-problems)

!!! unknown "References"
    P. E. Farrell, A. Birkisson, and S. W. Funke. **Deflation techniques for finding distinct solutions of nonlinear partial differential equations**. SIAM J. Sci. Comput., 2015.,

Assume you want to solve $F(x)=0$ with a Newton algorithm but you want to avoid the algorithm to return some already known solutions $x_i,\ i=1\cdots n$.

The idea proposed in the paper quoted above is to penalize these solutions by looking for the zeros of the function $G(x):={F(x)}{M(x)}$ where

$$M(x) = \prod_{i=1}^n\left(\|x - x_i\|^{-2p} + \alpha\right)$$

and $\alpha>0$. Obviously $F$ and $G$ have the same zeros away from the $x_i$s but the factor $M$ penalizes the residual of the Newton iterations of $G$, effectively producing zeros of $F$ different from $x_i$.

!!! tip "Tip"
    In some case, you may want to use a custom distance, in place of the squared norm $$\|\cdot\|^2$$. Please see [`DeflationOperator`](@ref) for how to do this.

## Encoding of the functional

A composite type [`DeflationOperator`](@ref) implements this functional. Given a deflation operator `M = DeflationOperator(p, dot, α, xis)`, you can build a deflated functional `pb = DeflatedProblem(F, J, M)` which you can use to access the values of $G$ by doing `pb(x)`. A Matrix-Free / Sparse linear solver is implemented which works on the GPU.

> the `dot` argument in `DeflationOperator` lets you specify a dot product from which the norm is derived in the expression of $M$.

See example [Snaking computed with deflation](@ref sh2dfd).

Note that you can add new solution `x0` to `M` by doing `push!(M, x0)`. Also `M[i]` returns `xi`.

## Computation with `newton`

Most newton functions can be used with a deflated problem, see for example [Snaking computed with deflation](@ref sh2dfd). The idea is to pass the deflation operator `M`. For example, we have the following overloaded method, which works on GPUs:

```julia
newton(prob::BifurcationKit.AbstractBifurcationProblem,
		defOp::DeflationOperator,
		options::NewtonPar,
		_linsolver = DefProbCustomLinearSolver();
		kwargs...)
```

We refer to the regular [`newton`](@ref) for more information. This newton penalises the roots saved in `defOp.roots`. 

Compared to [`newton`](@ref), the only different arguments are

- `defOp::DeflationOperator` deflation operator
- `linsolver` linear solver used to invert the Jacobian of the deflated functional.
    - custom solver `DefProbCustomLinearSolver()` with requires solving two linear systems `J⋅x = rhs`.
    - For other linear solvers `<: AbstractLinearSolver`, a matrix free method is used for the deflated functional.
    - if passed `Val(:autodiff)`, then `ForwardDiff.jl` is used to compute the jacobian of the deflated problem
    - if passed `Val(:fullIterative)`, then a full matrix free method is used.


## Simple example

In this basic example, we show how to get the different roots of `F`

```@example DEFNEWTON
using BifurcationKit, LinearAlgebra
F(x, p) = @. (x-1) * (x-2)
# define a deflation operator which deflates the 
# already know solution x = 1
deflationOp = DeflationOperator(2, dot, 0.1, [ones(1)])
# define a problem, this compute jacobian automatically
prob = BifurcationProblem(F, zeros(1), nothing)
# call deflated newton
sol = newton(prob, deflationOp, NewtonPar())
if BifurcationKit.converged(sol)
    println("We found the additional root: ", sol.u)
else
    println("Deflated newton did not converge!")
end
```

!!! tip "Tip"
    You can use this method for periodic orbits as well by passing the deflation operator `M` to the newton method
