# Linear solvers (LS)

```@contents
Pages = ["linearsolver.md"]
Depth = 3
```

> If you provide your own linear solver, it must be a subtype of `AbstractLinearSolver` otherwise `BifurcationKit.jl` will not recognize it. See example just below.

The linear solvers provide a way of inverting the Jacobian `J`, i.e. of solving $J\, x = \mathrm{rhs}$. Throughout the package, such a linear solver `linsolve` will be called like

```julia
sol, success, itnumber = linsolve(J, rhs; kwargs...)
```

The three returned quantities are the solution `sol`, a boolean `success` telling whether the solve succeeded (resp. converged for iterative solvers) and `itnumber` the number of iterations required by the computation (equal to `1` for a direct solver).

Here is an example of the simplest one (see `src/LinearSolver.jl` for the true implementation) to give you an idea, the backslash operator:

```julia
struct DefaultLS <: AbstractLinearSolver end

function (l::DefaultLS)(J, rhs; k...)
	return J \ rhs, true, 1
end
```

You can then call it as follows (and it will be called like this in [`newton`](@ref)):

```julia
ls = DefaultLS()
J = rand(2, 2) # example of linear operator
ls(J, rand(2))
```

!!! note "Several right hand sides"
    All linear solvers also accept **two** right hand sides: `linsolve(J, rhs1, rhs2; kwargs...)` returns `(sol1, sol2, success, itnumber)` with `itnumber` a tuple. This interface is used internally, e.g. by the bordered solvers of [Bordered linear solvers (BLS)](@ref), and it is where a direct solver can **reuse a factorization** for the second solve.

!!! note "Shifted systems"
    The linear solvers can also be called with the keyword arguments `a₀` and `a₁` to solve the shifted system $(a_0 I + a_1 J)\, x = \mathrm{rhs}$. Their defaults are `a₀ = 0` and `a₁ = 1`, i.e. one recovers $J\,x=\mathrm{rhs}$. This interface is used internally by the algorithms based on the Hopf problem, on Floquet multipliers or by the bordered solvers. If you provide your own linear solver, you can ignore these keywords as long as you do not need these features.

## List of implemented linear solvers

The detailed information for each one of them is located in the [API](@ref Library-LS).

1. **`DefaultLS`** — direct solver based on the default `\` operator, which dispatches to `LU`, `Cholesky` or a sparse factorization depending on the type of the Jacobian. It works for dense and sparse matrices. The option `useFactorization = true` (the default) caches the factorization so that a second linear solve with the same operator is very cheap; set it to `false` when your operator does not support `factorize` (e.g. some operators coming from `ApproxFun.jl`, or matrices such as `StaticArrays.MMatrix`). You can create one via `linsolver = DefaultLS()`.

2. **`GMRESIterativeSolvers`** — iterative solver based on `gmres` from [IterativeSolvers.jl](https://docs.juliahub.com/General/IterativeSolvers/stable/). You can create one via `linsolver = GMRESIterativeSolvers()` and tune it with the keyword arguments `abstol`, `reltol`, `restart`, `maxiter`, `N`, `verbose`, ... It supports Matrix-Free Jacobians (with `N` the dimension of the state space) and left/right preconditioners (fields `Pl`, `Pr`). The struct is mutable so that you can change the preconditioners on the fly.

3. **`GMRESKrylovKit`** — iterative solver based on `linsolve` from [KrylovKit.jl](https://jutho.github.io/KrylovKit.jl/stable/man/linear/#KrylovKit.linsolve). You can create one via `linsolver = GMRESKrylovKit()` and tune it with the keyword arguments `dim`, `atol`, `rtol`, `maxiter`, `verbose`, ...

    !!! tip "Different linear solvers"
        By tuning the options of [`GMRESKrylovKit`](@ref), you can select CG, GMRES... Indeed, `KrylovKit.jl` dispatches on the options `issymmetric`, `ishermitian`, `isposdef`, see [KrylovKit.jl](https://jutho.github.io/KrylovKit.jl/stable/man/linear/#KrylovKit.linsolve). A left preconditioner can be set through the field `Pl`.

4. **`KrylovLS` / `KrylovLSInplace`** — interfaces to the many solvers of [Krylov.jl](https://jso.dev/Krylov.jl). You can create one via `linsolver = KrylovLS()` and pass the required Krylov method through the keyword `KrylovAlg = :gmres` (the default): you have access to `cg`, `cr`, `gmres`, `symmlq`, `minres`, `cg_lanczos`, `cg_lanczos_shift_seq`, ... The extra keyword arguments are forwarded to the chosen Krylov method, e.g. `KrylovLS(atol = 1e-11, rtol = 1e-8)`. Left/right preconditioners are passed through the keywords `Pl`, `Pr` (see the Preconditioner section below). `KrylovLSInplace` pre-allocates the Krylov space (you then have to provide the workspace dimensions, see the docstring) which is very useful on the GPU or when many linear solves are performed. Both structs are mutable so that you can modify the preconditioners.

    !!! note "Other solvers"
        Thanks to the interface to `Krylov.jl`, you do not need to implement the Conjugate Gradients from IterativeSolvers.jl, `minres`, ... yourself: just select the algorithm with the keyword `KrylovAlg` as explained above.

## Preconditioner

Preconditioners should be considered when using Matrix Free methods such as GMRES. `GMRESIterativeSolvers` provides a very simple interface for using them (fields `Pl`, `Pr`). For `GMRESKrylovKit`, we implemented a left preconditioner (field `Pl`). Note that, for `GMRESKrylovKit`, you are not restricted to use `Vector`s anymore.

A curated list is provided by the package [`LinearSolve.jl`](https://docs.sciml.ai/LinearSolve/stable/basics/Preconditioners/#Curated-List-of-Pre-Defined-Preconditioners).

Finally, here are some packages to use preconditioners:

1. [`Krylov.jl`](https://jso.dev/Krylov.jl/stable/preconditioners/#Packages-that-provide-preconditioners) provides some preconditioners.
2. [IncompleteLU.jl](https://github.com/haampie/IncompleteLU.jl) an ILU like preconditioner
3. [AlgebraicMultigrid.jl](https://github.com/JuliaLinearAlgebra/AlgebraicMultigrid.jl) Algebraic Multigrid (AMG) preconditioners. This works especially well for symmetric positive definite matrices.
4. [Preconditioners.jl](https://github.com/mohamed82008/Preconditioners.jl) A convenient interface to conveniently called most of the above preconditioners using a single syntax.
5. We provide a preconditioner based on deflation of eigenvalues (also called preconditioner based on Leading Invariant Subspaces) using a partial Schur decomposition. There are two ways to define one *i.e.* [`PrecPartialSchurKrylovKit`](@ref) and [`PrecPartialSchurArnoldiMethod`](@ref).

!!! tip "Using Preconditioners"
    Apart from setting a preconditioner for a linear solver, it can be advantageous to change the preconditioner during computations, *e.g.* during a call to `continuation` or `newton`. This can be achieved by taking advantage of the callbacks to these methods. See the example [2d Ginzburg-Landau equation (finite differences, codim 2, Hopf aBS)](@ref cgl).
