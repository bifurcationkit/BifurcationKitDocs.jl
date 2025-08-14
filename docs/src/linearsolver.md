# Linear solvers (LS)

> If you provide your own linear solver, it must be a subtype of `AbstractLinearSolver` otherwise `BifurcationKit.jl` will not recognize it. See example just below. 

The linear solvers provide a way of inverting the Jacobian `J` or solving `J * x = rhs`. Such linear solver `linsolve` will be called like `sol, success, itnumber = linsolve(J, rhs; kwargs...)` throughout the package.

Here is an example of the simplest one (see `src/LinearSolver.jl` for the true implementation) to give you an idea, the backslash operator:

```julia
struct DefaultLS <: AbstractLinearSolver end

function (l::DefaultLS)(J, rhs; k...)
	return J \ rhs, true, 1
end
```

Note that for `newton` to work, the linear solver must return 3 arguments. The first one is the result, the second one is whether the computation was successful and the third is the number of iterations required to perform the computation.

You can then call it as follows (and it will be called like this in [`newton`](@ref))

```julia
ls = DefaultLS()
J = rand(2, 2) # example of linear operator
ls(J, rand(2))
```

## List of implemented linear solvers

The detailed information for each one of them is located in the [API](@ref Library-LS).

1. Default `\` solver based on `LU` or `Cholesky` depending on the type of the Jacobian. This works for sparse matrices as well. You can create one via `linsolver = DefaultLS()`.
2. GMRES from `IterativeSolvers.jl`. You can create one via `linsolver = GMRESIterativeSolvers()` and pass appropriate options.
3. GMRES from `KrylovKit.jl`. You can create one via `linsolver = GMRESKrylovKit()` and pass appropriate options.
4. All solvers in `Krylov.jl`. You can create one via `linsolver = KrylovLS()` or `KrylovLSInplace()` and pass appropriate options. This is available from `BifurcationKit@0.4.5`
    
!!! tip "Different linear solvers"
    By tuning the options of [`GMRESKrylovKit`](@ref), you can select CG, GMRES... see [KrylovKit.jl](https://jutho.github.io/KrylovKit.jl/stable/man/linear/#KrylovKit.linsolve).
    
!!! note "Other solvers"
    For `KrylovLS`, just pass the required solver, like `cg`.

    It is very straightforward to implement the Conjugate Gradients from [IterativeSolvers.jl](https://juliamath.github.io/IterativeSolvers.jl/dev/linear_systems/cg/) by copying the interface done for `gmres`. Same goes for `minres`,... Not needing them, I did not implement this.

## Preconditioner

Preconditioners should be considered when using Matrix Free methods such as GMRES. `GMRESIterativeSolvers` provides a very simple interface for using them. For `GMRESKrylovKit`, we implemented a left preconditioner. Note that, for `GMRESKrylovKit`, you are not restricted to use `Vector`s anymore. 

A curated list is provided by the package [`LinearSolve.jl`](https://docs.sciml.ai/LinearSolve/stable/basics/Preconditioners/#Curated-List-of-Pre-Defined-Preconditioners).

Finally, here are some packages to use preconditioners:

1. [`Krylov.jl`](https://jso.dev/Krylov.jl/stable/preconditioners/#Packages-that-provide-preconditioners) provides some preconditioners.
2. [IncompleteLU.jl](https://github.com/haampie/IncompleteLU.jl) an ILU like preconditioner
3. [AlgebraicMultigrid.jl](https://github.com/JuliaLinearAlgebra/AlgebraicMultigrid.jl) Algebraic Multigrid (AMG) preconditioners. This works especially well for symmetric positive definite matrices.
4. [Preconditioners.jl](https://github.com/mohamed82008/Preconditioners.jl) A convenient interface to conveniently called most of the above preconditioners using a single syntax.
5. We provide a preconditioner based on deflation of eigenvalues (also called preconditioner based on Leading Invariant Subspaces) using a partial Schur decomposition. There are two ways to define one *i.e.* [`PrecPartialSchurKrylovKit`](@ref) and [`PrecPartialSchurArnoldiMethod`](@ref). 

!!! tip "Using Preconditioners"
    Apart from setting a preconditioner for a linear solver, it can be advantageous to change the preconditioner during computations, *e.g.* during a call to `continuation` or `newton`. This can be achieved by taking advantage of the callbacks to these methods. See the example [2d Ginzburg-Landau equation (finite differences, codim 2, Hopf aBS)](@ref cgl).