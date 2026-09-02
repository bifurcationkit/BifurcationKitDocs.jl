# Eigen solvers (Eig)

```@contents
Pages = ["eigensolver.md"]
Depth = 3
```

> The eigen solvers must be subtypes of `AbstractEigenSolver`.

They provide a way of computing the eigen elements of the Jacobian `J` (or, more generally, of the linear operator governing the linearized dynamics). These eigen elements are used throughout the package to

- detect bifurcation points (Fold, Branch point, Hopf) and assess the stability of the solutions found along a branch,
- feed the algorithms that refine or continue these bifurcation points,
- compute Floquet exponents for periodic orbits (see [Floquet solvers](@ref Library-FLOQUET)).

Every eigen solver is a callable object. Such an eigen solver `eig` will be called like `ev, evecs, success, itnumber = eig(J, nev; kwargs...)` throughout the package, `nev` being the number of requested eigen elements and `kwargs...` being used to send information about the algorithm (shift, `which` eigen-elements to extract, tolerances, ...). The returned quantities are

- `ev` the computed eigenvalues,
- `evecs` the associated eigenvectors,
- `success::Bool` a flag indicating whether the computation converged,
- `itnumber` the number of internal iterations (or matrix-vector products) that were performed; it equals `1` for direct solvers.

!!! note "Full spectrum vs. few eigen-elements"
    The direct solver [`DefaultEig`](@ref) computes the *full* spectrum while the iterative solvers only compute a small number of eigen-elements. In both cases only the `nev` first ones (after sorting, see the warnings below) are returned.

Here is an example of the simplest of them (see `src/EigSolver.jl` for the true implementation) to give you an idea:

```julia
struct DefaultEig <: AbstractEigenSolver end

function (l::DefaultEig)(J, nev; kwargs...)
	# I put Array so we can call it on small sparse matrices
	F = eigen(Array(J))
	I = sortperm(F.values, by = real, rev = true)
	nev2 = min(nev, length(I))
	return F.values[I[1:nev2]], F.vectors[:, I[1:nev2]], true, 1
end
```

!!! warning "Eigenvalues ordering"
    The eigenvalues must be ordered by decreasing real part for the detection of bifurcations to work properly. The eigen solvers shipped with `BifurcationKit` perform this sorting automatically (see the keyword `by` / `which` of each solver below). If you implement your own solver, you are responsible for it.

!!! warning "Eigenvectors"
    Several routines (for example [`newton_hopf`](@ref) and the normal form computations) need to extract the eigenvector number `i` from the returned eigenvectors. The default implementation

    ```julia
    geteigenvector(eigsolver, eigenvectors, i::Int) = eigenvectors[:, i]
    ```

    assumes that the eigenvectors are stored as the columns of a matrix. If your eigen solver stores them differently (this is the case of `EigKrylovKit` which returns a vector of vectors), you have to implement the method `geteigenvector(eigsolver, eigenvectors, i::Int)`.

## Methods for computing eigenvalues

Like for the linear solvers, computing the spectrum of operators $A$ associated to PDEs is a highly non trivial task because of the clustering of eigenvalues. Most methods are based on the so-called [power method](https://en.wikipedia.org/wiki/Power_iteration) but this only yields the eigenvalues with largest modulus. In the case of the Laplacian operator, this can be disastrous and it is better to apply the power method to $(\sigma I-A)^{-1}$ instead.

This method, called **Shift-invert**, is built-in for the solvers [`EigArpack`](@ref) and [`EigArnoldiMethod`](@ref), see below. It is mostly used to compute interior eigenvalues, i.e. the eigenvalues lying close to a given shift $\sigma$. For the solver [`EigKrylovKit`](@ref), which does not accept a shift directly, one can use the generic wrapper [`ShiftInvert`](@ref): it builds the operator $(\sigma I - A)^{-1}$ with the help of a linear solver (for example `GMRESKrylovKit`) and applies an eigen solver to it.

In some cases, it may be advantageous to consider the **Cayley transform** $(\sigma I-A)^{-1}(\tau I+A)$ to focus on a specific part of the spectrum. As it is mathematically equivalent to the Shift-invert method, we did not implement it.

## List of implemented eigen solvers

The detailed information for each one of them is located in the [API](@ref Library-EIG).

1. **`DefaultEig`** — direct solver based on the Julia function `eigen` for **matrices**. You can create it via `eig = DefaultEig()`. Because it computes the full spectrum of a dense matrix, it should only be used for moderate dimensions (the matrix is densified internally, so small sparse matrices are also accepted). You can specify how the eigenvalues are ordered through the keyword `which`, a function of an eigenvalue; by default `which = real`, so the eigen-elements are returned **by decreasing real part** (e.g. `DefaultEig(which = abs)` returns them by decreasing modulus). You can then compute the 3 eigen-elements of `J` of largest real part like `eig(J, 3)`. This solver does **not** accept Matrix-Free (functional) Jacobians.

2. **`EigArpack`** — iterative solver based on [Arpack.jl](https://github.com/JuliaLinearAlgebra/Arpack.jl). You can create one via `eigsolver = EigArpack()` and pass appropriate options (see [Arpack.jl](https://github.com/JuliaLinearAlgebra/Arpack.jl)). The first two (optional, positional) arguments are the shift `sigma` and the selection criterion `which`, *e.g.* `EigArpack(σ, :LM)`, the remaining options (`tol`, `maxiter`, `ritzvec`, `v0`, ...) being passed as keyword arguments and forwarded to `Arpack.eigs`. By default (`which = :LR`) it computes the eigen-elements of largest real part; it then re-orders them by decreasing real part (the sorting function can be changed with the keyword `by`). This solver works for (sparse) matrices as well as Matrix-Free Jacobians `dx -> J(dx)`. In the latter case you need to tell the eigensolver the dimension of the state space by giving an example of vector: `eig = EigArpack(v0 = zeros(10))`; you can then compute 3 eigen-elements using `eig(dx -> J(dx), 3)`.

    !!! note "Shift-Invert with `EigArpack`"
        You can compute the eigen-elements close to a shift $\sigma$ by passing it as the first argument: `EigArpack(σ, :LM)`. Arpack then applies a Shift-Invert strategy internally, and one looks for the eigen-elements of largest magnitude of the shifted operator (`which = :LM`) so that the eigen-elements closest to $\sigma$ are returned.

3. **`EigKrylovKit`** — iterative solver based on `KrylovKit.jl`. You create one via `eig = EigKrylovKit()` and pass appropriate options (see [KrylovKit.jl](https://github.com/Jutho/KrylovKit.jl)): the Krylov dimension `dim`, the tolerances `tol`, the maximum number of iterations `maxiter`, the verbosity `verbose`, ... The eigen-elements to extract are selected through the keyword `which` (e.g. `:LR` largest real part — the default, `:LM` largest modulus, `:SR` smallest real part, ...). Contrarily to `EigArpack`/`EigArnoldiMethod`, the eigen-elements are **not** re-ordered internally: keep `which = :LR` for a correct bifurcation detection. This solver works for (sparse) matrices as well as Matrix-Free Jacobians `dx -> J(dx)`. In the latter case you need to tell the eigensolver the dimension of the state space by giving an example of vector: `eig = EigKrylovKit(x₀ = zeros(10))`; you can then compute 3 eigen-elements using `eig(dx -> J(dx), 3)`. There is no built-in shift; use the wrapper [`ShiftInvert`](@ref) instead.

4. **`EigArnoldiMethod`** — iterative solver based on [ArnoldiMethod.jl](https://github.com/haampie/ArnoldiMethod.jl). You create one via `eig = EigArnoldiMethod()` and pass appropriate options (see [ArnoldiMethod.jl](https://github.com/haampie/ArnoldiMethod.jl)): the constructor reads `EigArnoldiMethod(; sigma = nothing, which = ArnoldiMethod.LR(), x₀ = nothing, kwargs...)` where the keyword arguments (`tol`, `mindim`, `maxdim`, `restarts`, ...) are forwarded to `ArnoldiMethod.partialschur`. The eigen-elements to extract are selected through the keyword `which`, which takes values like `LR()`, `LM()`, `SR()`, ... (the returned eigen-elements are re-ordered by decreasing real part, the sorting function being changed with the keyword `by`). This solver works for (sparse) matrices as well as Matrix-Free Jacobians `dx -> J(dx)`. In the latter case you need to tell the eigensolver the dimension of the state space by giving an example of vector: `eig = EigArnoldiMethod(x₀ = zeros(10))`; you can then compute 3 eigen-elements using `eig(dx -> J(dx), 3)`.

    !!! note "Shift-Invert with `EigArnoldiMethod`"
        You can compute the eigen-elements close to a shift $\sigma$ by passing the keyword `sigma = σ` and asking for the eigen-elements of largest magnitude of the shifted operator, e.g. `which = LM()`. In that case the matrix case is solved by factorizing $\sigma I - J$; note that this Shift-Invert strategy is not available for Matrix-Free Jacobians.

5. **`ShiftInvert`** — a general Shift-Invert wrapper which can be combined with *any* eigen solver, e.g. `ShiftInvert(sigma, ls, eig)` where `ls` is a linear solver used to apply $(\sigma I - J)^{-1}$ and `eig` the eigen solver applied to this operator. This is the recommended way to compute interior eigenvalues with `EigKrylovKit`. See [`ShiftInvert`](@ref).

!!! tip "Slow computations "
    This is probably due to iterative refinement conducted by `SuiteSparse` as explained in this blog [post](https://discourse.julialang.org/t/some-eigenpairs-from-a-large-sparse-nonsymmetric-matrix-julia-vs-matlab/93742). You can disable this using

    ```julia
    using SuiteSparse
    SuiteSparse.UMFPACK.umf_ctrl[8] = 0
    ```

## Generalized eigen solvers (GEV)

Some problems require solving the *generalized* eigen problem $A\,x = \lambda\, B\,x$. This is the case, for example, of DAE problems $M\,u' = f(u)$ for which the eigen-elements of the Jacobian are obtained by solving $J\,x = \lambda\, M\,x$, or of the Floquet problem based on a collocation method.

Associated to an eigensolver `eig` in `(DefaultEig, EigArnoldiMethod, EigArpack)`, a GEV is provided which can be called as

```julia
    gev(eig, A, B, nev; k...)
```

where `A` is the Jacobian-like operator, `B` the mass matrix/operator and `nev` the number of requested eigen-elements.

Another way is to rely on `EigenMassMatrix`, which wraps an eigen solver and a mass matrix `B` into a plain eigen solver `EigenMassMatrix(B, eig)`. This is convenient for DAE problems. See [`EigenMassMatrix`](@ref).
