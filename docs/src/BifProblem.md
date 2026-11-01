# Bifurcation problems

```@contents
Pages = ["BifProblem.md"]
Depth = 3
```

The idea behind `BifurcationKit` is to compute bifurcation diagrams in memory limited environments where the device can barely hold the current continuation state. You can thus disable the default saving of all solutions along the branch and all eigenvectors (see [`ContinuationPar`](@ref) to change this behavior). Still, one needs to save a few solution indicators, like for plotting. This is the reason for the function `record_from_solution` (see below).

## Generic bifurcation problem

[`BifurcationProblem`](@ref) is the basic / generic structure for encoding a bifurcation problem ; it holds the following fields:

- the vector field
- an initial guess
- a set of parameters
- a parameter axis

as well as user defined functions for 

- plotting, `plot_solution`
- recording (`record_from_solution`) indicators about the solution when this one is too large to be saved at every continuation step.
- update function `update!`

### Example

```julia
f(x,p) = @. sin(x * p.a)
u0 = zeros(100_000_000) 
params = (a = 1.0, b = 2.0)

# record a few components / indicators about x 
myRecord(x,p;k...) = (x1 = x[1], max = maximum(x), nrm = norm(x, Inf))

prob = BifurcationProblem(f, u0, params, (@optic _.a);
	record_from_solution = myRecord
	)
```


## Problem modification

In case you want to modify an existing problem, you should use the following method

```@docs
re_make(prob::BifurcationKit.AbstractBifurcationProblem;
		u0 = prob.u0,
		params = prob.params,
		lens = prob.lens,
		record_from_solution = prob.record_from_solution,
		plot_solution = prob.plot_solution,
		J = missing,
		Jᵗ = missing,
		d2F = missing,
		d3F = missing)
```

### Example

```@example 
using BifurcationKit
F(x,p) = @. p.a + x^2
# parameters
par = (a = 0., b = 2)
prob = BifurcationProblem(F, zeros(3), par, (@optic _.a))
# change u0
prob2 = BifurcationKit.re_make(prob, u0 = rand(3))
```

## Problem update

Some problems must be modified **in place during the continuation**, for instance when the
discretization is adaptive (mesh, moving Poincaré section) or when an expensive cache /
preconditioner must follow the solution. This is the purpose of the hook
`BifurcationKit.update!`.

It is called after each **successful** Newton correction (and eigenvalue computation), before
the step size is adapted and the next predictor is computed; the solution is recorded/saved
just after. It has the signature

```julia
BifurcationKit.update!(prob, iter, state) -> Bool
```

where `prob` is the current problem, `iter::ContIterable` and `state::ContState` describe the
current continuation step (solution `getx(state)`, parameter `getp(state)`, step number, stable
manifold, ...). Returning `false` stops the continuation. The default implementation
(`BifurcationKit.update_default`) returns `true`.

For a custom `AbstractBifurcationProblem`, implement the method above. Note that
`update!` is called on **converged** states, including bisection steps: guard with
`BifurcationKit.in_bisection(state)` if your update modifies the discretization.

### Built-in problems

For [`BifurcationProblem`](@ref), [`ODEBifProblem`](@ref) and [`PDEBifProblem`](@ref), pass the
function with the keyword `update!`. It is stored in the problem and called as
`myupdate(iter, state)`; the current problem is available as `BifurcationKit.getprob(iter)`:

```julia
function myupdate(iter, state)
    prob = BifurcationKit.getprob(iter)      # current problem
    x = getx(state)                          # current solution
    p = getparams(iter, state)               # full parameter set
    # mutate a cache / mesh / preconditioner stored in prob using (x, p)
    return true
end

prob = BifurcationProblem(F, u0, params, lens; update! = myupdate)
```

!!! note "`update!` vs `re_make` vs `finalise_solution`"
    - [`re_make`](@ref) rebuilds a problem (`u0`, `params`, `lens`, derivatives, ...) and returns a new object.
    - `update!` mutates the current problem in place and is called by the continuation.
    - `finalise_solution` is a post-step callback `(z, τ, step, contResult)`; it is convenient to post-process or to stop the continuation, see [`faq.md`](faq.md).

### Saving and restoring the problem state

When `update!` changes the problem (e.g. the mesh of a collocation method), the solutions
saved along the branch must remember that state, otherwise post-processing (normal forms,
codim 2 continuation, ...) uses the wrong discretization. This is the role of:

- `save_solution(prob, x, p)`: builds the object stored in `br.sol` (defaults to `x`);
- `saved_solution(saved)`: extracts the state vector from the saved object (identity by default);
- `restore_problem!(prob, saved, pars)`: restores the problem to the state it had when `saved`
  was recorded (no-op by default).

The collocation method for periodic orbits uses this mechanism to store the mesh together with
the orbit (`POSavedSolutionAndState`) and restore it later; see
[`periodicOrbitCollocation.md`](periodicOrbitCollocation.md) for an application to mesh
adaptation.

For the built-in problems, the restoration can be customised with the keyword
`restore_problem!`, called as `myrestore(prob, x_saved, pars)`:

```julia
prob = BifurcationProblem(F, u0, params, lens; restore_problem! = myrestore)
```

## Minimal interface of `AbstractBifurcationProblem`

[`BifurcationProblem`](@ref) / [`ODEBifProblem`](@ref) already implement the whole interface. A custom type
`struct MyPB <: BifurcationKit.AbstractBifurcationProblem …` must provide the methods below,
**depending on the functionalities** it will be used with.

### Core (always required)

- `residual(pb, x, p)` — the vector field $F(x,p)$
- `jacobian(pb, x, p)` — its derivative $J(x,p)=d_xF$
- `getparams(pb)`, `getlens(pb)` — parameters and continuation axis (or the fields `params`, `lens`)
- `getu0(pb)` — initial guess (defaults to the field `u0`)

### Required per functionality

| Functionality | Additional methods |
|---|---|
| Newton / 1-parameter continuation | `residual`, `jacobian` (+ `getlens` so `setparam` works) |
| Stability & detection (Hopf / Fold …) | `jacobian` returning a matrix / operator accepted by the eigen solver (see [`eigensolver.md`](eigensolver.md)); a Hopf pair requires `nev ≥ 2` and complex support |
| Adjoint / jvp / matrix-free paths | `dF(pb, x, p, dx)` (jvp, used by `apply_jacobian`), `jacobian_adjoint`, and the traits `has_adjoint`, `is_symmetric` (non-symmetric problems fall back to `transpose(J)` when `J` is an array) |
| Fold / Hopf refinement (MA) and codim-2 curves | bordered solves on `J` / `Jᵗ`, `getdelta` for finite differences, `d2F` (or automatic AD/FD when `usehessian = true`) |
| Normal forms (Hopf, Cusp, Bogdanov–Takens, …) | parameter derivative (`R01`) and the jet `d2F`, `d3F` (ForwardDiff is used when the residual is smooth) |
| [DAE](dae.md) (mass matrix) | `getmassmatrix` (and `is_mass_matrix_constant`) |
| Inplace / GPU (optional) | `residual!`, `jacobian!`, `isinplace` |
| Plotting / recording (optional) | `plot_solution`, `record_from_solution`, `save_solution`, `restore_problem!`, `update!` |

### Defaults and remarks

Most generic helpers assume the **conventional fields** `u0`, `params`, `lens`,
`recordFromSolution`, `plotSolution` (and `VF` for the jet), see
`getu0`, `getparams`, `setparam`, `re_make`, `apply_jacobian` in the source.
If your type stores these data differently, override the corresponding accessors,
as well as `re_make` (used by the Fold/Hopf MA formulations and codim-2
continuation).

There is **no generic fallback** for `getdelta`, `is_symmetric`, `has_adjoint`,
`isinplace`, `save_solution`, `update!` and the jet methods on an arbitrary
`AbstractBifurcationProblem`. They are provided for the built-in problems built
on a `BifFunction` (subtypes of `AbstractAllJetBifProblem`) and for
[`BifurcationProblem`](@ref) / [`ODEBifProblem`](@ref). For a custom type,
either implement them or wrap your vector field in a `BifFunction` / subtype
`AbstractAllJetBifProblem` so that the jet (`dF`, `d2F`, `d3F`, …) and the
associated traits are obtained automatically.

When using [`BifurcationProblem`](@ref), the functions `J`, `Jᵗ`, `d2F`,
`d3F` are optional: if not provided, `BifFunction` differentiates `residual`
with **ForwardDiff**, so `F` must be smooth and compatible with `ForwardDiff`
(for large-scale problems, provide an analytic `J`).

## 🚧🚧 Automatic option setting (work in progress) 🚧🚧

Setting the continuation options can be difficult for new comers. In the case of small ODE, we suggest to use `ODEBifProblem` instead of `BifurcationProblem`. Indeed, in this case, most continuation options will be set up automatically so that very good performance is achieved.

In the case of `BifurcationProblem`, the user has to set up these options manually. This can be useful for large scale problems where the specificities of the problem have to be used.