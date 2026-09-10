# Differential-Algebraic Equations (DAE)

A (parametrised) *differential-algebraic equation* (DAE)
is written as

$$M(x, p)\cdot \dot x = F(x, p),$$

where $x \in \mathbb R^{n}$ is the state, $p$ the parameters and
$M(x, p)$ is the **mass matrix** (possibly singular for index-1 DAEs).
Denoting by $x_{eq}$ an equilibrium, i.e. a solution of the algebraic system
$F(x_{eq}, p) = 0$, the dynamics around $x_{eq}$ is governed by the
**generalized eigenproblem**

$$J(x_{eq}, p)\cdot v = \lambda\, M(x_{eq}, p)\cdot v,\qquad J \equiv d_xF,$$

so that stability is decided by the generalized eigenvalues $(\lambda)$.
More details on the generalized eigen solvers are given on the page
[`eigensolver.md`](eigensolver.md) (section *Generalized eigen problem*).

## Functionalities and support

The following table summarises which functionalities are supported depending
on the type of mass matrix. A *constant* mass matrix is independent of
$(x, p)$ (it can be singular, i.e. index-1); a *non-constant* one depends on
the state / parameters, $M = M(x,p)$.

| Functionality | Constant $M$ | Non-constant $M(x,p)$ |
|---|---|---|
| Equilibrium continuation & Newton | ✔ | ⚠ partial |
| Stability — generalized eigenvalues $(J, M)$ | ✔ | ⚠ partial |
| Fold / Hopf location & refinement (`newton_fold`, `newton_hopf`) | ✔ | ✘ |
| Hopf normal form ([`hopf_normal_form`](@ref)) | ✔ | ✘ |
| Codim 2 curves of Fold / Hopf + BT / Cusp detection | ✔ | ✘ |
| Detailed codim 2 normal forms (BT, Cusp, Bautin, ZH, HH) | ✘ | ✘ |
| Periodic orbits — Shooting (SciML `mass_matrix`) | ✔ | ⚠ |
| Periodic orbits — [`Trapeze`](@ref) + `massmatrix` | ✔ (dense Floquet; matrix-free in progress) | ✘ |

Here ✔ means supported, ⚠ means only partially supported (the mass matrix is
evaluated once at the current solution and its derivatives are not taken into
account), and ✘ means not supported (a clear error is raised).

## Wrapping a problem with a mass matrix

The equilibrium machinery of `BifurcationKit` acts on
[`BifurcationProblem`](@ref) (or [`ODEBifProblem`](@ref)) and treats the case
$M = I$ (identity). A problem with a general mass matrix is encoded with
[`BifurcationKit.DAEMassBifProblem`](@ref) which wraps an
[`ODEBifProblem`](@ref) together with the mass matrix $M$:

```julia
prob   = ODEBifProblem(F, u0, params, (@optic _.μ))
daeprob = BifurcationKit.DAEMassBifProblem(prob, Mass)
```

The mass matrix $M$ can be
- a constant matrix, e.g. `M = [1. 0; 0 0.]` (or a sparse matrix),
- a `UniformScaling` (e.g. `LA.I`, `α * LA.I`): it is then represented by the
  marker [`BifurcationKit.IdentityOperator`](@ref),
- a function $M(x, p)$: this case is **only partially supported** (basically only simple continuation).

The kind of mass matrix (constant vs. state dependent) is tracked by the
first type parameter, [`BifurcationKit.ConstantMass`](@ref) being the default
for constant matrices. Useful accessors are

| call | role |
|---|---|
| `BifurcationKit.getmassmatrix(pb, x, p)` | evaluate the mass matrix |
| `BifurcationKit.is_mass_matrix_constant(pb)` | is the mass matrix constant? |
| `re_make(pb; M = newMass, u0 = …)` | rebuild with another mass matrix |


## A small example

We consider the (polynomial) Stuart-Landau model with the invertible constant
mass matrix $M = \mathrm{diag}(2, 1)$:

```@example DAE1
using BifurcationKit, LinearAlgebra

# the vector field F(x, p)
function Fsl2(x, p)
    (; r, μ, ν, c3) = p
    u1, u2 = x
    ua = u1^2 + u2^2
    return [r * u1 - ν * u2 - ua * (c3 * u1 - μ * u2),
            r * u2 + ν * u1 - ua * (c3 * u2 + μ * u1)]
end

par_sl = (r = -0.1, μ = 0.132, ν = 1.0, c3 = 1.123)

# constant (diagonal) mass matrix
Mass = Diagonal([2.0, 1.0])

prob    = ODEBifProblem(Fsl2, zeros(2), par_sl, (@optic _.r))
daeprob = BifurcationKit.DAEMassBifProblem(prob, Mass)

# continuation with detection of the Hopf bifurcation
opts = ContinuationPar(ds = 0.01, dsmax = 0.02, dsmin = 1e-3,
    p_min = -0.3, p_max = 0.1)
br = continuation(daeprob, PALC(), opts; verbosity = 0, normC = norminf)

ind_h = findfirst(pt -> pt.type == :hopf, br.specialpoint)
println("Hopf point located at r = ", br.specialpoint[ind_h].param)

# Hopf normal form (mass matrix handled with the bordered vectors)
hp = BifurcationKit.hopf_normal_form(daeprob, br, ind_h; start_with_eigen = Val(false))
```

Note that the eigenvalue / stability computation during `continuation` is
performed on the pencil $(J, M)$, the eigen solver being automatically wrapped
into [`BifurcationKit.EigenDAE`](@ref) when the problem is a `DAEMassBifProblem`.

## Bifurcations points and normal forms

- **Fold / Hopf points** of the DAE are located and refined with the usual
  [`newton_fold`](@ref) / [`newton_hopf`](@ref) (Minimally Augmented
  formulation), the bordered systems now involving the mass matrix
  (e.g. $(J - i\omega M)$ for the Hopf problem).
- The **Hopf normal form** is available through
  [`get_normal_form`](@ref). Only `start_with_eigen = Val(false)` is
  supported for a mass matrix (the eigen-computation would otherwise ignore
  $M$): the right/left eigenvectors are obtained by solving bordered linear
  systems on the pencil $(J, M)$ and the left one is normalized by
  $\langle \zeta^\star, M \zeta\rangle = 1$.
- **Codimension 2**: the continuation of Fold / Hopf curves (2 parameters)
  and the *detection* of codim 2 points (BT, Cusp, Bautin, …) work for a
  constant mass matrix. The detailed **codim 2 normal forms** (e.g. for the
  Bogdanov-Takens point) are however **not implemented yet** for problems with
  a mass matrix and raise an error.

## Periodic orbits of DAEs

Periodic orbits of a DAE can be tackled in two ways:

1. **Shooting** the underlying SciML `ODEProblem` with a `mass_matrix` option,
   see the tutorial *Colpitts-type oscillator* (Tutorials → DAE examples).
2. **Finite differences** with the [`Trapeze`](@ref) discretization equipped
   with a `massmatrix`, see the page [`periodicOrbitTrapeze.md`](periodicOrbitTrapeze.md).

## Linear algebra internals

The shifted linear systems which appear with a mass matrix, e.g.
$(a_0 M + a_1 J)\,x = b$, are encoded by wrapping the pair $(M, J)$ in a
`BifurcationKit.MassAndJacobian` inside a `BifurcationKit.ShiftedOperator`.
Both direct and iterative / matrix-free (bordered) linear solvers take this
into account: see the pages [`linearsolver.md`](linearsolver.md) and
[`borderedlinearsolver.md`](borderedlinearsolver.md).

## API

```@docs
BifurcationKit.DAEMassBifProblem
BifurcationKit.ConstantMass
BifurcationKit.IdentityOperator
BifurcationKit.EigenDAE
```
