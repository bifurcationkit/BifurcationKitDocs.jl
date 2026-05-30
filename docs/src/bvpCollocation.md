# BVP based on orthogonal collocation

```@contents
Pages = ["bvpCollocation.md"]
Depth = 3
```

!!! warning "Work in progress"
    The BVP module is under active development. We isolate features in the BVP module. Hence `BVP.Collocation` is the analog of `Collocation` (for periodic orbits) but we have not yet performed merge of the two implementations. The API may change in future releases.

The Collocation method discretizes the BVP using orthogonal collocation on `Ntst` mesh intervals with polynomials of degree `m`. This is the most precise discretization method and supports automatic mesh adaptation during continuation.

## Mathematical formulation

We look for solutions of

$$u'(t) = F(u(t), p),\quad u\in\mathbb R^n,\ t\in[t_0,t_f]$$

with boundary conditions $g(u(t_0), u(t_f), p)=0$.

The time domain is partitioned into $N_{tst}$ intervals. On each interval, the solution is approximated by a polynomial of degree $m$ and the residual is enforced at $m$ Gauss-Legendre collocation points.

## Usage

```julia
using BifurcationKit
const BK = BifurcationKit

F(u, p) = [u[2], -p.ω^2 * u[1]]
g(u0, u1, p) = [u0[1], u1[1]]

model = BK.BVP.BVPModel(F, g; n=2)
disc = BK.BVP.Collocation(Ntst=20, m=4)
bvp = BK.BVP.discretize(model, disc)

x0 = BK.BVP.generate_solution(bvp, t -> [cos(t), sin(t)])
prob = BK.BVP.BVPBifProblem(bvp, x0, (ω=1.0,), (@optic _.ω))

br = continuation(prob, PALC(), ContinuationPar())
```

## Mesh adaptation

The mesh can be automatically adapted during continuation to better resolve the solution:

```julia
disc = BK.BVP.Collocation(Ntst=20, m=4, meshadapt=true,
    K=100.0, update_every_step=1)
```

The adaptation works by equidistributing the collocation error across mesh intervals. It is triggered every `update_every_step` steps when the Newton solver has converged and the continuation step is not inside a bisection.

## Jacobians

Two analytical Jacobian types are available:

- `DenseAnalytical()`: dense analytical Jacobian, suitable for small to medium systems
- `FullSparse()`: sparse analytical Jacobian, suitable for large sparse systems

```julia
prob_ana = BK.BVP.BVPBifProblem(bvp, x0, params, (@optic _.a);
    jacobian = BK.DenseAnalytical())
```

## References

[^Dankowicz]: Dankowicz, Harry, and Frank Schilder. Recipes for Continuation. Computational Science & Engineering. Philadelphia, PA: Society for Industrial and Applied Mathematics, 2013.

[^Doedel]: Doedel, Eusebius J. “Lecture Notes on Numerical Analysis of Nonlinear Equations.” In Numerical Continuation Methods for Dynamical Systems, edited by Bernd Krauskopf, Hinke M. Osinga, and Jorge Galán-Vioque, 1–49. Understanding Complex Systems. Dordrecht: Springer Netherlands, 2007.
