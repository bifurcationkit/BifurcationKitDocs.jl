# BVP based on Trapezoidal rule

```@contents
Pages = ["bvpTrapeze.md"]
Depth = 3
```

!!! warning "Work in progress"
    The BVP module is under active development. We isolate features in the BVP module. Hence `BVP.Trapeze` is the analog of `Trapeze` (for periodic orbits) but we have not yet performed merge of the two implementations. The API may change in future releases.

The Trapezoid method discretizes the BVP using finite differences based on a trapezoidal rule. The discretization is implemented in the structure `BVP.Trapeze` which wraps the periodic orbit `Trapeze` solver.

## Mathematical formulation

We look for solutions $(x(t_0), \dots, x(t_f))$ of

$$u'(t) = F(u(t), p),\quad u\in\mathbb R^n,\ t\in[t_0,t_f]$$

subject to boundary conditions

$$g(u(t_0), u(t_f), p)=0.$$

After discretization on $M$ time slices with mesh $\{0=t_1<\dots<t_M=1\}$ normalized to $[0,1]$, the trapezoidal rule gives:

$$\begin{array}{l}
u_{i+1} - u_i = \frac{h_i}{2}\big(F(u_i) + F(u_{i+1})\big),\quad i=1,\dots,M-1\\
g(u_1, u_M, p) = 0
\end{array}$$

where $h_i = \delta_T\cdot w_i$ with $\delta_T = t_f-t_0$ and $w_i$ are the normalized mesh weights.

## Usage

```julia
using BifurcationKit
const BK = BifurcationKit

F(u, p) = [u[2], -p.ω^2 * u[1]]
g(u0, u1, p) = [u0[1], u1[1]]

model = BK.BVP.BVPModel(F, g; n=2)
disc = BK.BVP.Trapeze(M=100)
bvp = BK.BVP.discretize(model, disc)

x0 = BK.BVP.generate_solution(bvp, t -> [cos(t), sin(t)])
prob = BK.BVP.BVPBifProblem(bvp, x0, (ω=1.0,), (@optic _.ω))

br = continuation(prob, PALC(), ContinuationPar())
```

## Non-uniform mesh

The mesh can be non-uniform by passing a vector of step sizes:

```julia
mesh = BK.TimeMesh([0.2, 0.3, 0.5])
disc = BK.BVP.Trapeze(M=4, mesh=mesh)
```

## Jacobian

The Jacobian is computed analytically using a block-banded structure. The default uses `AutoDiffDense` for the vector field Jacobian. You can also use `BK.Dense` for a dense analytical Jacobian:

```julia
J = BK.BVP.bvp_jacobian(bvp, BK.Dense(), x0, p)
```

## References

[^Uecker]: Uecker, Hannes. Hopf Bifurcation and Time Periodic Orbits with Pde2path – Algorithms and Applications. Communications in Computational Physics 25, no. 3 (2019)

[^Lust]: Lust, Kurt, Numerical Bifurcation Analysis of Periodic Solutions of Partial Differential Equations, PhD thesis, 1997.
