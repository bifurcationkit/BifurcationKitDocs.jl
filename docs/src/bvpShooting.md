# BVP based on Shooting

```@contents
Pages = ["bvpShooting.md"]
Depth = 3
```

!!! warning "Work in progress"
    The BVP module is under active development. We isolate features in the BVP module. Hence `BVP.Shooting` is the analog of `Shooting` (for periodic orbits) but we have not yet performed merge of the two implementations. The API may change in future releases.

The Shooting method solves the BVP by integrating the ODE between shooting points and imposing matching conditions at the boundaries.

## Mathematical formulation

We look for solutions of

$$u'(t) = F(u(t), p),\quad u\in\mathbb R^n,\ t\in[t_0,t_f]$$

with boundary conditions $g(u(t_0), u(t_f), p)=0$.

The time domain is divided into $M$ intervals. On each interval, the ODE is integrated from the shooting point $u_i$ using the flow:

$$\varphi_i = \text{evolve}(F, u_i, p, \Delta t_i)$$

The residual enforces matching between intervals and the boundary condition:

$$\begin{array}{l}
\varphi_i - u_{i+1} = 0,\quad i=1,\dots,M-1\\
g(u_1, \varphi_M, p) = 0
\end{array}$$

## Usage

```julia
using BifurcationKit, OrdinaryDiffEq
const BK = BifurcationKit

F(u, p) = [u[2], -p.ω^2 * u[1]]
g(u0, u1, p) = [u0[1], u1[1]]

model = BK.BVP.BVPModel(F, g; n=2)

# Shooting requires an ODE algorithm from OrdinaryDiffEq
disc = BK.BVP.Shooting(M=4, alg=Tsit5())
bvp = BK.BVP.discretize(model, disc)

x0 = BK.BVP.generate_solution(bvp, t -> [cos(t), sin(t)])
prob = BK.BVP.BVPBifProblem(bvp, x0, (ω=1.0,), (@optic _.ω))

br = continuation(prob, PALC(), ContinuationPar())
```

## Simple vs Multiple shooting

- **Simple shooting** ($M=1$): single integration from $t_0$ to $t_f$, cheaper but less stable
- **Multiple shooting** ($M>1$): multiple integration intervals, more stable for stiff problems

## Parallel shooting

Multiple shooting can be parallelized across the shooting intervals:

```julia
disc = BK.BVP.Shooting(M=8, alg=Tsit5(), parallel=true)
```

## Jacobian

The analytical Jacobian for shooting is computed using monodromy matrices from the variational equations.


