# Boundary Value Problems

```@contents
Pages = ["bvp.md"]
Depth = 3
```

!!! warning "Work in progress"
    The BVP module is under active development. We isolate features in the BVP module. Hence `BVP.Trapeze` is the analog of `Trapeze` (for periodic orbits) but we have not yet performed merge of the two implementations. The API may change in future releases.

Consider the Boundary Value Problem (BVP)

$$\frac{du}{dt}=F(u,p),\quad t\in[t_0,t_f]$$

with boundary conditions

$$g(u(t_0), u(t_f), p)=0.$$

A periodic orbit is a special case where $g(u_0,u_1)=u_0-u_1$ and $t_f-t_0=T$ is the period.

The BVP module provides a unified interface for solving such problems using the same discretization methods as periodic orbits:

- `BVPModel`: mathematical formulation (vector field + boundary conditions)
- `Trapeze`, `Collocation`, `Shooting`: discretization methods
- `discretize`: combine model + method into a `DiscretizedBVP`
- `BVPBifProblem`: wrap for continuation

## Basic Usage

```julia
using BifurcationKit
const BK = BifurcationKit

# 1. Define the vector field
F(u, p) = [u[2], -p.ω^2 * u[1]]

# 2. Define boundary conditions
g(u0, u1, p) = [u0[1], u1[1]]

# 3. Create a BVP model
model = BK.BVP.BVPModel(F, g; n=2)

# 4. Choose discretization
disc = BK.BVP.Trapeze(M=100)

# 5. Discretize
bvp = BK.BVP.discretize(model, disc)

# 6. Generate initial guess
x0 = BK.BVP.generate_solution(bvp, t -> [cos(t), sin(t)])

# 7. Create BVP bifurcation problem
prob = BK.BVP.BVPBifProblem(bvp, x0, (ω=1.0,), (@optic _.ω))

# 8. Continue
br = BK.continuation(prob, PALC(), ContinuationPar())
```

## Important notes

### 1. Relation with periodic orbits

The BVP module is a generalization of the periodic orbit solvers. It uses the same underlying discretizations (`Trapeze`, `Collocation`, `Shooting`) but exposes a cleaner interface where the boundary conditions are user-defined rather than hardcoded to $u(0)=u(T)$.

### 2. Accessing the solution

After a computation, you can obtain the time-domain solution using `get_solution_bvp`:

```julia
sol = BK.BVP.get_solution_bvp(bvp, x, p)
# sol.t contains the time mesh
# sol.u contains the solution at each time slice
```

For periodic orbits, you can also use `get_periodic_orbit`.

### 3. Record from solution

You can pass your own function to `continuation`. It is called as `record_from_solution(x, opt; k...)` where `opt.p` is the current continuation parameter.

### 4. Plot solution

Similarly, the method is called like `plot_solution(x, opt; k...)` where `opt.p` is the current value of the continuation parameter.

## Methods

We provide three discretization methods:

### Trapeze method

The Trapezoid method discretizes the BVP using a finite difference scheme based on the trapezoidal rule. It is usually fast and works well for small to medium sized problems. See [BVP Trapeze](bvpTrapeze.md).

### Collocation method

The Collocation method uses orthogonal collocation on time intervals. It is the most precise method and supports mesh adaptation during continuation. See [BVP Collocation](bvpCollocation.md).

### Shooting method

The Shooting method solves the BVP by integrating the ODE between shooting points and matching conditions. It is well suited for large-scale problems and supports parallel integration. See [BVP Shooting](bvpShooting.md).
