# Periodic orbits based on the shooting method

```@contents
Pages = ["periodicOrbitShooting.md"]
Depth = 3
```

A set of shooting algorithms is provided which are called either *Simple Shooting (SS)* if a single shooting is used and *Multiple Shooting (MS)* otherwise. For the exposition, we follow the PhD thesis[^Lust] and also [^Umbria].

We aim at finding periodic orbits for the Cauchy problem 

$$\tag{1} \frac{d x}{d t}=f(x)$$

and we write $\phi^t(x_0)$ the associated flow (or semigroup of solutions).

!!! tip "Tip about convenience functions"
    For convenience, we provide the following helper functions:
    - `get_periodic_orbit(prob, sol, par)` recomputes the periodic orbit from a shooting solution `sol` and returns it as a `BVPSolution`, which can then be plotted or post-processed;
    - `plot_periodic_shooting` plots a shooting solution (only available with the **Plots** backend);
    - during continuation, the default `record_from_solution` of a shooting problem records the maximum, the minimum, the amplitude and the period of the orbit.
    See the tutorials for examples of use.

!!! note "Cauchy problem with a mass matrix"
    The shooting method integrates the Cauchy problem (1) through an ODE solver. If the Cauchy problem involves a mass matrix, as in $M(x)\dot x = f(x)$, it must be encoded in the `ODEProblem`/`DAEProblem` passed to [`Shooting`](@ref) (or [`PoincareShooting`](@ref)) and is handled by the underlying ODE solver from `DifferentialEquations.jl`. There is no specific treatment of the mass matrix in the shooting discretization itself.

## Standard Shooting
### Simple shooting
A periodic orbit is found when we have a couple $(x, T)$ such that $\phi^T(x) = x$ and the trajectory is non constant. Therefore, we want to solve the equations $G(x,T)=0$ given by

$$\tag{SS}\begin{array}{l}{\phi^T(x)-x=0} \\ {s(x,T)=0}\end{array}.$$

The section $s(x,T)=0$ is a phase condition to remove the indeterminacy of the point on the limit cycle. Simple shooting corresponds to a **single** section; in practice, this is the case `M = 1` of the multiple shooting below.

### Multiple shooting
This case is similar to the previous one but more sections are used. To this end, we partition the unit interval with $m+1$ points $0=s_{0}<s_{1}<\cdots<s_{m-1}<s_{m}=1$
and consider the equations $G(x_1,\cdots,x_m,T)=0$ given by:

$$\begin{aligned}
\phi^{\delta s_1T}(x_{1})-x_{2} &=0 \\ 
\phi^{\delta s_2T}(x_{2})-x_{3} &=0 \\ & \vdots \\ 
\phi^{\delta s_{m-1}T}(x_{m-1})-x_{m} &=0 \\ 
\phi^{\delta s_mT}(x_{m})-x_{1} &=0 \\ s(x_{1}, T) &=0. 
\end{aligned}\tag{MS}$$

where $\delta s_i:=s_{i+1}-s_i$. The unknowns are the $m$ points $x_i$ together with the period $T$, so the system (MS) has $m\cdot N+1$ unknowns for $N$ the dimension of the state space. The Jacobian of the system of equations *w.r.t.* $(x,T)$ is given by 

$$\mathcal{J}=\left(\begin{array}{cc}{\mathcal J_c} & {\partial_TG} \\ {\star} & {d}\end{array}\right)$$

where the cyclic matrix $\mathcal J_c$ is

$$\mathcal J_c := 
\left(\begin{array}{ccccc}
{M_{1}} & {-I} & {} & {} \\ 
{} & {M_{2}} & {-I} & {}\\ 
{} & {} & {\ddots} & {-I}\\ 
{-I} & {} & {} & {M_{m}}\\ 
\end{array}\right)$$

and $M_i=\partial_x\phi^{\delta s_i T}(x_i)$.

### Section

The periodic orbits solutions of (SS) or (MS) are not uniquely defined because of the phase invariance. A section $s(x,T)=0$ (resp. $s(x_1,T)=0$) for (SS) (resp. (MS)) must be provided. The default is the same for both
$$ s(x,T) = T\cdot \langle x-x_\pi, \phi\rangle.$$

### Encoding of the functional

The functional is encoded in the composite type [`Shooting`](@ref). In particular, the user can pass its own time stepper or one can use the different ODE solvers in  [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl) which makes it very easy to choose a solver tailored for a specific problem. See the link [`Shooting`](@ref) for more information ;  for example on how to access the underlying functional, its jacobian...

### Jacobians

We provide many different jacobian types to take advantage of the formulations or the dimensionality. The jacobian is selected through the keyword `jacobian` in the constructor of [`Shooting`](@ref) (resp. [`PoincareShooting`](@ref)). For example, you can pass `jacobian = AutoDiffDense()` (the default for `Shooting`) or `jacobian = MatrixFree()` for a matrix-free implementation which is well suited to large dimensional problems. Note that all the internal linear solvers and jacobians are set up automatically so you don't need to do anything. See [`Shooting`](@ref) for the different `jacobian` available.

## Poincaré shooting

The algorithm is based on the one described in [^Sanchez] and [^Waugh].

We look for periodic orbit solutions of (1) using the hyperplanes

$$\Sigma_i = \{ x \mid \langle x - x_i^c, n_i \rangle = 0 \},\qquad i = 1, \dots, M,$$

which an initial periodic orbit guess intersects transversally.  Let $\Pi_i : \Sigma_i \to \Sigma_{i+1}$ (with the convention $\Sigma_{M+1}=\Sigma_1$) denote the partial Poincaré return map from the section $\Sigma_i$ to the next one $\Sigma_{i+1}$.  The key observation is that when each $x_i \in \mathbb{R}^N$ is constrained to lie in $\Sigma_i$, the problem becomes $(N-1)M$-dimensional rather than $NM$-dimensional.  Enforcing these hyperplane constraints explicitly is therefore necessary for the Newton iteration to converge reliably.

We thus need to parametrize these hyperplanes.

To this end, we introduce the projection operator $R_i:\mathbb R^N\to \mathbb R^{N-1}$ such that 

$$R_{i}\left(x_{1}, x_{2}, \ldots, x_{k_i-1}, x_{k_i}, x_{k_i+1}, \ldots, x_{N}\right)=\left(x_{1}, x_{2}, \ldots, x_{k_i-1}, x_{k_i+1}, \ldots, x_{N}\right)$$

where $k_i:=\arg\max_p |n_{i,p}|$. The inverse operator $E_i:\mathbb R^{N-1}\to\Sigma_i$ is defined by (where $\bar x:=R_i(x)$)

$$E_{i}(\bar x) := E_{i}\left(x_{1}, x_{2}, \ldots, x_{k_i-1}, x_{k_i+1}, \ldots, x_{N}\right)=
\left(x_{1}, x_{2}, \ldots, x_{k_i-1}, x^c_{i,k_i}-\frac{\bar{n}_i \cdot\left(\overline{x}-\overline{x}^c_{i}\right)}{n_{i,k_i}}, x_{k_i+1}, \ldots, x_{N}\right).$$ 

We note that $R_i\circ E_i = I_{\mathbb R^{N-1}}$ and that $E_i\circ R_i$ is the identity on the hyperplane $\Sigma_i$.

We then look for the fixed points $(\bar x_1,\cdots,\bar x_M)\in(\mathbb R^{N-1})^M$ of the composed maps:

$$\begin{aligned} 
\bar x_1 - R_1\Pi_M(E_M(\bar x_M))&=0 \\ 
\bar x_2 - R_2\Pi_1(E_1(\bar x_1))&=0 \\ & \vdots \\ 
\bar x_M - R_M\Pi_{M-1}(E_{M-1}(\bar x_{M-1}))&=0. 
\end{aligned}$$



### Encoding of the functional

The functional is encoded in the composite type [`PoincareShooting`](@ref). In particular, the user can pass their own time stepper or he can use the different ODE solvers in  [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl) which makes it very easy to choose a tailored solver: the partial Poincaré return maps are implemented using **callbacks**. See the link [`PoincareShooting`](@ref) for more information, in particular on how to access the underlying functional, its jacobian...

## Floquet multipliers computation


### Standard shooting
The Floquet multipliers are the eigenvalues of the monodromy matrix $M=M_M\cdots M_1$.

> Unlike the case with [Finite differences](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/periodicOrbitTrapeze/), the matrices $M_i$ are not sparse.

### Poincaré shooting
The (non trivial) Floquet exponents are eigenvalues of the Poincaré return map $\Pi:\Sigma_1\to\Sigma_1$. We have $\Pi = \Pi_M\circ\Pi_{M-1}\circ\cdots\circ\Pi_2\circ\Pi_1$. Its differential is thus

$$d\Pi(x)\cdot h = d\Pi_M(x_{M})d\Pi_{M-1}(x_{M-1})\cdots d\Pi_1(x_1)\cdot h$$

### Numerical method

We provide two methods to compute the Floquet coefficients.

- A **not very precise** algorithm for computing the Floquet multipliers is provided. The method, dubbed Quick and Dirty (QaD), is not numerically very precise for large / small Floquet exponents. 
It amounts to computing the eigenvalues of $M=M_M\cdots M_1$ (resp. $d\Pi$) for the Standard (resp. Poincaré) Shooting.
The method allows, nevertheless, to detect bifurcations of periodic orbits. It seems to work reasonably well for the tutorials considered here. For more information, have a look at [`FloquetQaD`](@ref).
- The state of the art method is based on a Periodic Schur decomposition. It is available through the package [PeriodicSchurBifurcationKit.jl](https://github.com/bifurcationkit/PeriodicSchurBifurcationKit.jl). For more information, have a look at `FloquetPQZ`.

## Computation with `newton`

We provide a simplified call to `newton` to locate the periodic orbit. Have a look at the tutorial [Continuation of periodic orbits (Standard Shooting)](@ref) for a simple example on how to use the above methods. 

The docs for this specific `newton` are located at [`newton`](@ref).

## Computation with `newton` and deflation

We also provide a simplified call to `newton` to locate the periodic orbit with a deflation operator:

```@docs
newton(prob::BifurcationKit.AbstractShootingDiscretization,
				orbitguess,
				options::NewtonPar;
				lens::BifurcationKit.OpticType = nothing,
				kwargs...)
```

and

```
newton(prob::BifurcationKit.AbstractShootingDiscretization,
				orbitguess::vectype,
				defOp::DeflationOperator{Tp, Tdot, T, vectype},
				options::NewtonPar{T, S, E};
				lens::OpticType = nothing,
				kwargs...,
			) where {T, Tp, Tdot, vectype, S, E}
```

## Continuation

Have a look at the [Continuation of periodic orbits (Standard Shooting)](@ref) example for the Brusselator.

In order to plot the orbit during continuation, one has to provide a `plot_solution` function (or rely on the default one) which recomputes the orbit from the current guess; this is simplified by the function `get_periodic_orbit` which returns the periodic orbit to be plotted. We refer to [Period doubling in Lur'e problem](@ref pdlure) for an example of use.

The docs for this specific `continuation` are located at [`continuation`](@ref).

```@docs
continuation(probPO::BifurcationKit.AbstractShootingDiscretization, orbitguess,
						alg::BifurcationKit.AbstractContinuationAlgorithm,
						contParams::ContinuationPar,
						linear_algo::BifurcationKit.AbstractBorderedLinearSolver;
						δ = convert(VI.scalartype(orbitguess), getdelta(probPO)),
						kwargs...,
					)
```

## References

[^Lust]:> **Numerical Bifurcation Analysis of Periodic Solutions of Partial Differential Equations**, Lust Kurt, 1997. 

[^Umbria]:> J. S. Umbría and M. Net. **Numerical continuation methods for large-scale dissipative dynamical systems**. The European Physical Journal Special Topics, 225(13):2465–2486, 2016.

[^Sanchez]:> Sánchez, J., M. Net, B. Garcı́a-Archilla, and C. Simó. “Newton–Krylov Continuation of Periodic Orbits for Navier–Stokes Flows.” Journal of Computational Physics 201, no. 1 (November 20, 2004): 13–33. https://doi.org/10.1016/j.jcp.2004.04.018.

[^Waugh]:> Waugh, Iain, Simon Illingworth, and Matthew Juniper. “Matrix-Free Continuation of Limit Cycles for Bifurcation Analysis of Large Thermoacoustic Systems.” Journal of Computational Physics 240 (May 2013): 225–47. https://doi.org/10.1016/j.jcp.2012.12.034.

