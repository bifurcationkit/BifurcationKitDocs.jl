# Periodic orbits based on orthogonal collocation

We compute `Ntst` time slices of a periodic orbit using orthogonal collocation. This is implemented in the structure `PeriodicOrbitOCollProblem`.

!!! warning "Large scale"
    The current implementation is not yet optimized for large scale problems. This will be improved in the future.    

The general method is very well exposed in [^Dankowicz],[^Doedel] and we adopt the notations of [^Dankowicz]. However our implementation is based on [^Doedel] because it is more economical (less equations) when it enforces the continuity of the solution.

We look for periodic orbits as solutions $(x(0), T)$ of

$$\dot x = T\cdot F(x),\ x(0)=x(1)\in\mathbb R^n.$$

We focus on the differential equality and consider a partition of the time domain

$$0=\tau_{1}<\cdots<\tau_{j}<\cdots<\tau_{N_{tst}+1}=1$$

where the points are referred to as **mesh points**. On each mesh interval $[\tau_j,\tau_{j+1}]$ for $j=1,\cdots,N_{tst}$, we define the affine transformation

$$\tau=\tau^{(j)}(\sigma):=\tau_{j}+\frac{(1+\sigma)}{2}\left(\tau_{j+1}-\tau_{j}\right), \sigma \in[-1,1].$$

The functions $x^{(j)}$ defined on $[-1,1]$ by $x^{(j)}(\sigma) \equiv x(\tau_j(\sigma))$ satisfies the following equation on $[-1,1]$:

$$\dot x^{(j)} = T\frac{\tau_{j+1}-\tau_j}{2}\cdot F(x^{(j)})\tag{$E_j$}$$

with the continuity equation $x^{(j+1)}(-1) = x^{(j)}(1)$.

We now aim at  solving $(E_j)$ by using an approximation with a polynomial of degree $m$. Following [^Dankowicz], we define a (uniform) partition:

$$-1=\sigma_{1}<\cdots<\sigma_{i}<\cdots<\sigma_{m+1}=1.$$

The points $\tau_{i,j} = \tau^{(i)}(\sigma_j)$ are called the **base points**: they serve as collocation points.

The associated $m+1$ Lagrange polynomials of degree $m$ are:

$$\mathcal{L}_{i}(\sigma):=\prod_{k=1, k \neq i}^{m+1} \frac{\sigma-\sigma_{k}}{\sigma_{i}-\sigma_{k}}, i=1, \ldots, m+1.$$

We then introduce the approximation $p_j$ of $x^{(j)}$:

$$\mathcal p_j(\sigma)\equiv \sum\limits_{k=1}^{m+1}\mathcal L_k(\sigma)x_{j,k}$$

and the problem to be solved at the **nodes** $z_l$, $l=1,\cdots,m$:

$$\forall 1\leq l\leq m,\quad 1\leq j\leq N_{tst},\quad \dot p_j(z_l) = T\frac{\tau_{j+1}-\tau_j}{2}\cdot F(p_j(z_l))\tag{$E_j^2$}.$$

The **nodes** $(z_l)$ are associated with a Gauss–Legendre quadrature.

In order to have a unique solution, we need to remove the phase freedom. This is done by imposing a *phase* condition.

## Number of unknowns

Putting the period unknown aside, we have to find the $x_{j,k}$ which gives $n\times N_{tst}\times (m+1)$ unknowns. 

The equations $E_j^2$ provides $n\times N_{tst}\times m$ plus the $(N_{tst}-1)\times n$ equations for the continuity equations. This makes a total of $(N_{tst}-1)\times m\times n+n\times N_{tst}\times m = n[N_{tst}(m+1)-1]$ equations to which we add the $n$ equations for the periodic boundary condition. In total, we have

$$n\times N_{tst}\times (m+1)$$

equations which matches the number of unknowns.

## Phase condition

To ensure uniqueness of the solution to the functional, we add the following phase condition

$$\frac{1}{T} \int_{0}^{T}\left\langle x(s), \dot x_0(s)\right\rangle d s =0$$

> During continuation at step $k$, we use $\frac{1}{T} \int_{0}^{T}\left\langle x(s), \dot x_{k-1}(s)\right\rangle d s$

## Interpolation

```@docs
BifurcationKit.POSolution
```

## Mesh adaptation
    
The goal of this method[^Russell] is to adapt the mesh $\tau_i$ in order to minimize the error.

## Encoding of the functional

The functional is encoded in the composite type [`PeriodicOrbitOCollProblem`](@ref). See the link for more information, in particular on how to access the underlying functional, its jacobian...

## Floquet multipliers computation

We provide three methods to compute the Floquet coefficients.

- The algorithm (Default) `FloquetColl` is based on the condensation of parameters described in [^Doedel]. It is the fastest method.
- The algorithm `FloquetCollGEV` is a simplified version of the procedure described in [^Fairgrieve]. It boils down to solving a large generalized eigenvalue problem. There is clearly room for improvements here but this can be used to check the results of the previous method.

These methods allow to detect bifurcations of periodic orbits. It seems to work reasonably well for the tutorials considered here. However they may be imprecise[^Lust].

- The state of the art method is based on a Periodic Schur decomposition. It is available through the package [PeriodicSchurBifurcationKit.jl](https://github.com/bifurcationkit/PeriodicSchurBifurcationKit.jl). For more information, have a look at `FloquetPQZ`.


## Computation with `newton`

We provide a simplified call to `newton` to locate the periodic orbits. Compared to the regular `newton` function, there is an additional option `jacobianPO` to select one of the many ways to deal with the jacobian of the above problem. The default solver `jacobianPO` is `:autodiffDense`.

The docs for this specific `newton` are located at [`newton`](@ref).

```@docs
newton(prob::PeriodicOrbitOCollProblem, orbitguess, options::NewtonPar; kwargs...)
```

## Continuation

We refer to [`continuation`](@ref) for more information regarding the arguments.

## References

[^Dankowicz]:> Dankowicz, Harry, and Frank Schilder. Recipes for Continuation. Computational Science and Engineering Series. Philadelphia: Society for Industrial and Applied Mathematics, 2013.

[^Doedel]:> Doedel, Eusebius, Herbert B. Keller, and Jean Pierre Kernevez. “NUMERICAL ANALYSIS AND CONTROL OF BIFURCATION PROBLEMS (II): BIFURCATION IN INFINITE DIMENSIONS.” International Journal of Bifurcation and Chaos 01, no. 04 (December 1991): 745–72.

[^Fairgrieve]:> Fairgrieve, Thomas F., and Allan D. Jepson. “O. K. Floquet Multipliers.” SIAM Journal on Numerical Analysis 28, no. 5 (October 1991): 1446–62. https://doi.org/10.1137/0728075.

[^Russell]:> Russell, R. D., and J. Christiansen. “Adaptive Mesh Selection Strategies for Solving Boundary Value Problems.” SIAM Journal on Numerical Analysis 15, no. 1 (February 1978): 59–80. https://doi.org/10.1137/0715004.

[^Lust]:> Lust, Kurt. “Improved Numerical Floquet Multipliers.” International Journal of Bifurcation and Chaos 11, no. 09 (September 2001): 2389–2410. https://doi.org/10.1142/S0218127401003486.

