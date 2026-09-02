# Predictors - Correctors

```@contents
Pages = ["Predictors.md"]
Depth = 3
```

The continuation method works with the following pattern (see [^Allgower1990]):

1. compute **tangent**
2. call **predictor** (based on tangent, mostly)
3. call **corrector**


There are several couples predictor-tangent/corrector which can be used in **BifurcationKit.jl** as we now explain. The tangent computation is formally included in the predictor whereas it is a distinct function in the code.

!!! info "Corrector"
    Note that setting the predictor also sets the corresponding corrector: it selects the couple predictor-corrector. You don't have (in fact cannot) set them independently.

## 1. Natural, zeroth order predictor

This is the simplest predictor based on the formula $(x_1,p_1) = (x_0, p_0 + ds)$ with plain Newton corrector ; it fails at Turning points because the Jacobian of $F$ becomes singular there. This is set by the algorithm `Natural()` in [`continuation`](@ref). With matrix based jacobians, this predictor brings little saving compared to the pseudo-arclength one: both require a (direct) factorization of the Jacobian at each continuation step, the PALC corrector merely adding a bordering of small size. For matrix-free (iterative) linear solvers, the Natural predictor can be faster than the following ones because the tangent computation is trivial, but only until it hits a Turning point.

## 2. First order predictor

This predictor is based on a computation of the tangent $\tau = (dx,dp)$ to the curve of solutions, it is given by $(x_1,p_1) = (x_0,p_0) + ds\cdot \tau$. This predictor passes Turning points when used with PALC Newton corrector.
**BifurcationKit.jl** provides two ways to compute the tangent $(dx, dp)$.

### 2a. Secant predictor
This predictor is called **secant** and is selected by the option `PALC(tangent = Secant())` in [`continuation`](@ref), see [`Secant`](@ref). It is computed from the last two points on the branch $\tau = (x_1, p_1) - (x_0, p_0)$ and then normalized with the weighted norm
$$\|(x, p)\|^2_\theta := \frac{\theta}{\mathtt{length}(x)} \langle x,x\rangle + (1 - \theta)\cdot p^2,\qquad 0<\theta<1.$$

!!! warning "Parameter `θ`"
    The parameter `θ` in the norm above is the field `θ` of [`PALC`](@ref) (default `θ = 0.5`). It is very important and should be tuned for the continuation to work properly, especially for large problems where the $\langle x - x_0, dx_0\rangle$ component in the [Pseudo arclength  continuation](@ref) constraint might be favored too much. Indeed, the term $1-\theta$ multiplies $p$ so that large `θ`s favour a parametrization by `p` and small `θ`s favour a parametrization by `x`.

### 2b. Bordered predictor
This predictor departs from the previous one in the way the tangent $\tau$ is estimated.
It computes $\tau:=(dx, dp)$ by solving the bordered linear system
$$\begin{bmatrix} F_x & F_p	\\ \frac{\theta}{\mathtt{length}(x)}dx_0 & (1-\theta)dp_0\end{bmatrix}\begin{bmatrix}dx \\  dp\end{bmatrix} =\begin{bmatrix}0 \\ 1\end{bmatrix}$$

where $\tau_0:=(dx_0, dp_0)$ is the tangent at the previous continuation step.

The predictor is selected by the option `PALC(tangent = Bordered())` in [`continuation`](@ref), see [`Bordered`](@ref). The linear solver for the bordered system in $(dx, dp)$ is set by the field `bls` of [`PALC`](@ref): it is one of [Bordered linear solvers (BLS)](@ref). It is also the linear solver used by the PALC corrector (see [Pseudo arclength  continuation](@ref)).

## 3. Polynomial predictor

The polynomial predictor is based on a fit (least square regression) of an $n$-th-order polynomial $P$ on the last $k$ solution vectors, where $n < k$. The arclength $s$ is used for the polynomial which then fits the solution $(x_i,p_i,s_i)$ as $P(s_i)\approx (x_i,p_i)$. To keep $s$ in suitable range (see [^Waugh]), we rescale it as $s\to \frac{s-\bar s}{\sigma}$ where $\sigma$ is the standard deviation of the $s_i$.

This predictor is selected by `PALC(tangent = Polynomial(pred, n, k, v0))` where `pred::AbstractTangentComputation` is the tangent predictor used for the first $k$ steps before the polynomial regression is operational and `v0` is an example of guess (used to set the types). More information is available in [`Polynomial`](@ref).

## 4. Multiple predictor (aka `pmcont` in `pde2path`)

The predictor is designed [^Uecker2014] to avoid spurious branch switching and pass singular points especially in PDE where branch point density can be quite high. It is based on the use of many predictors with increasing "jumps"
$$(x_i,p_i) = (x_0,p_0) + i\cdot ds\cdot \tau,\ i\leq nb$$
and use a corrector (PALC Newton) with the following twist. The criterion is that in each Newton step, the residual has to decrease by a factor $0<\alpha<1$:

$$\| F(u_n,p_n)\|\leq \alpha \| F(u_{n-1},p_{n-1}) \|$$

otherwise the corrector fails. The solution that is returned is the one for the highest $i$ for which the corrector succeeded. We refer to [^Uecker2014] for an exposition of the step size adaption strategy.

This algorithm is selected by `alg = Multiple(pred, x0, α, nb)` where `pred::PALC` is the underlying predictor-corrector used (the default value is `pred = PALC()`) and `x0` is an initial guess whose type is used to set the types of the tangents. More information is available in [`Multiple`](@ref).

## References

[^Allgower1990]: > Allgower and Georg, Numerical Continuation Methods, 1990

[^Uecker2014]: > 1.Uecker, H. pde2path - A Matlab Package for Continuation and Bifurcation in 2D Elliptic Systems. NMTMA 7, 58–106 (2014).

[^Waugh]: > Waugh, Illingworth, and Juniper, “Matrix-Free Continuation of Limit Cycles for Bifurcation Analysis of Large Thermoacoustic Systems.”
