# Normal form of the Zero-Hopf bifurcation

We follow the paper[^Kuznetsov],[^Kuznetsov2] and consider a Cauchy problem

$$\dot x=\mathbf F(x,p).$$

We denote by $\mathbf L$ the jacobian of $\mathbf F$ at the bifurcation point $(x_0,p_0)$. We choose a basis such that:

$$\mathbf L q_0=0, \quad \mathbf L q_1=i \omega_{0} q_1.$$

Under some conditions, $x(t)\approx x_0+2\Re w_1(t)q_1+w_0(t)q_0$ where $w_i$ satisfy the normal form:

$$\left\{\begin{aligned}
\dot{w}_0= & \frac{1}{2} G_{200} w_0^2+G_{011}\left|w_1\right|^2+\frac{1}{6} G_{300} w_0^3  +G_{111} w_0\left|w_1\right|^2+O\left(\left\|\left(w_0, w_1, \bar{w}_1\right)\right\|^4\right) \\
\dot{w}_1= & i \omega_0 w_1+G_{110} w_0 w_1+\frac{1}{2} G_{210} w_0^2 w_1+\frac{1}{2} G_{021} w_1\left|w_1\right|^2  +O\left(\left\|\left(w_0, w_1, \bar{w}_1\right)\right\|^4\right) .
\end{aligned}\right.\tag{E}$$

> This normal form is usually computed in order to branch from a Zero-Hopf bifurcation point to curves of Neimark-Sacker bifurcations of periodic orbits (see [^Kuznetsov2]). The flag `hasNS` (see below) tells whether such a curve emanates from the point. Passing `detailed = Val(false)` returns only the data needed for this branching procedure while `detailed = Val(true)` (the default through [`get_normal_form`](@ref)) returns all the coefficients of (E).

## Normal form computation

The normal form (E) can be automatically computed as follows

```julia
get_normal_form(br, ind_bif;
    verbose = false, lens = getlens(br),
    detailed = Val(true),            # full normal form
    autodiff = true,                 # use ForwardDiff for the differentiations
    start_with_eigen = Val(true),    # Val(false): kernel basis via bordered systems
    bls = MatrixBLS(), bls_adjoint = bls)
```

`br` is a branch computed after a call to [`continuation`](@ref) with detection of bifurcation points enabled and `ind_bif` is the index of the bifurcation point on the branch `br`. The above call returns a point with information needed to compute the bifurcated branch. For more information about the optional parameters (`nev`, `ζs`, `scaleζ`, ...), we refer to [`get_normal_form`](@ref). The result returns an object of type `ZeroHopf`.

!!! info "Note"
    You should not need to call `get_normal_form` except if you need the full information about the branch point.

### Returned object

The call `get_normal_form(br, ind_bif)` returns a `ZeroHopf` point with the following fields

- `x0`, `params`, `lens`: the bifurcation point, the full parameter set and the two parameter axes,
- `ζ = (; q0, q1)` (resp. `ζ★ = (; p0, p1)`): the real null vector $q_0$ and the complex Hopf vector $q_1$ (resp. the left vectors), normalized so that $\langle q_0, p_0\rangle = \langle q_1, p_1\rangle = 1$,
- `nf`: a named tuple holding
  - `λ0`: the (real, null) eigenvalue $\mathbf L q_0 = 0$, ideally $\approx 0$,
  - `ω`: the imaginary frequency of the Hopf pair ($\omega = \operatorname{imag}(\lambda_1)$ where $\mathbf L q_1 = \lambda_1 q_1$, $\lambda_1 \approx i\omega_0$),
  - the coefficients of (E): `G200`, `G110`, `G011`, `G111`, `G021` together with the derived quantities `g110 = G110`, `f011 = G011` and the flag `hasNS` indicating whether a curve of Neimark-Sacker points emanates from the point,
  - the remaining data used by the predictors: the homological-equation terms `h200, h110, h020, h011, h00010, h00001` and `v10, v01, β1, β2, τ1, τ2, x`.

### Returned object (non-detailed)

Passing `detailed = Val(false)` returns only `nf = (; ω, λ0, dFp)` where `ω` is the (complex) eigenvalue $\lambda_1$ of the Hopf pair, together with `ζ` and `ζ★`. This is the minimal information required to start the continuation of the curves of equilibria (fold or Hopf curves) emanating from the point.

## Predictors

The predictor for a non trivial guess at distance $\delta p$ from the bifurcation point is provided by the methods

```@docs
BifurcationKit.predictor(zh::BifurcationKit.ZeroHopf, ::Val{:HopfCurve}, ds::T; verbose = false, ampfactor = T(1)) where T
```

```@docs
BifurcationKit.predictor(zh::BifurcationKit.ZeroHopf, ::Val{:FoldCurve}, ds::T; verbose = false, ampfactor = T(1)) where T
```

```@docs
BifurcationKit.predictor(zh::BifurcationKit.ZeroHopf, ::Val{:NS}, ϵ::T; verbose = false, ampfactor = T(1)) where T
```

!!! tip "Detailed normal form"
    The `:NS` predictor gives an approximation of the Neimark-Sacker curve of periodic orbits and requires the **detailed** normal form (`detailed = Val(true)`, the default). The predictors `:HopfCurve` and `:FoldCurve` only need the non-detailed data.

## References


[^Kuznetsov]:> Kuznetsov, Yu. A. “Numerical Normalization Techniques for All Codim 2 Bifurcations of Equilibria in ODE’s.” SIAM Journal on Numerical Analysis 36, no. 4 (January 1, 1999): 1104–24. https://doi.org/10.1137/S0036142998335005.

[^Kuznetsov2]:> Kuznetsov, Yu A., H. G. E. Meijer, W. Govaerts, and B. Sautois. “Switching to Nonhyperbolic Cycles from Codim 2 Bifurcations of Equilibria in ODEs.” Physica D: Nonlinear Phenomena 237, no. 23 (December 2008): 3061–68. https://doi.org/10.1016/j.physd.2008.06.006.
