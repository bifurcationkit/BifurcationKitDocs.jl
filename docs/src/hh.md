# Normal form of the Hopf-Hopf bifurcation

We follow the paper[^Kuznetsov],[^Kuznetsov2] and consider a Cauchy problem

$$\dot x=\mathbf F(x,p).$$

We denote by $\mathbf L$ the jacobian of $\mathbf F$ at the bifurcation point $(x_0,p_0)$. We choose a basis such that:

$$\mathbf L q_1=i \omega_{1} q_1, \quad \mathbf L q_2=i \omega_{2} q_2.$$

Under some conditions, $x(t)\approx x_0+2\Re w_1(t)q_1+2\Re w_2(t)q_2$ where $w_i$ satisfy the normal form:

$$\left\{\begin{aligned}
\dot{w}_1= & i \omega_1 w_1+\frac{1}{2} G_{2100} w_1\left|w_1\right|^2+G_{1011} w_1\left|w_2\right|^2 
 +\frac{1}{12} G_{3200} w_1\left|w_1\right|^4+\frac{1}{2} G_{2111} w_1\left|w_1\right|^2\left|w_2\right|^2+\frac{1}{4} G_{1022} w_1\left|w_2\right|^4 \\
& +O\left(\left\|\left(w_1, \bar{w}_1, w_2, \bar{w}_2\right)\right\|^6\right) \\
\dot{w}_2= & i \omega_2 w_2+G_{1110} w_2\left|w_1\right|^2+\frac{1}{2} G_{0021} w_2\left|w_2\right|^2 +\frac{1}{4} G_{2210} w_2\left|w_1\right|^4+\frac{1}{2} G_{1121} w_2\left|w_1\right|^2\left|w_2\right|^2+\frac{1}{12} G_{0032} w_2\left|w_2\right|^4 \\
& +O\left(\left\|\left(w_1, \bar{w}_1, w_2, \bar{w}_2\right)\right\|^6\right)
\end{aligned}\right.\tag{E}$$

> This normal form is usually computed in order to branch from a Hopf-Hopf bifurcation point to curves of Neimark-Sacker bifurcations of periodic orbits (see [^Kuznetsov2]). Not all coefficients in (E) are required for this branching procedure, that is why only a subset of the $G_{ijkl}$ is returned.

## Normal form computation

The normal form (E) can be automatically computed as follows

```julia
get_normal_form(br::ContResult, ind_bif::Int ; verbose = false, ζs = nothing, lens = getlens(br))
```

`br` is a branch computed after a call to [`continuation`](@ref) with detection of bifurcation points enabled and `ind_bif` is the index of the bifurcation point on the branch `br`. The above call returns a point with information needed to compute the bifurcated branch. For more information about the optional parameters, we refer to [`get_normal_form`](@ref). The result returns an object of type `HopfHopf`.

!!! info "Note"
    You should not need to call `get_normal_form` except if you need the full information about the branch point.

## Predictors

The predictor for a non trivial guess at distance $\delta p$ from the bifurcation point is provided by the methods

```@docs
BifurcationKit.predictor(hh::BifurcationKit.HopfHopf, ::Val{:HopfCurve}, ds::T; verbose = false, ampfactor = T(1)) where T
```

```@docs
BifurcationKit.predictor(hh::BifurcationKit.HopfHopf, ::Val{:NS}, ϵ::T; verbose = false, ampfactor = T(1)) where T
```

## References


[^Kuznetsov]:> Kuznetsov, Yu. A. “Numerical Normalization Techniques for All Codim 2 Bifurcations of Equilibria in ODE’s.” SIAM Journal on Numerical Analysis 36, no. 4 (January 1, 1999): 1104–24. https://doi.org/10.1137/S0036142998335005.

[^Kuznetsov2]:> Kuznetsov, Yu A., H. G. E. Meijer, W. Govaerts, and B. Sautois. “Switching to Nonhyperbolic Cycles from Codim 2 Bifurcations of Equilibria in ODEs.” Physica D: Nonlinear Phenomena 237, no. 23 (December 2008): 3061–68. https://doi.org/10.1016/j.physd.2008.06.006.


