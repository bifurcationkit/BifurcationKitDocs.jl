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

> This normal form is usually computed in order to branch from a Hopf-Hopf bifurcation point to curves of Neimark-Sacker bifurcations of periodic orbits (see [^Kuznetsov2]). Passing `detailed = Val(false)` returns only the data needed for this branching procedure while `detailed = Val(true)` (the default through [`get_normal_form`](@ref)) returns the cubic coefficients $G_{2100}, G_{0021}, G_{1011}, G_{1110}$ of (E) together with the coefficients required by the predictors.

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

`br` is a branch computed after a call to [`continuation`](@ref) with detection of bifurcation points enabled and `ind_bif` is the index of the bifurcation point on the branch `br`. The above call returns a point with information needed to compute the bifurcated branch. For more information about the optional parameters (`nev`, `ζs`, `scaleζ`, ...), we refer to [`get_normal_form`](@ref). The result returns an object of type `HopfHopf`.

!!! info "Note"
    You should not need to call `get_normal_form` except if you need the full information about the branch point.

### Returned object

The call `get_normal_form(br, ind_bif)` returns a `HopfHopf` point with the following fields

- `x0`, `params`, `lens`: the bifurcation point, the full parameter set and the two parameter axes,
- `ζ = (; q1, q2)` (resp. `ζ★ = (; p1, p2)`): the two complex eigenvectors of the Hopf pairs (resp. the left vectors), normalized by `scaleζ` and such that $\langle p_i, q_i\rangle = 1$. The pairs are ordered so that $\operatorname{imag}(\lambda_1) \geq \operatorname{imag}(\lambda_2) > 0$,
- `nf`: a named tuple holding
  - `ω0`: the frequency recorded by the continuation (i.e. that of the pair from which the Hopf-Hopf point was detected),
  - `λ1`, `λ2`: the eigenvalues of the two Hopf pairs (with the ordering above),
  - the cubic coefficients of (E): `G2100`, `G0021` (self couplings) and `G1011`, `G1110` (cross couplings),
  - the data required by the `:NS` predictor: `γ₁₁₀, γ₁₀₁, γ₂₁₀, γ₂₀₁`, the $2\times 2$ matrix `Γ`, the homological-equation terms `h₁₁₀₀, h₀₀₁₁, h₂₀₀₀, h₀₀₂₀`, the parameter terms `h₀₀₀₀₁₀, h₀₀₀₀₀₁` and the two sets `ns1 = (; dω1, dω2, α)` and `ns2 = (; dω1, dω2, α)`.

## Predictors

The predictor for a non trivial guess at distance $\delta p$ from the bifurcation point is provided by the methods

```@docs
BifurcationKit.predictor(hh::BifurcationKit.HopfHopf, ::Val{:HopfCurve}, ds::T; verbose = false, ampfactor = T(1)) where T
```

```@docs
BifurcationKit.predictor(hh::BifurcationKit.HopfHopf, ::Val{:NS}, ϵ::T; verbose = false, ampfactor = T(1)) where T
```

!!! tip "Detailed normal form"
    The `:NS` predictor gives an approximation of the two curves of Neimark-Sacker points of periodic orbits and requires the **detailed** normal form (`detailed = Val(true)`, the default). The `:HopfCurve` predictor (used to continue the curve of Hopf points) only needs `λ1`, `λ2`, `ω0` and the eigenvectors.

## References


[^Kuznetsov]:> Kuznetsov, Yu. A. “Numerical Normalization Techniques for All Codim 2 Bifurcations of Equilibria in ODE’s.” SIAM Journal on Numerical Analysis 36, no. 4 (January 1, 1999): 1104–24. https://doi.org/10.1137/S0036142998335005.

[^Kuznetsov2]:> Kuznetsov, Yu A., H. G. E. Meijer, W. Govaerts, and B. Sautois. “Switching to Nonhyperbolic Cycles from Codim 2 Bifurcations of Equilibria in ODEs.” Physica D: Nonlinear Phenomena 237, no. 23 (December 2008): 3061–68. https://doi.org/10.1016/j.physd.2008.06.006.
