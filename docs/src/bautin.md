# Normal form of the Bautin bifurcation

```@contents
Pages = ["bautin.md"]
Depth = 2
```

We follow the paper[^Kuznetsov] and consider a Cauchy problem

$$\dot x=\mathbf F(x,p).$$

We denote by $\mathbf L$ the jacobian of $\mathbf F$ at the bifurcation point $(x_0,p_0)$. We choose a basis such that:

$$\mathbf L q=i \omega_{0} q, \quad \mathbf L^{T} p=-i \omega_{0} p, \quad \langle p, q\rangle=1.$$

Under some conditions, $x(t)\approx x_0+2\Re w(t)q$ where $w$ satisfies the normal form:

$$\dot{w}=i \omega_{0} w+\frac{1}{2} G_{21} w|w|^{2}+\frac{1}{12} G_{32} w|w|^{4}+O\left(|w|^{6}\right).\tag{E}$$

The second Lyapunov coefficient is 

$$l_2:=\frac{1}{12} \operatorname{Re} G_{32}.$$ 

## Normal form computation

The normal form (E) can be automatically computed as follows

```julia
get_normal_form(br, ind_bif;
    verbose = false, lens = getlens(br),
    detailed = Val(true),            # full normal form
    start_with_eigen = Val(true),    # Val(false): kernel basis via bordered systems
    bls = MatrixBLS(), bls_adjoint = bls)
```

`br` is a branch computed after a call to [`continuation`](@ref) with detection of bifurcation points enabled and `ind_bif` is the index of the bifurcation point on the branch `br`. The above call returns a point with information needed to compute the bifurcated branch. For more information about the optional parameters (`nev`, `ζs`, `scaleζ`, ...), we refer to [`get_normal_form`](@ref). The result returns an object of type `Bautin`.

!!! info "Note"
    You should not need to call `get_normal_form` except if you need the full information about the branch point.

### Returned object

The call `get_normal_form(br, ind_bif)` returns a `Bautin` point with the following fields

- `x0`, `params`, `lens`: the bifurcation point, the full parameter set and the two parameter axes,
- `ζ` (resp. `ζ★`): the complex right (resp. left) eigenvector of the Hopf pair, normalized by `scaleζ` and such that $\langle \zeta*, \zeta\rangle = 1$,
- `nf`: a named tuple holding
  - `ω`: the frequency of the Hopf pair,
  - `G21`, `G32`: the complex cubic and quintic coefficients of (E),
  - `l1 = G21/2`, `l2 = real(G32)/12`: the first and second Lyapunov coefficients (the second Lyapunov coefficient is that of (E)),
  - in detailed mode, the additional data required to branch to the curve of folds of periodic orbits: `h₂₀₀₀, h₁₁₀₀, h₀₀₁₀, h₀₀₀₁` (homological equation terms), `γ₁₁₀, γ₁₀₁, γ₂₁₀, γ₂₀₁` and `α` (unfolding coefficients, see [^Kuznetsov]).

The predictor below relies on the detailed normal form: call `get_normal_form(br, ind_bif; detailed = Val(true))` (the default).

## Predictor

The predictor for a non trivial guess at distance $\delta p$ from the bifurcation point is provided by the method

```@docs
predictor(gh::BifurcationKit.Bautin, ::Val{:FoldPeriodicOrbitCont}, ϵ::T; verbose = false, ampfactor = T(1)) where T
```

## References


[^Kuznetsov]:> Kuznetsov, Yu. A. “Numerical Normalization Techniques for All Codim 2 Bifurcations of Equilibria in ODE’s.” SIAM Journal on Numerical Analysis 36, no. 4 (January 1, 1999): 1104–24. https://doi.org/10.1137/S0036142998335005.
