# Normal form of the Bautin bifurcation

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
get_normal_form(br::ContResult, ind_bif::Int;
    verbose = false, 
    ζs = nothing, 
    lens = getlens(br),
    kwargs...)
```

`br` is a branch computed after a call to [`continuation`](@ref) with detection of bifurcation points enabled and `ind_bif` is the index of the bifurcation point on the branch `br`. The above call returns a point with information needed to compute the bifurcated branch. For more information about the optional parameters, we refer to [`get_normal_form`](@ref). The result returns an object of type `Bautin`.

!!! info "Note"
    You should not need to call `get_normal_form` except if you need the full information about the branch point.

## Predictor

The predictor for a non trivial guess at distance $\delta p$ from the bifurcation point is provided by the method

```@docs
predictor(gh::BifurcationKit.Bautin, ::Val{:FoldPeriodicOrbitCont}, ϵ::T; verbose = false, ampfactor = T(1)) where T
```

## References


[^Kuznetsov]:> Kuznetsov, Yu. A. “Numerical Normalization Techniques for All Codim 2 Bifurcations of Equilibria in ODE’s.” SIAM Journal on Numerical Analysis 36, no. 4 (January 1, 1999): 1104–24. https://doi.org/10.1137/S0036142998335005.
