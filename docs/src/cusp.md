# Normal form of the Cusp bifurcation

We follow the paper[^Kuznetsov] and consider a Cauchy problem

$$\dot x=\mathbf F(x,p).$$

We denote by $\mathbf L$ the jacobian of $\mathbf F$ at the bifurcation point $(x_0,p_0)$. We choose a basis such that:

$$\mathbf L q=0, \quad \mathbf L^{T} p=0, \quad \langle p, q\rangle=1.$$


Under some conditions, $x(t)\approx x_0+ w(t)q$ where $w$ satisfies the normal form:

$$\dot{w}=c w^{3}+O\left(w^{4}\right).\tag{E}$$

## Normal form computation

The normal form (E) can be automatically computed as follows

```julia
get_normal_form(br::ContResult, ind_bif::Int ; 
    ζs = nothing,
    lens = getlens(br)
    )
```

`br` is a branch computed after a call to [`continuation`](@ref) with detection of bifurcation points enabled and `ind_bif` is the index of the bifurcation point on the branch `br`. The above call returns a point with information needed to compute the bifurcated branch. For more information about the optional parameters, we refer to [`get_normal_form`](@ref). The result returns an object of type `Cusp`.

!!! info "Note"
    You should not need to call `get_normal_form` except if you need the full information about the branch point.

## References

[^Kuznetsov]:> Kuznetsov, Yu. A. “Numerical Normalization Techniques for All Codim 2 Bifurcations of Equilibria in ODE’s.” SIAM Journal on Numerical Analysis 36, no. 4 (January 1, 1999): 1104–24. https://doi.org/10.1137/S0036142998335005.
