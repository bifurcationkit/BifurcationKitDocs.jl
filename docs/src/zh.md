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

> This normal form is usually computed in order to branch from a Zero-Hopf bifurcation point to curves of Neimark-Sacker bifurcations of periodic orbits (see [^Kuznetsov2]). Not all coefficients in (E) are required for this branching procedure, that is why only a subset of the $G_{ijk}$ is returned.

## Normal form computation

The normal form (E) can be automatically computed as follows

```julia
get_normal_form(br::ContResult, ind_bif::Int ;  
    ζs = nothing, 
    lens = getlens(br))
```

`br` is a branch computed after a call to [`continuation`](@ref) with detection of bifurcation points enabled and `ind_bif` is the index of the bifurcation point on the branch `br`. The above call returns a point with information needed to compute the bifurcated branch. For more information about the optional parameters, we refer to [`get_normal_form`](@ref). The result returns an object of type `ZeroHopf`.


!!! info "Note"
    You should not need to call `get_normal_form` except if you need the full information about the branch point.

## References


[^Kuznetsov]:> Kuznetsov, Yu. A. “Numerical Normalization Techniques for All Codim 2 Bifurcations of Equilibria in ODE’s.” SIAM Journal on Numerical Analysis 36, no. 4 (January 1, 1999): 1104–24. https://doi.org/10.1137/S0036142998335005.

[^Kuznetsov2]:> Kuznetsov, Yu A., H. G. E. Meijer, W. Govaerts, and B. Sautois. “Switching to Nonhyperbolic Cycles from Codim 2 Bifurcations of Equilibria in ODEs.” Physica D: Nonlinear Phenomena 237, no. 23 (December 2008): 3061–68. https://doi.org/10.1016/j.physd.2008.06.006.


