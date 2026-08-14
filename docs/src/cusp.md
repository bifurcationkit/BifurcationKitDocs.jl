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
get_normal_form(br, ind_bif;
    verbose = false, lens = getlens(br),
    start_with_eigen = Val(true),    # Val(false): kernel basis via bordered systems
    bls = MatrixBLS(), bls_adjoint = bls)
```

`br` is a branch computed after a call to [`continuation`](@ref) with detection of bifurcation points enabled and `ind_bif` is the index of the bifurcation point on the branch `br`. The above call returns a point with information needed to compute the bifurcated branch. For more information about the optional parameters (`nev`, `ζs`, `scaleζ`, ...), we refer to [`get_normal_form`](@ref). The result returns an object of type `Cusp`.

!!! info "Note"
    You should not need to call `get_normal_form` except if you need the full information about the branch point.

### Returned object

The call `get_normal_form(br, ind_bif)` returns a `Cusp` point with the following fields

- `x0`, `params`, `lens`: the bifurcation point, the full parameter set and the two parameter axes,
- `ζ` (resp. `ζ★`): the real right (resp. left) vector spanning the kernel of $\mathbf L$ (resp. of $\mathbf L^{T}$), normalized so that $\langle \zeta,\zeta*\rangle = 1$,
- `nf`: a named tuple holding
  - `nf.c`: the coefficient $c$ of (E).

For completeness, the coefficient is computed as

$$c = \frac{1}{6}\left(\langle p, \mathbf C(q,q,q)\rangle + 3\langle p, \mathbf B(q, H_2)\rangle\right)$$

where $\mathbf B$, $\mathbf C$ are the second and third derivatives of $\mathbf F$ at $(x_0,p_0)$ and $H_2$ solves the bordered system $\mathbf L H_2 = \langle p, \mathbf B(q,q)\rangle\, q - \mathbf B(q,q)$ with $\langle p, H_2\rangle = 0$.

## References

[^Kuznetsov]:> Kuznetsov, Yu. A. “Numerical Normalization Techniques for All Codim 2 Bifurcations of Equilibria in ODE’s.” SIAM Journal on Numerical Analysis 36, no. 4 (January 1, 1999): 1104–24. https://doi.org/10.1137/S0036142998335005.
