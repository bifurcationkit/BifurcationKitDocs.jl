# Normal form of the Bogdanov-Takens bifurcation

We follow the book[^Haragus] and consider a Cauchy problem

$$\dot x=\mathbf F(x,p).$$

We denote by $\mathbf L$ the jacobian of $\mathbf F$ at the bifurcation point $(x_0,p_0)$. We choose a basis such that:

$$\mathbf{L} \zeta_{0}=0, \quad \mathbf{L} \zeta_{1}=\zeta_{0}.$$

We can also select a basis:

$$\mathbf{L}^{*} \zeta_{1}^{*}=0, \quad \mathbf{L}^{*} \zeta_{0}^{*}=\zeta_{1}^{*}$$

such that

$$\left\langle\zeta_{0}, \zeta_{0}^{*}\right\rangle=1, \quad\left\langle\zeta_{1}, \zeta_{0}^{*}\right\rangle=0, \quad\left\langle\zeta_{0}, \zeta_{1}^{*}\right\rangle=0, \quad\left\langle\zeta_{1}, \zeta_{1}^{*}\right\rangle=1.$$

Under some conditions, $x(t)\approx x_0+A(t)\zeta_0 + B(t)\zeta_1$ where $A,B$ satisfy the normal form:

$$\begin{aligned}
&\frac{d A}{d t}=B \\
&\frac{d B}{d t}=\alpha_{1}(\mu)+\alpha_{2}(\mu) A+\alpha_{3}(\mu) B+b A B+a A^{2}\widetilde{\rho}(A, B, \mu)
\end{aligned}\tag{E}$$

where $p = p_0+\mu$ and with coefficients

$$\begin{aligned}
&a=\left\langle\mathbf{F}_{20}\left(\zeta_{0}, \zeta_{0}\right), \zeta_{1}^{*}\right\rangle \\
&b=\left\langle 2 \mathbf{F}_{20}\left(\zeta_{0}, \zeta_{1}\right)-2 \Psi_{200}, \zeta_{1}^{*}\right\rangle.
\end{aligned}$$

The $\Psi$s satisfy

$$\begin{aligned}
a \zeta_{1} &=\mathbf{L} \Psi_{200}+\mathbf{F}_{20}\left(\zeta_{0}, \zeta_{0}\right) \\
b \zeta_{1}+2 \Psi_{200} &=\mathbf{L} \Psi_{110}+2 \mathbf{F}_{20}\left(\zeta_{0}, \zeta_{1}\right) \\
\Psi_{110} &=\mathbf{L} \Psi_{020}+\mathbf{F}_{20}\left(\zeta_{1}, \zeta_{1}\right)
\end{aligned}$$

which gives

$$0=\left\langle\Psi_{200}, \zeta_{1}^{*}\right\rangle + \left\langle\mathbf{F}_{20}\left(\zeta_{0}, \zeta_{0}\right), \zeta_{0}^{*}\right\rangle.$$

We conclude that

$$\begin{aligned}
&a=\left\langle\mathbf{F}_{20}\left(\zeta_{0}, \zeta_{0}\right), \zeta_{1}^{*}\right\rangle \\
&b=2\left\langle  \mathbf{F}_{20}\left(\zeta_{0}, \zeta_{1}\right), \zeta_{1}^{*}\right\rangle + 2\left\langle\mathbf{F}_{20}\left(\zeta_{0}, \zeta_{0}\right), \zeta_{0}^{*}\right\rangle.
\end{aligned}$$


### Computation of the basis

To build the basis $\left\{\zeta_{0}, \zeta_{1}\right\}$, we follow the procedure described in [^AlHdaibat] on page 972.

### Computation of the parameter transform

To invert the mapping $\mu\to (\alpha_{1}(\mu),\alpha_{2}(\mu),\alpha_{3}(\mu))$, we follow the procedure described in [^AlHdaibat] on page 956 forward.

## Normal form computation

The normal form (E) can be automatically computed as follows

```julia
get_normal_form(br, ind_bif;
    verbose = false, lens = getlens(br),
    detailed = Val(true),            # full normal form
    autodiff = true,                 # use ForwardDiff for the differentiations
    start_with_eigen = Val(true),    # Val(false): kernel basis via bordered systems
    bls = MatrixBLS(), bls_adjoint = bls, bls_block = bls)
```

`br` is a branch computed after a call to [`continuation`](@ref) with detection of bifurcation points enabled and `ind_bif` is the index of the bifurcation point on the branch `br`. The option `detailed` controls the computation of a simplified version of the normal form. `autodiff` controls the use of `ForwardDiff` during the normal form computation.

The above call returns a point with information needed to compute the bifurcated branch. For more information about the optional parameters (`nev`, `ζs`, `scaleζ`, ...), we refer to [`get_normal_form`](@ref). The result returns an object of type `BogdanovTakens`.

!!! info "Note"
    You should not need to call `get_normal_form` except if you need the full information about the branch point.

### Returned object

The call `get_normal_form(br, ind_bif)` returns a `BogdanovTakens` point with the following fields

- `x0`, `params`, `lens`: the bifurcation point, the full parameter set and the two parameter axes,
- `ζ = (; q0, q1)` (resp. `ζ★ = (; p0, p1)`): the real right (resp. left) generalized eigenvectors forming the Jordan chain $\mathbf L q_0 = 0,\ \mathbf L q_1 = q_0$ (resp. $\mathbf L^{T} p_1 = 0,\ \mathbf L^{T} p_0 = p_1$),
- `nf`: a named tuple holding
  - `a`, `b`: the quadratic coefficients of (E) computed above,
  - in detailed mode, the additional coefficients `γ, c, K10, K11, K2, d, e, a1, b1` of the truncated normal form and of the parameter transform $\mu \mapsto (\alpha_1,\alpha_2,\alpha_3)$ (see the two sections above and [^AlHdaibat]) as well as the center manifold terms `H0001, H0010, H0002, H1001, H2000` (homological equations, analogous to the $\Psi$s of the simple Hopf case).

The predictors below rely on the detailed normal form: call `get_normal_form(br, ind_bif; detailed = Val(true))` (the default).


## Predictors

The predictor for a non trivial guess at distance $\delta p$ from the bifurcation point is provided by the method

```@docs
BifurcationKit.predictor(bt::BifurcationKit.BogdanovTakens, ::Val{:HopfCurve}, ds::T; verbose = false, ampfactor = T(1)) where T
```

```@docs
BifurcationKit.predictor(bt::BifurcationKit.BogdanovTakens, ::Val{:FoldCurve}, ds::T; verbose = false, ampfactor = T(1)) where T
```

```@docs
BifurcationKit.predictor(bt::BifurcationKit.BogdanovTakens, ::Val{:HomoclinicCurve}, ds::T; verbose = false, ampfactor = one(T)) where T
```

## References

[^Haragus]:> Haragus, Mariana, and Gérard Iooss. Local Bifurcations, Center Manifolds, and Normal Forms in Infinite-Dimensional Dynamical Systems. London: Springer London, 2011. https://doi.org/10.1007/978-0-85729-112-7.


[^AlHdaibat]:> Al-Hdaibat, B., W. Govaerts, Yu. A. Kuznetsov, and H. G. E. Meijer. “Initialization of Homoclinic Solutions near Bogdanov--Takens Points: Lindstedt--Poincaré Compared with Regular Perturbation Method.” SIAM Journal on Applied Dynamical Systems 15, no. 2 (January 2016): 952–80. https://doi.org/10.1137/15M1017491.
