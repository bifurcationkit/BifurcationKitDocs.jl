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
computeNormalForm(F, dF, d2F, d3F, br::ContResult, ind_bif::Int ; δ = 1e-8,
	nev = 5, Jᵗ = nothing, verbose = false, ζs = nothing, autodiff = true, detailed = true)
```

where `dF, d2F, d3F` are the differentials of `F`. `br` is a branch computed after a call to [`continuation`](@ref) with detection of bifurcation points enabled and `ind_bif` is the index of the bifurcation point on the branch `br`. The option `detailed` controls the computation of a simplified version of the normal form. `autodiff` controls the use of `ForwardDiff` during the normal form computation.


The above call returns a point with information needed to compute the bifurcated branch. For more information about the optional parameters, we refer to [`computeNormalForm`](@ref). The result returns the following:

```julia
mutable struct BogdanovTakens{Tv, Tpar, Tlens, Tevr, Tevl, Tnf, Tnf2} <: AbstractBifurcationPoint
	"Bogdanov-Takens point"
	x0::Tv

	"Parameters used by the vector field."
	params::Tpar

	"Parameter axis used to compute the branch on which this BT point was detected."
	lens::Tlens

	"Right eigenvectors"
	ζ::Tevr

	"Left eigenvectors"
	ζstar::Tevl

	"Normal form coefficients (basic)"
	nf::Tnf

	"Normal form coefficients (detailed)"
	nfsupp::Tnf2

	"Type of bifurcation"
	type::Symbol
end
```

!!! info "Note"
    You should not need to call `computeNormalForm` except if you need the full information about the branch point.

## References

[^Haragus]:> Haragus, Mariana, and Gérard Iooss. Local Bifurcations, Center Manifolds, and Normal Forms in Infinite-Dimensional Dynamical Systems. London: Springer London, 2011. https://doi.org/10.1007/978-0-85729-112-7.


[^AlHdaibat]:> Al-Hdaibat, B., W. Govaerts, Yu. A. Kuznetsov, and H. G. E. Meijer. “Initialization of Homoclinic Solutions near Bogdanov--Takens Points: Lindstedt--Poincaré Compared with Regular Perturbation Method.” SIAM Journal on Applied Dynamical Systems 15, no. 2 (January 2016): 952–80. https://doi.org/10.1137/15M1017491.
