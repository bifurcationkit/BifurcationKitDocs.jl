# Simple Hopf point


At a Hopf branch point $(x_0,p_0)$ for the problem $F(x,p)=0$, the spectrum of the linear operator $dF(x_0,p_0)$ contains two purely imaginary $\pm i\omega,\ \omega > 0$ which are simple. At such point, we can compute the **normal form** to transform the initial Cauchy problem

$$\dot x = \mathbf{F}(x,p)$$

in large dimensions to a **complex** polynomial vector field ($\delta p\equiv p-p_0$):

$$\dot z = z\left(a \cdot\delta p + i\omega + l_1|z|^2\right)\quad\text{(E)}$$

whose solutions give access to the solutions of the Cauchy problem in a neighborhood of $(x,p)$.

More precisely, if $\mathbf{J} \equiv d\mathbf{F}(x_0,p_0)$, then we have $\mathbf{J}\zeta = i\omega\zeta$ and $\mathbf{J}\bar\zeta = -i\omega\bar\zeta$ for some complex eigenvector $\zeta$. It can be shown that $x(t) \approx x_0 + 2\Re(z(t)\zeta)$ when $p=p_0+\delta p$.

!!! tip "Coefficient $l_1$"
    The coefficient $l_1$ above is called the **Lyapunov** coefficient

### Expression of the coefficients

The coefficients $a,b$ above are computed as follows[^Haragus]:

$$a=\left\langle\mathbf{F}_{11}(\zeta)+2 \mathbf{F}_{20}\left(\zeta, \Psi_{001}\right), \zeta^{*}\right\rangle,$$

$$l_1=\left\langle 2 \mathbf{F}_{20}\left(\zeta, \Psi_{110}\right)+2 \mathbf{F}_{20}\left(\bar{\zeta}, \Psi_{200}\right)+3 \mathbf{F}_{30}(\zeta, \zeta, \bar{\zeta}), \zeta^{*}\right\rangle.$$

where

$$\begin{aligned}
-\mathbf{J} \Psi_{001} &=\mathbf{F}_{01} \\
(2 i \omega-\mathbf{J}) \Psi_{200} &=\mathbf{F}_{20}(\zeta, \zeta) \\
-\mathbf{J} \Psi_{110} &=2 \mathbf{F}_{20}(\zeta, \bar{\zeta}).
\end{aligned}$$

## Normal form computation

The normal form (E) is automatically computed as follows

```julia
getNormalForm(br::ContResult, ind_bif::Int ;
	verbose = false, ζs = nothing, lens = br.param_lens)
```

where `prob` is a bifurcation problem. `br` is a branch computed after a call to `continuation` with detection of bifurcation points enabled and `ind_bif` is the index of the bifurcation point on the branch `br`. The above call returns a point with information needed to compute the bifurcated branch. For more information about the optional parameters, we refer to [`getNormalForm`](@ref). The above call returns a point with information needed to compute the bifurcated branch.

```julia
mutable struct Hopf{Tv, T, Tω, Tevr, Tevl, Tnf} <: BifurcationPoint
	"Hopf point"
	x0::Tv

	"Parameter value at the Hopf point"
	p::T

	"Frequency of the Hopf point"
	ω::Tω

	"Right eigenvector"
	ζ::Tevr

	"Left eigenvector"
	ζstar::Tevl

	"Normal form coefficient (a = 0., b = 1 + 1im)"
	nf::Tnf

	"Type of Hopf bifurcation"
	type::Symbol
end
```

!!! info "Note"
    You should not need to call `getNormalForm ` except if you need the full information about the branch point.

## References

[^Haragus]: > Haragus, Mariana, and Gérard Iooss. Local Bifurcations, Center Manifolds, and Normal Forms in Infinite-Dimensional Dynamical Systems. London: Springer London, 2011. https://doi.org/10.1007/978-0-85729-112-7.
