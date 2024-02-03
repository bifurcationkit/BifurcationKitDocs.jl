# Simple bifurcation branch point

!!! unknown "References"
    The general method is exposed in Golubitsky, Martin, David G Schaeffer, and Ian Stewart. **Singularities and Groups in Bifurcation Theory**. New York: Springer-Verlag, 1985, VI.1.d page 295

A simple branch point $(x_0,p_0)$ for the problem $F(x,p)=0$ satisfies $\dim \ker dF(x_0,p_0) = 1$. At such point, we can apply **Lyapunov-Schmidt** reduction to transform the initial problem in large dimensions to a **scalar** polynomial ($\delta p \equiv p-p_0$): 

$$a\delta p + z\left(b_1\delta p + \frac{b_2}{2}z + \frac{b_3}{6}z^2\right) = 0 \tag{E}$$

whose solutions give access to all solutions in a neighborhood of $(x,p)$.

More precisely, if $\ker dF(x_0,p_0) = \mathbb R\zeta$, one can show that $x_0+z\zeta$ is close to a solution on a new branch, thus satisfying $F(x_0+z\zeta,p_0+\delta p)\approx 0$.

In the above scalar equation,

- if $a\neq 0$, this is a *Saddle-Node* bifurcation
- if $a=0,b_2\neq 0$, the bifurcation point is *Transcritical* and the bifurcated branch exists on each side of $p_0$.
- if $a=0,b_2=0, b_3\neq 0$, the bifurcation point is a *Pitchfork* and the bifurcated branch only exists on one side of $p_0$.

## Normal form computation

The reduced equation (E) can be automatically computed as follows

```julia
get_normal_form(br::ContResult, ind_bif::Int ;
	verbose = false, ζs = nothing, lens = getlens(br))
```

where `prob` is the bifurcation problem. `br` is a branch computed after a call to [`continuation`](@ref) with detection of bifurcation points enabled and `ind_bif` is the index of the bifurcation point on the branch `br`. The above call returns a point with information needed to compute the bifurcated branch. For more information about the optional parameters, we refer to [`get_normal_form`](@ref). The result returns an object of type `BranchPoint`.

!!! info "Note"
    You should not need to call `get_normal_form` except if you need the full information about the branch point.

## Predictor

The predictor for a non trivial guess at distance $\delta p$ from the bifurcation point is provided by the methods (depending on the type of the bifurcation point)

```@docs
BifurcationKit.predictor(bp::Union{BifurcationKit.Transcritical, BifurcationKit.TranscriticalMap}, ds::T; verbose = false, ampfactor = T(1)) where T
```

```@docs
BifurcationKit.predictor(bp::Union{BifurcationKit.Pitchfork, BifurcationKit.PitchforkMap}, ds::T; verbose = false, ampfactor = T(1)) where T
```
