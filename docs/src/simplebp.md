# Simple bifurcation branch point

!!! unknown "References"
    The general method is exposed in Golubitsky, Martin, David G Schaeffer, and Ian Stewart. **Singularities and Groups in Bifurcation Theory**. New York: Springer-Verlag, 1985, VI.1.d page 295

A simple branch point $(x_0,p_0)$ for the problem $F(x,p)=0$ satisfies $\dim \ker dF(x_0,p_0) = 1$. At such point, we can apply **Lyapunov-Schmidt** reduction to transform the initial problem in large dimensions to a **scalar** polynomial ($\delta p \equiv p-p_0$): 

$$a_{01}\delta p + \frac{a_{02}}{2}\delta p^2 + z\left(b_{11}\delta p + \frac{b_{20}}{2}z + \frac{b_{30}}{6}z^2\right) = 0 \tag{E}$$

whose solutions give access to all solutions in a neighborhood of $(x,p)$.

More precisely, if $\ker dF(x_0,p_0) = \mathbb R\zeta$, one can show that $x_0+z\zeta$ is close to a solution on a new branch, thus satisfying $F(x_0+z\zeta,p_0+\delta p)\approx 0$.

In the above scalar equation,

- if $a_{01}\neq 0$, this is a *Saddle-Node* bifurcation
- if $a_{01}=0, b_{20}\neq 0$, the bifurcation point is *Transcritical* and the bifurcated branch exists on each side of $p_0$.
- if $a_{01}=0,b_{20}=0, b_{30}\neq 0$, the bifurcation point is a *Pitchfork* and the bifurcated branch only exists on one side of $p_0$.

## Reduced equation computation

The reduced equation (E) can be automatically computed as follows

```julia
get_normal_form(br, ind_bif;
	verbose = false, lens = getlens(br),
	start_with_eigen = Val(true),    # Val(false): kernel basis via bordered systems
	bls = MatrixBLS(), bls_adjoint = bls)
```

where `prob` is the bifurcation problem. `br` is a branch computed after a call to [`continuation`](@ref) with detection of bifurcation points enabled and `ind_bif` is the index of the bifurcation point on the branch `br`. The above call returns a point with information needed to compute the bifurcated branch. For more information about the optional parameters (`nev`, `ζs`, `scaleζ`, ...), we refer to [`get_normal_form`](@ref).

!!! info "Note"
    You should not need to call `get_normal_form` except if you need the full information about the branch point.

!!! info "Note"
    Strictly speaking, this is not a normal form but a reduced equation. However, to keep the API simple, we use the same name especially for `get_normal_form` because it also returns the Hopf/BT/... normal forms which are true normal forms.

### Returned object

Depending on the coefficients of the reduced equation (E), the returned object is a `Fold` (saddle-node, $a_{01}\neq 0$), a `Transcritical` / `Pitchfork` ($a_{01}\approx 0$, with the corresponding `...Map` variant for maps) or a `BranchPoint` point. In all cases it has the following fields

- `x0`, `p`, `params`, `lens`: the bifurcation point, its parameter value, the full parameter set and the parameter axis,
- `ζ` (resp. `ζ★`): the real right (resp. left) vector spanning the kernel of the jacobian (resp. of its adjoint), normalized so that $\langle\zeta,\zeta*\rangle = 1$,
- `nf`: the coefficients of the reduced equation (E)
  - `a01`, `a02`: coefficients of $\delta p$ and $\delta p^2/2$,
  - `b11`, `b20`, `b30`: coefficients of $z\,\delta p$, $z^2/2$ and $z^3/6$,
  - `Ψ01`, `Ψ20`: the second order terms of the Lyapunov-Schmidt parametrization of the solutions $x \approx x_0 + z\zeta + \delta p\,\Psi_{01} + \frac{z^2}{2}\Psi_{20}$.

## Predictor

The predictor for a non trivial guess at distance $\delta p$ from the bifurcation point is provided by the methods (depending on the type of the bifurcation point)

```@docs
BifurcationKit.predictor(bp::Union{BifurcationKit.Transcritical, BifurcationKit.TranscriticalMap}, ds::T; verbose = false, ampfactor = T(1)) where T
```

```@docs
BifurcationKit.predictor(bp::Union{BifurcationKit.Pitchfork, BifurcationKit.PitchforkMap}, ds::T; verbose = false, ampfactor = T(1)) where T
```
