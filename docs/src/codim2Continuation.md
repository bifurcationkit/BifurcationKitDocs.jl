# Fold / Hopf Continuation

In this page, we explain how to perform continuation of Fold / Hopf points and detect the associated bifurcations.

For this to work best, it is necessary to have an analytical expression for the jacobian. See the tutorial [Temperature model (Simplest example)](@ref) for more details.

A quite complete example for detection of codim 2 bifurcations of equilibria is [CO-oxydation (codim 2)](@ref) although it is for ODEs.

### List of detected codim 2 bifurcation points
|Bifurcation|symbol used|
|---|---|
| Bogdanov-Takens | bt |
| Bautin | gh |
| Cusp | cusp |
| Zero-Hopf | zh |
| Hopf-Hopf | hh |

In a nutshell, all you have to do (see below) is to call `continuation(br, ind_bif)` to continue the bifurcation point stored in `br.specialpoint[ind_bif]` and set proper options. 

## Fold continuation

The continuation of Fold bifurcation points is based on a **Minimally Augmented**[^Govaerts] formulation which is an efficient way to detect singularities. The continuation of Fold points is based on the formulation $G(u,p) = (F(u,p), g(u,p))\in\mathbb R^{n+1}$ where the test function $g$ is solution of

$$\left[\begin{array}{cc}
dF(u,p) & w \\
v^{\top} & 0
\end{array}\right]\left[\begin{array}{c}
r \\
g(u,p)
\end{array}\right]=\left[\begin{array}{c}0_{n} \\1\end{array}\right]\quad\quad (M_f)$$

where $w,v$ are chosen in order to have a non-singular matrix $(M_f)$. More precisely, $v$ (resp. $w$) should be close to a null vector of `dF(u,p)` (resp. `dF(u,p)'`). During continuation, the vectors $w,v$ are updated so that the matrix $(M_f)$ remains non-singular ; this is controlled with the argument `updateMinAugEveryStep` (see below).

> note that there are very simplified calls for this, see **Newton refinement** below. In particular, you don't need to set up the Fold Minimally Augmented problem yourself. This is done in the background.

!!! warning "Linear Method"
    You can pass the bordered linear solver to solve $(M_f)$ using the option `bdlinsolver ` (see below). Note that the choice `bdlinsolver = BorderingBLS()` can lead to singular systems. Indeed, in this case, $(M_f)$ is solved by inverting `dF(u,p)` which is singular at Fold points.
    
### Detection of codim 2 bifurcation points

You can detect the following codim 2 bifurcation points by using the option `detectCodim2Bifurcation` in the method `continuation` (see [Codim 2 continuation](@ref)). Under the hood, the detection of these bifurcations is done by using Event detection as explained in [Event Handling](@ref).

- the detection of Cusp (Cusp) is done by the detection of Fold bifurcation points along the curve of Folds by monitoring the parameter component of the tangent.
- the detection of Bogdanov-Takens (BT) is performed using the test function[^Bindel] $\psi_{BT}(p) = \langle w(p),v(p)\rangle$
- the detection of Zero-Hopf (ZH) is performed by monitoring the number of eigenvalues $\lambda$ such that $\Re\lambda > \min\limits_{\nu\in\Sigma(dF)}|\Re\nu|$ and $\Im\lambda > \epsilon$ where $\epsilon$ is the Newton tolerance.

## Hopf continuation

The continuation of Fold bifurcation points is based on a **Minimally Augmented** (see [^Govaerts] p. 87) formulation which is an efficient way to detect singularities. The continuation of Hopf points is based on the formulation $G(u,\omega,p) = (F(u,\omega,p), g(u,\omega,p))\in\mathbb R^{n+2}$ where the test function $g$ is solution of

$$\left[\begin{array}{cc}
dF(u,p)-i\omega I_n & w \\
v^{\top} & 0
\end{array}\right]\left[\begin{array}{c}
r \\
g(u,\omega,p)
\end{array}\right]=\left[\begin{array}{c}
0_{n} \\
1
\end{array}\right]\quad\quad (M_h)$$

where $w,v$ are chosen in order to have a non-singular matrix $(M_h)$. More precisely, $w$ (resp. $v$) should be a left (resp. right) approximate null vector of $dF(u,p)-i\omega I_n$. During continuation, the vectors $w,v$ are updated so that the matrix $(M_h)$ remains non-singular ; this is controlled with the argument `updateMinAugEveryStep ` (see below).

> note that there are very simplified calls to this, see **Newton refinement** below. In particular, you don't need to set up the Hopf Minimally Augmented problem yourself. This is done in the background.

!!! warning "Linear Method"
    You can pass the bordered linear solver to solve $(M_h)$ using the option `bdlinsolver ` (see below). Note that the choice `bdlinsolver = BorderingBLS()` can lead to singular systems. Indeed, in this case, $(M_h)$ is solved by inverting `dF(u,p)-i?? I_n` which is singular at Hopf points.
    
### Detection of codim 2 bifurcation points

You can detect the following codim 2 bifurcation points by using the option `detectCodim2Bifurcation` in the method `continuation` (see [Codim 2 continuation](@ref)). Under the hood, the detection of these bifurcations is done by using Event detection as explained in [Event Handling](@ref).

- the detection of Bogdanov-Takens (BT) is performed using the test function[^Bindel] $\psi_{BT}(p) = 	\langle w(p),v(p)\rangle$
- the detection of Bautin (GH) is based on the test function $\psi_{GH}(p) = \Re(l_1(p))$ where $l_1$ is the Lyapunov coefficient defined in [Simple Hopf point](@ref).
- the detection of Zero-Hopf (ZH) is performed by monitoring the eigenvalues.
- the detection of Hopf-Hopf (HH) is performed by monitoring the eigenvalues.

> The continuation of Hopf points is stopped at BT and when $\omega<100\epsilon$ where $\epsilon$ is the newton tolerance. 

## Newton refinement

Once a Fold / Hopf point has been detected after a call to `br = continuation(...)`, it can be refined using `newton` iterations. Let us say that `ind_bif` is the index in `br.specialpoint` of a Fold / Hopf point. This guess can be refined as follows:

```julia
outfold = newton(br::AbstractBranchResult, ind_bif::Int;  
	normN = norm, options = br.contparams.newtonOptions, 
	bdlinsolver = BorderingBLS(options.linsolver),
	startWithEigen = false, kwargs...)
```

For the options parameters, we refer to [Newton](@ref).

It is important to note that for improved performances, a function implementing the expression of the **hessian** should be provided. This is by far the fastest. Reader interested in this advanced usage should look at the code `example/chan.jl` of the tutorial [Temperature model (Simplest example)](@ref). 

## Codim 2 continuation

To compute the codim 2 curve of Fold / Hopf points, one can call [`continuation`](@ref) with the following options

```@docs
 continuation(br::BifurcationKit.AbstractBranchResult, ind_bif::Int64,
				lens2::Lens, options_cont::ContinuationPar = br.contparams ;
				kwargs...)
```

where the options are as above except with have an additional parameter axis `lens2` which is used to locate the bifurcation points. 


See [Temperature model (Simplest example)](@ref) for an example of use. 

## Advanced use

Here, we expose the solvers that are used to perform newton refinement or codim 2 continuation in case the above methods fails. This is useful in case it is too involved to expose the linear solver options. An example of advanced use is the continuation of Folds of periodic orbits, see [Continuation of Fold of periodic orbits](@ref).

```@docs
newtonFold
```

```@docs
newtonHopf
```


```@docs
continuationFold
```

```@docs
continuationHopf
```

## References

[^Govaerts]: > Govaerts, Willy J. F. Numerical Methods for Bifurcations of Dynamical Equilibria. Philadelphia, Pa: Society for Industrial and Applied Mathematics, 2000.

[^Blank]: > Blank, H. J. de, Yu. A. Kuznetsov, M. J. Pekk??r, and D. W. M. Veldman. ???Degenerate Bogdanov???Takens Bifurcations in a One-Dimensional Transport Model of a Fusion Plasma.??? Physica D: Nonlinear Phenomena 331 (September 15, 2016): 13???26. https://doi.org/10.1016/j.physd.2016.05.008.

[^Bindel]: > Bindel, D., M. Friedman, W. Govaerts, J. Hughes, and Yu.A. Kuznetsov. ???Numerical Computation of Bifurcations in Large Equilibrium Systems in Matlab.??? Journal of Computational and Applied Mathematics 261 (May 2014): 232???48. https://doi.org/10.1016/j.cam.2013.10.034.



