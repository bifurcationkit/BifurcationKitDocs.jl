# Bogdanov-Takens refinement

```@contents
Pages = ["codim3Continuation.md"]
Depth = 2
```

In this page, we explain how to perform precise localisation of Bogdanov-Takens (BT) points. This is an unusual feature of numerical continuation libraries. We chose to implement it because the localisation of the BT points on the Hopf bifurcation curves is rather imprecise.


## Method

The continuation of BT bifurcation points is based on a **Minimally Augmented**[^Govaerts],[^Blank],[^Bindel] formulation which is an efficient way to detect singularities. The continuation of BT points is based on the formulation

$$G(u,p) = (F(u,p), g_1(u,p), g_2(u,p))\in\mathbb R^{n+2}\quad\quad (F_{bt})$$

where the test functions $g_1,g_2$ are solutions of

$$\left[\begin{array}{cc}
dF(u,p) & w \\
v^{\top} & 0
\end{array}\right]\left[\begin{array}{c}
v_1 \\
g_1(u,p)
\end{array}\right]=\left[\begin{array}{c}0_{n} \\1\end{array}\right]\quad\quad (M_{bt})$$

and

$$\left[\begin{array}{cc}
dF(u,p) & w \\
v^{\top} & 0
\end{array}\right]\left[\begin{array}{c}
v_2 \\
g_2(u,p)
\end{array}\right]=\left[\begin{array}{c}v_1 \\0\end{array}\right]\quad\quad (M_{bt})$$

and where $w,v$ are chosen in order to have a non-singular matrix $(M_{bt})$. More precisely, $v$ (resp. $w$) should be close to a null vector of `dF(u,p)` (resp. `dF(u,p)'`).

> note that there are very simplified calls for this, see **Newton refinement** below. In particular, you don't need to set up the Minimally Augmented problem yourself. This is done in the background.

!!! warning "Linear Method"
    You can pass the bordered linear solver to solve $(M_{bt})$ using the option `bdlinsolver ` (see below). Note that the choice `bdlinsolver = BorderingBLS()` can lead to singular systems. Indeed, in this case, $(M_{bt})$ is solved by inverting `dF(u,p)` which is singular at Fold points.

## Setting the jacobian

In order to apply the newton algorithm to $F_{bt}$, one needs to invert the jacobian. This is not completely trivial as one must compute this jacobian and then invert it. You can select the following jacobians for your computations (see below):

- `jacobian_ma = AutoDiff()`: [Default] automatic differentiation is applied to $F_{bt}$ and the matrix is then inverted using the provided linear solver. In particular, the jacobian is formed. This is very well suited for small dimensions  (say < 100)
- `jacobian_ma = MinAug()`: a specific procedure for evaluating the jacobian $F_{bt}$ and inverting it (without forming the jacobian!) is used. This is well suited for large dimensions.

## Example

```@example CODIM3
using Revise, BifurcationKit
Fbt(x, p) = [x[2], p.β1 + p.β2 * x[2] + p.a * x[1]^2 + p.b * x[1] * x[2]]
par = (β1 = 0.01, β2 = -0.3, a = -1., b = 1.)
prob  = BifurcationProblem(Fbt, [0.01, 0.01], par, (@optic _.β1))
opts_br = ContinuationPar(p_max = 0.5, p_min = -0.5, detect_bifurcation = 3, nev = 2)

br = continuation(prob, PALC(), opts_br; bothside = true)

# compute branch of Hopf points
hopf_codim2 = continuation(br, 3, (@optic _.β2), ContinuationPar(opts_br, max_steps = 40) ;
	detect_codim2_bifurcation = 2,
	update_minaug_every_step = 1,
	bothside = true,
	)

# refine BT point
solbt = BifurcationKit.newton_bt(hopf_codim2, 2)
solbt.u
```


## Newton refinement

Once a Bogdanov-Takens point has been detected after a call to `br = continuation(...)`, it can be refined using `newton` iterations. Let us say that `ind_bif` is the index in `br.specialpoint` of a Bogdanov-Takens point. This guess can be refined as follows:

```julia
outfold = newton(br::AbstractBranchResult, ind_bif::Int;  
	normN = norm,
	options = br.contparams.newton_options,
	bdlinsolver = BorderingBLS(options.linsolver),
	jacobian_ma = AutoDiff(),
	start_with_eigen = false, kwargs...)
```

For the options parameters, we refer to [newton](@ref).

It is important to note that for improved performances, a function implementing the expression of the **hessian** should be provided. This is by far the fastest. `BifurcationProblem` provides it by default using AD though.

## Advanced use

Here, we expose the solvers that are used to perform newton refinement. This is useful in case it is too involved to expose the linear solver options.

```@docs
BifurcationKit.newton_bt
```

## References

[^Govaerts]: > Govaerts, Willy J. F. Numerical Methods for Bifurcations of Dynamical Equilibria. Philadelphia, Pa: Society for Industrial and Applied Mathematics, 2000.

[^Blank]: > Blank, H. J. de, Yu. A. Kuznetsov, M. J. Pekkér, and D. W. M. Veldman. “Degenerate Bogdanov–Takens Bifurcations in a One-Dimensional Transport Model of a Fusion Plasma.” Physica D: Nonlinear Phenomena 331 (September 15, 2016): 13–26. https://doi.org/10.1016/j.physd.2016.05.008.

[^Bindel]: > Bindel, D., M. Friedman, W. Govaerts, J. Hughes, and Yu.A. Kuznetsov. “Numerical Computation of Bifurcations in Large Equilibrium Systems in Matlab.” Journal of Computational and Applied Mathematics 261 (May 2014): 232–48. https://doi.org/10.1016/j.cam.2013.10.034.
