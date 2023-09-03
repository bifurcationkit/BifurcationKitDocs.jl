# Continuation of Period-doubling (PD) bifurcations of periodic orbits

```@contents
Pages = ["PDContinuationPO.md"]
Depth = 2
```

!!! warning "Work in progress"
    Continuation of PD of POs is a recent addition. It has not been thoroughly tested.  

In this page, we explain how to perform continuation of PD points of periodic orbits and detect the following codim 2 bifurcations.

### List of detected codim 2 bifurcation points
|Bifurcation|symbol used|
|---|---|
| Strong resonance 1:2 bifurcation | R2 |
| Fold / Flip| foldFlip |
| Period-Doubling / Neimark-Sacker | pdNS |
| Generalized Period-Doubling | gpd |


In a nutshell, all you have to do (see below) is to call `continuation(br, ind_bif, lens2)` to continue the bifurcation point stored in `br.specialpoint[ind_bif]` and set proper options.

## PD continuation

The continuation of PD bifurcation points is based on a **Minimally Augmented**[^Govaerts] formulation which is an efficient way to detect singularities (see [Fold / Hopf Continuation](@ref)). All the methods (see [Periodic orbits computation](@ref)), except the Trapezoid one, for computing periodic orbits are compatible with this algorithm. In particular, you can perform these computations in large dimensions.


## Detection of codim 2 bifurcation points

You can detect the following codim 2 bifurcation points by using the option `detect_codim2_bifurcation` in the method `continuation` 

- the detection of Generalised Period-Doubling bifurcation is done by computing the PD normal form
- the detection the above bifurcation points is done by monitoring the number of eigenvalues $\lambda$ such that $\Re\lambda > \min\limits_{\nu\in\Sigma(dF)}|\Re\nu|$ and $\Im\lambda > \epsilon$ where $\epsilon$ is the Newton tolerance.

## Setting the jacobian

In order to apply the newton algorithm to the PD functional, one needs to invert the jacobian. This is not completely trivial as one must compute this jacobian and then invert it. You can select the following jacobians for your computations (see below):

- [Default] for `jacobian_ma = :autodiff`, automatic differentiation is applied to the PD functional and the matrix is then inverted using the provided linear solver. In particular, the jacobian is formed. This is very well suited for small dimensions  (say < 100)
- for `jacobian_ma = :finiteDifferences`, same as `jacobian_ma = :autodiff` but the jacobian is computed using finite differences.
- for `jacobian_ma = :minaug`, a specific procedure for evaluating the jacobian and inverting it (without forming the jacobian!) is used. This is well suited for large dimensions.

## Codim 2 continuation

To compute the codim 2 curve of PD points of periodic orbits, one can call [`continuation`](@ref) with the following options

```@docs
 continuation(br::BifurcationKit.AbstractBranchResult, ind_pd::Int64,
				lens2::Lens, options_cont::ContinuationPar = br.contparams ;
				kwargs...)
```

where `br` is a branch of periodic orbits and the options are as above except with have an additional parameter axis `lens2` which is used to locate the bifurcation points.

## Algorithmic details

The continuation of PD points is based on the formulation

$$G(u,p,\omega) = (F_{po}(u,p), \sigma(u,p))\in\mathbb R^{n+1}\quad\quad (\mathcal F_{pd})$$

where $F_{po}$ is the functional for locating periodic orbits and the test function $\sigma$ is solution of

$$\left[\begin{array}{cc}
N(u,p) & w \\
v^{\top} & 0
\end{array}\right]\left[\begin{array}{c}
r \\
\sigma(u,p)
\end{array}\right]=\left[\begin{array}{c}
0_{n} \\
1
\end{array}\right].$$

In the case of Multiple Standard Shooting, the matrix $N$ is based on the monodromy $M(x_i,T_i)$

$$N:=\left(\begin{array}{cccccc}
{M_1} & -I & {0} & {\cdots} & 0 \\
0 & {M_2} & -I & {\cdots} & {0} \\
{\vdots} & {\cdots} & {\ddots} & {\ddots} & {\vdots} \\
{0} & {\cdots} & {\cdots} & {\ddots} & -I \\
I & {\cdots} & {\cdots} & 0 & {M_{m}} \\
\end{array}\right).$$

The jacobian of the PD functional to use for the Newton algorithm

$$\left[\begin{array}{cc}
\partial_{u}F_{po} & \partial_pF_{po} \\
\partial_u\sigma & \partial_p\sigma
\end{array}\right].$$

## References

[^Govaerts]:> Govaerts, Willy J. F. Numerical Methods for Bifurcations of Dynamical Equilibria. Philadelphia, Pa: Society for Industrial and Applied Mathematics, 2000.
