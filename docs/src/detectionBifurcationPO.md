# Detection of bifurcation points of periodic orbits

The bifurcations are detected during a call to `br = continuation(prob, alg, contParams::ContinuationPar;kwargs...)` by turning on the following flags:

- `contParams.detect_bifurcation = 2`

The bifurcation points are located by looking at the spectrum **e.g.** by monitoring the unstable eigenvalues. The Floquet exponent λ is declared unstable if `real(λ) > contParams.tol_stability`. The located bifurcation points are then returned in `br.specialpoint`. 
    
## Precise detection of bifurcation points using bisection    

Note that the bifurcation points detected when `detect_bifurcation = 2` can be rather *crude*  localization of the true bifurcation points. Indeed, we only signal that, in between two continuation steps *which can be large*, a (several) bifurcation has been detected. Hence, we only have a rough idea of where the bifurcation is located, unless your `dsmax` is very small... This can be improved as follows.

If you choose `detect_bifurcation = 3`, a **bisection algorithm** is used to locate the bifurcation points more precisely. It means that we recursively track down the change in stability. Some options in [`ContinuationPar`](@ref) control this behavior:

- `n_inversion`: number of sign inversions in the bisection algorithm
- `max_bisection_steps` maximum number of bisection steps
- `tol_bisection_eigenvalue` tolerance on real part of Floquet exponent to detect bifurcation points in the bisection steps

!!! tip "Bisection mode"
    During the bisection, the eigensolvers are called like `eil(J, nev; bisection = true)` in order to be able to adapt the solver precision.

## Large scale computations

The user must specify the number of eigenvalues to be computed (like `nev = 10`) in the parameters `::ContinuationPar` passed to `continuation`. Note that `nev` is automatically incremented whenever a bifurcation point is detected [^1]. Also, there is an option in `::ContinuationPar` to save (or not) the eigenvectors. This can be useful in memory limited environments (like on GPUs).
    
[^1]: In this case, the Krylov dimension is not increased because the eigensolver could be a direct solver. You might want to increase this dimension using the callbacks in [`continuation`](@ref). 

## List of detected bifurcation points
|Bifurcation|index used|
|---|---|
| Bifurcation point (single eigenvalue stability change, Fold or branch point) | bp |
| Neimark-Sacker | ns |
| Period doubling | pd |
| Not documented | nd |

## Eigensolver

The user must provide an eigensolver by setting `newton_options.eigsolver` where `newton_options` is located in the parameter `::ContinuationPar` passed to continuation. See [`NewtonPar`](@ref) and [`ContinuationPar`](@ref) for more information on the composite type of the options passed to `newton` and `continuation`.

The eigensolver is highly problem dependent and this is why the user should implement / parametrize its own eigensolver through the abstract type `AbstractEigenSolver` or select one among [List of implemented eigen solvers](@ref).

!!! danger "Floquet multipliers computation"
    The computation of Floquet multipliers is necessary for the detection of bifurcations of periodic orbits (which is done by analyzing the Floquet exponents obtained from the Floquet multipliers). Hence, the eigensolver needs to compute the eigenvalues with largest modulus (and not with largest real part which is their default behavior). This can be done by changing the option `which = :LM` of the eigensolver. Nevertheless, note that for most implemented eigensolvers in `BifurcationKit`, the proper option is automatically set.   

## Generic bifurcation

By this we mean a change in the dimension of the Jacobian kernel. The detection of Branch point is done by analysis of the spectrum of the Jacobian.

The detection is triggered by setting `detect_bifurcation > 1` in the parameter `::ContinuationPar` passed to `continuation`. 

## Fold bifurcation
The detection of **Fold** point is done by monitoring  the monotonicity of the parameter.

The detection is triggered by setting `detect_fold = true` in the parameter `::ContinuationPar` passed to `continuation`. When a **Fold** is detected on a branch `br`, a point is added to `br.foldpoint` allowing for later refinement using the function `newton_fold`.

## Neimark-Sacker bifurcation

The detection of Neimark-Sacker point is done by analysis of the spectrum of the Jacobian.

The detection is triggered by setting `detect_bifurcation > 1` in the parameter `::ContinuationPar` passed to `continuation`. When a **Neimark-Sacker point** is detected, a point is added to `br.specialpoint`.

## Period-doubling bifurcation

The detection of Period-doubling point is done by analysis of the spectrum of the Jacobian.

The detection is triggered by setting `detect_bifurcation > 1` in the parameter `::ContinuationPar` passed to `continuation`. When a **Period-doubling point** is detected, a point is added to `br.specialpoint`.

