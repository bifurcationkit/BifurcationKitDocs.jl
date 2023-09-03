# FAQ

### How can I save a solution every n steps, or at specific parameter values?

You can use the callback `finalise_solution` in the function call `continuation`. For example, you can use something like this to save all steps

```julia
function mySave(u, tau, step, contResult, personaldata)
	push!(personaldata, u)
end
```
and pass it like `continuation(prob, alg, opts; finalise_solution = (z, tau, step, contResult; k...) -> mySave(z, tau, step, contResult, myData))`

### The Fold / Hopf Continuation does not work, why?

This requires some precise computations. Have you tried passing the expression of the Jacobian instead of relying on finite differences.

### What is the parameter `theta` about in `ContinuationPar`?

See the description of [continuation](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/library/#Continuation-1) on the page [Library](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/library/).

### How can I change the preconditioner during computations?

The easiest way to achieve this is by using the callbacks provided by `newton` and `continuation`. See the documentation about these two methods. See also the example [2d Ginzburg-Landau equation (finite differences, codim 2, Hopf aBS)](@ref cgl)

### How can I implement my own bifurcation detection method?

You can use the callback `finalise_solution` but the best way is probably to use the [Iterator Interface](@ref) to inject your code anywhere in the continuation procedure.

### How do I dissociate the computation of eigenvalues from the jacobian that I passed?

Sometimes, for example when implementing boundary conditions, you pass a jacobian `J` but the eigenvalues, and the bifurcation points are not simply related to `J`. One way to bypass this issue is to define a new eigensolver `<: AbstractEigenSolver` and pass it to the `NewtonPar` field `eigsolver`. This is done for example in `example/SH2d-fronts-cuda.jl`.

### How can I print the eigenvalues during `continuation`?

You can print the eigenvalues using the following callback:

```juliaw
finalise_solution = (z, tau, step, contResult; k...) -> begin
		BK.haseigenvalues(contResult) && Base.display(contResult.eig[end].eigenvals)
		return true
	end,
```

### How can I reject a Newton Step?

You can reject a newton step by passing to `continuation` the argument `callback_newton`

```julia
function mycallback((x, f, J, res, iteration, itlinear, options); kwargs...)
	# stop Newton algo if residual too large
	if res > 1e2
		@warn "Reject Newton step!!"
		return false
	end
	return true
end
```

### How do I stop `continuation`?

Using the argument `finalise_solution` in `continuation`. Simply make this function `finalise_solution` return false.

### How do I compute both sides of a branch?

Instead of using two calls to `continuation`, you can pass the keyword `bothside=true` to `continuation`

### How do I compute period orbits for non-autonomous problems

The package does not yet allow to compute periodic orbits solutions of non-autonomous Cauchy problems like

$$\frac{du}{dt}  = F(u, par, t).$$

On certain cases, one can still go away with the following trick. Say one is interested (dummy example!) to study

$$\dot u = cos(u) + cos(\omega \cdot t).$$

Then one can use the following autonomous vector field

```julia
function vectorField(U, par)
	u, x, y = U
	out = similar(U)
	out[1] = cos(u) + x
	x2 = x^2+y^2
	out[2] = x + par.ω * y - x * x2
	out[3] = y - par.ω * x - y * x2
	out
end
```

### Arpack is slow in computing eigenvalues

This is probably due to iterative refinement conducted by `SuiteSparse` as explained in this blog [post](https://discourse.julialang.org/t/some-eigenpairs-from-a-large-sparse-nonsymmetric-matrix-julia-vs-matlab/93742). You can disable this using

```julia
using SuiteSparse
SuiteSparse.UMFPACK.umf_ctrl[8] = 0
```

### Should I use CVODE_BDF?

SciML is now able to match the performance of the Sundials solver `CVODE_BDF`. Check the [news](https://sciml.ai/news/2021/05/24/QNDF/) for more information.
