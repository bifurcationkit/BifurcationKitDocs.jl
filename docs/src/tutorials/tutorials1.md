# [🟡 Temperature model (codim 2)](@id temperature)

```@contents
Pages = ["tutorials1.md"]
Depth = 3
```

> This is a classical example from the **Trilinos** library.

This is a simple example in which we aim at solving $\Delta T+\alpha N(T,\beta)=0$ with boundary conditions $T(0) = T(1)=\beta$. This example is coded in `examples/chan.jl`. We start with some imports:

```@example TUT1
using Revise, BifurcationKit, LinearAlgebra, Plots
const BK = BifurcationKit

N(x; a = 0.5, b = 0.01) = 1 + (x + a*x^2)/(1 + b*x^2)
nothing #hide
```

We then write our functional:

```@example TUT1
function F_chan(x, p)
	(;α, β) = p
	f = similar(x)
	n = length(x)
	f[1] = x[1] - β
	f[n] = x[n] - β
	for i=2:n-1
		f[i] = (x[i-1] - 2 * x[i] + x[i+1]) * (n-1)^2 + α * N(x[i], b = β)
	end
	return f
end
nothing #hide
```
We want to call a Newton solver. We first need an initial guess:

```@example TUT1
n = 101
sol0 = [(i-1)*(n-i)/n^2+0.1 for i=1:n]

# set of parameters
par = (α = 3.3, β = 0.01)
nothing #hide
```

Finally, we need to provide some parameters for the Newton iterations. This is done by calling

```@example TUT1
optnewton = NewtonPar(tol = 1e-11, verbose = true)
nothing #hide
```

We call the Newton solver:

```@example TUT1
prob = BifurcationProblem(F_chan, sol0, par, (@optic _.α),
	# function to plot the solution
	plot_solution = (x, p; k...) -> plot!(x; ylabel="solution", label="", k...))
sol = BK.solve(prob, Newton(), @set optnewton.verbose=false) # hide
sol = @time BK.solve( prob, Newton(), optnewton)
nothing #hide
```

Note that, in this case, we did not give the Jacobian. It was computed internally using Automatic Differentiation.

We can perform numerical continuation w.r.t. the parameter $\alpha$. This time, we need to provide additional parameters, but now for the continuation method:

```@example TUT1
optcont = ContinuationPar(dsmin = 0.01, dsmax = 0.2, ds= 0.1, p_min = 0., p_max = 4.2,
	newton_options = NewtonPar(max_iterations = 10, tol = 1e-9))
nothing #hide
```

Next, we call the continuation routine as follows.

```@example TUT1
br = continuation(prob, PALC(), optcont; plot = true)
nothing #hide		
```

The parameter axis `lens = @optic _.α` is used to extract the component of `par` corresponding to `α`. Internally, it is used as `get(par, lens)` which returns `3.3`.

!!! tip "Tip"
    We don't need to call `newton` first in order to use `continuation`.

You should see

```@example TUT1
scene = title!("") #hide		
```

The left figure is the norm of the solution as function of the parameter $p=\alpha$, the *y-axis* can be changed by passing a different `recordFromSolution` to `BifurcationProblem `. The top right figure is the value of $\alpha$ as function of the iteration number. The bottom right is the solution for the current value of the parameter. This last plot can be modified by changing the argument `plotSolution` to `BifurcationProblem `.

!!! note "Bif. point detection"
    Two Fold points were detected. This can be seen by looking at `br.specialpoint`, by the black	dots on the continuation plots when doing `plot(br, plotfold=true)` or by typing `br` in the REPL. Note that the bifurcation points are located in `br.specialpoint`.


## Continuation of Fold points

We get a summary of the branch by doing

```@example TUT1
br
```

We can take the first Fold point, which has been guessed during the previous continuation run and locate it precisely. However, this only works well when the jacobian is computed analytically. We use automatic differentiation for that

```@example TUT1
# index of the Fold bifurcation point in br.specialpoint
indfold = 2

outfold = newton(
	#index of the fold point
	br, indfold)
BK.converged(outfold) && printstyled(color=:red, "--> We found a Fold Point at α = ", outfold.u.p, ", β = 0.01, from ", br.specialpoint[indfold].param,"\n")
```

We can finally continue this fold point in the plane $(α, β)$ by performing a Fold Point continuation. In the present case, we find a Cusp point.

!!! tip "Tip"
    We don't need to call `newton` first in order to use `continuation` for the codim 2 curve of bifurcation points.

```@example TUT1
outfoldco = continuation(br, indfold,
	# second parameter axis to use for codim 2 curve
	(@optic _.β),
	# we disable the computation of eigenvalues, it makes little sense here
	ContinuationPar(optcont, detect_bifurcation = 0))
scene = plot(outfoldco, plotfold = true, legend = :bottomright)
```

!!! tip "Tip"
    The performances for computing the curve of Fold is not that great. It is because we use the default solver tailored for ODE. If you pass `jacobian_ma = :minaug` to the last `continuation` call, you should see a great improvement in performances.

## Using GMRES or another linear solver

We continue the previous example but now using Matrix Free methods. The user can pass its own solver by implementing a version of `LinearSolver`. Some linear solvers have been implemented from `KrylovKit.jl` and `IterativeSolvers.jl` (see [Linear solvers (LS)](@ref) for more information), we can use them here. Note that we can also use preconditioners as shown below. The same functionality is present for the eigensolver.

```@example TUT1
# derivative of N
dN(x; a = 0.5, b = 0.01) = (1-b*x^2+2*a*x)/(1+b*x^2)^2

# Matrix Free version of the differential of F_chan
# Very easy to write since we have F_chan.
# We could use Automatic Differentiation as well
function dF_chan(x, dx, p)
	(;α, β) = p
	out = similar(x)
	n = length(x)
	out[1] = dx[1]
	out[n] = dx[n]
	for i=2:n-1
		out[i] = (dx[i-1] - 2 * dx[i] + dx[i+1]) * (n-1)^2 + α * dN(x[i], b = β) * dx[i]
	end
	return out
end

# we create a new linear solver
ls = GMRESKrylovKit(dim = 100)

# and pass it to the newton parameters
optnewton_mf = NewtonPar(verbose = true, linsolver = ls, tol = 1e-10)

# we change the problem with the new jacobian
prob = re_make(prob;
	# we pass the differential a x,
	# which is a linear operator in dx
	J = (x, p) -> (dx -> dF_chan(x, dx, p))
	)

# we can then call the newton solver
out_mf = BK.solve(prob, Newton(), @set optnewton_mf.verbose = false) # hide
out_mf = @time BK.solve(prob, Newton(), optnewton_mf)
nothing #hide
```

We can improve this computation, *i.e.* reduce the number of `Linear-Iterations`, by using a preconditioner

```@example TUT1
using SparseArrays

# define preconditioner which is basically Δ
P = spdiagm(0 => -2 * (n-1)^2 * ones(n), -1 => (n-1)^2 * ones(n-1), 1 => (n-1)^2 * ones(n-1))
P[1,1:2] .= [1, 0.];P[end,end-1:end] .= [0, 1.]

# define gmres solver with left preconditioner
ls = GMRESIterativeSolvers(reltol = 1e-4, N = length(sol.u), restart = 10, maxiter = 10, Pl = lu(P))
optnewton_mf = NewtonPar(verbose = true, linsolver = ls, tol = 1e-10)
out_mf = BK.solve(prob, Newton(), @set optnewton_mf.verbose = false) # hide
out_mf = @time BK.solve(prob, Newton(), optnewton_mf)
nothing #hide
```
