# 🟢 Temperature model

> This is a classical example from the **Trilinos** library.

This is a simple example in which we aim at solving $\Delta T+\alpha N(T,\beta)=0$ with boundary conditions $T(0) = T(1)=\beta$. This example is coded in `examples/chan.jl`. We start with some imports:

```@example TUT1ODE
using BifurcationKit, Plots
const BK = BifurcationKit

N(x; a = 0.5, b = 0.01) = 1 + (x + a*x^2)/(1 + b*x^2)
nothing #hide
```

We then write our functional:

```@example TUT1ODE
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

```@example TUT1ODE
n = 101
sol0 = [(i-1)*(n-i)/n^2+0.1 for i=1:n]

# set of parameters
par = (α = 3.3, β = 0.01)
nothing #hide
```

Finally, we need to provide some parameters for the Newton iterations. This is done by calling

```@example TUT1ODE
optnewton = NewtonPar(tol = 1e-9, max_iterations = 10)
nothing #hide
```

We call the Newton solver:

```@example TUT1ODE
prob = BifurcationProblem(F_chan, sol0, par, (@optic _.α),
	# function to plot the solution
	plot_solution = (x, p; k...) -> plot!(x; ylabel="solution", label="", k...))
sol = solve(prob, Newton(), optnewton) # hide
# we set verbose to true to see the newton iterations
sol = @time solve(prob, Newton(), @set optnewton.verbose = true)
nothing #hide
```

Note that, in this case, we did not give the Jacobian. It was computed internally using Automatic Differentiation.

We can perform numerical continuation w.r.t. the parameter $\alpha$. This time, we need to provide additional parameters, but now for the continuation method:

```@example TUT1ODE
optcont = ContinuationPar(max_steps = 150,
	p_min = 0., p_max = 4.2,
	newton_options =optnewton)
nothing #hide
```

Next, we call the continuation routine as follows.

```@example TUT1ODE
br = continuation(prob, PALC(), optcont; plot = true)
nothing #hide		
```

The parameter axis `lens = @optic _.α` is used to extract the component of `par` corresponding to `α`. Internally, it is used as `get(par, lens)` which returns `3.3`.

!!! tip "Tip"
    We don't need to call `newton` first in order to use `continuation`.

You should see

```@example TUT1ODE
scene = title!("") #hide		
```

The left figure is the norm of the solution as function of the parameter $p=\alpha$, the *y-axis* can be changed by passing a different `record_from_solution` to `BifurcationProblem `. The top right figure is the value of $\alpha$ as function of the iteration number. The bottom right is the solution for the current value of the parameter. This last plot can be modified by changing the argument `plot_solution` to `BifurcationProblem `.

!!! note "Bif. point detection"
    Two Fold points were detected. This can be seen by looking at `show(br)` or `br.specialpoint`, by the black	dots on the continuation plots when doing `plot(br, plotfold=true)` or by typing `br` in the REPL. Note that the bifurcation points are located in `br.specialpoint`.

What if we want to compute to continue both ways in one call?

```@example TUT1ODE
br = continuation(prob, PALC(), optcont; bothside = true)
plot(br)
```