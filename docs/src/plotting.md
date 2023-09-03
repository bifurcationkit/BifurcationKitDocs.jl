# Plotting

```@contents
Pages = ["plotting.md"]
Depth = 3
```

## Standard plots using the plot recipe from Plots.jl

Plotting is provided by calling recipes to `Plots.jl`. It means that to plot a branch `br`, you just need to call

```julia
#]add Plots # You need to install Plots.jl before your first time using it!
using Plots
plot(br)
```

where `br` is a branch computed after a call to `br = continuation(...)`. Plots can be customized using
[all the keyword arguments provided by Plots.jl](https://docs.juliaplots.org/stable/generated/attributes_plot/).
For example, we can change the plotting backend to the GR package and put a title
on the plot by doing:

```julia
gr()
plot!(br, title = "I have a branch!")
```
or you can use a scatter plot

```julia
scatter(br)
```

Then to save the plot, use `savefig`, for example:

```julia
savefig("myplot.png")
```

## Specific plotting keyword arguments

The available arguments specific to our plotting methods are

- `plotfold = true`: plot the fold points with black dots
- `putspecialptlegend = true`: display the legend corresponding to the bifurcation points
- `vars = nothing`: see below
- `plotstability = true`: display the stability of the branch
- `plotspecialpoints = true`: plot the special (bifurcation) points on the branch
- `branchlabel = "fold branch"`: assign label to a branch which is printed in the legend
- `linewidthunstable`: set the linewidth for the unstable part of the branch
- `linewidthstable`: set the linewidth for the stable part of the branch
- `plotcirclesbif = false` use circles to plot bifurcation points
- `applytoX = identity` apply transformation `applytoX` to x-axis
- `applytoY = identity` apply transformation `applytoY` to y-axis

If you have several branches `br1, br2`, you can plot them in the same figure by doing

```julia
plot(br1, br2)
```

in place of

```julia
plot(br1)
plot!(br2)
```

!!! warn "Plot of bifurcation points"
    The bifurcation points for which the bisection was successful are indicated with circles and with squares otherwise.

Note that the plot recipes use the parameter axis as `xlabel`, and the passed variable as `ylabel`.

## Choosing Variables

You can select which variables to plot using the keyword argument `vars`, for example:

```julia
plot(br, vars = (:param, :x))
```
The available symbols are `:x, :param, :itnewton, :itlinear, :ds, :θ, :n_unstable, :n_imag, :stable, :step`,... and:

- `x` if `record_from_solution` (see [`continuation`](@ref)) returns a `Number`.
- `x1, x2,...` if `record_from_solution` returns a `Tuple`.
- the keys of the `NamedTuple` returned by `record_from_solution`.

The available symbols are provided by calling `propertynames(br.branch)`.

## Plotting bifurcation diagrams

To do this, you just need to call

```julia
plot(diagram)
```

where `diagram` is a branch computed after a call to `diagram = bifurcationdiagram(...)`. You can use the keywords provided by `Plots.jl` and the different backends. You can thus call `scatter(diagram)`. In addition to the options for plotting branches (see above), there are specific arguments available for bifurcation diagrams

- `code` specify the part of the bifurcation diagram to plot. For example `code = (1,1,)` plots the part after the first branch of the first branch of the root branch.
- `level = (-Inf, Inf)` restrict the branching level for plotting.

## Plotting without the plot recipe

What if you don't want to use Plots.jl? You can define your own plotting functions using the internal fields of `br` which is of type [`ContResult`](@ref). For example, in PyPlot, Gadfly, GR, etc., you can do the following to plot the branch (like the plot recipe `plot(br, vars = (:param, :x))`):

```julia
plot(br.branch.param, br.branch.x)
```

You can also have access to the stability of the points by using `br.stable`. More information concerning the fields can be found in [`ContResult`](@ref). For example, you can change the color depending on the stability:

```julia
col = [stb ? :green : :red for stb in br.stable]
plot(br.param, br.x, color=col)
```

You can also plot the spectrum at a specific continuation `step::Int` by calling

```julia
# get the eigenvalues
eigvals = br.eig[step].eigenvals

# plot them in the complex plane
scatter(real.(eigvals), imag.(eigvals))
```

## Standard plots using the Makie.jl [Experimental]

Plotting is also provided by calling recipes to `Makie.jl`. It means that to plot a branch `br`, you just need to call

```julia
#]add GLMakie # You need to install GLMakie.jl before your first time using it!
using GLMakie
BifurcationKit.plot(br)
```

The keyword arguments to `BifurcationKit.plot` are the same as decribed above in the page. You can also combine diagrams with `BifurcationKit.plot(br1, br2)` or use `BifurcationKit.plot!(ax, br)` to add a branch to an existing plot.

### Example

```julia
using Revise, GLMakie, BifurcationKit
Makie.inline!(true)
q = 1/0
k = 2
F(x, p) = (@. p + x - x^(k+1)/(k+1))
prob = BifurcationProblem(F, [0.8], 1., (@lens _); record_from_solution = (x,p) -> x[1])
opts = ContinuationPar(dsmax = 0.1, dsmin = 1e-3, ds = -0.001, p_min = -1., p_max = 1.)
br = continuation(prob, PALC(), opts)
BifurcationKit.plot(br)
```