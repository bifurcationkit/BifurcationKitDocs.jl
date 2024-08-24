# From equilibria to equilibria

```@contents
Pages = ["abs-from-eq.md"]
Depth = 3
```

## From simple branch point to equilibria

> See [Branch switching (branch point)](@ref) for the precise method definition


You can perform automatic branch switching by calling `continuation` with the following options:

```julia
continuation(br::ContResult, ind_bif::Int, optionsCont::ContinuationPar; kwargs...)
```

where `br` is a branch computed after a call to `continuation` with detection of bifurcation points enabled. This call computes the branch bifurcating from the `ind_bif `th bifurcation point in `br`. An example of use is provided in [2d Bratu–Gelfand problem](@ref gelfand).

> See [Branch switching (branch point)](@ref) precise method definition

### Simple example: transcritical bifurcation

```@example TUT1_ABS_EQ_EQ
using BifurcationKit, Plots

# vector field of transcritical bifurcation
F(x, p) = [x[1] * (p.μ - x[1])]

# parameters of the vector field
par = (μ = -0.2, )

# problem (automatic differentiation)
prob = BifurcationProblem(F, [0.], par, (@optic _.μ); record_from_solution = (x, p) -> x[1])

# compute branch of trivial equilibria and detect a bifurcation point
br = continuation(prob, PALC(), ContinuationPar())
	
# perform branch switching on both sides of the bifurcation point
br1 = continuation(br, 1; bothside = true )

scene = plot(br, br1; branchlabel = ["br", "br1"], legend = :topleft)
```

### Simple example: pitchfork bifurcation

```@example TUT1b_ABS_EQ_EQ
using BifurcationKit, Plots

# vector field of pitchfork bifurcation
F(x, p) = [x[1] * (p.μ - x[1]^2)]

# parameters of the vector field
par = (μ = -0.2, )

# problem (automatic differentiation)
prob = BifurcationProblem(F, [0.], par, (@optic _.μ); record_from_solution = (x, p) -> x[1])

# compute branch of trivial equilibria and 
# detect a bifurcation point with improved precision
br = continuation(prob, PALC(), ContinuationPar(n_inversion = 6))
	
# perform branch switching on both sides of the bifurcation point
br1 = continuation(br, 1; bothside = true )

scene = plot(br, br1; branchlabel = ["br", "br1"], legend = :topleft)
```

### Algorithm
- the normal form is computed and non-trivial zeros are used to produce guesses for points on the bifurcated branch.


## [From non simple branch point to equilibria](@id abs-nonsimple-eq)

We provide an automatic branch switching method in this case. The underlying method is to first compute the reduced equation (see [Non-simple branch point](@ref)) and its zeros. These zeros are then seeded as initial guess for [`continuation`](@ref). Hence, you can perform automatic branch switching by calling `continuation` with the following options:

```julia
continuation(br::ContResult, 
	ind_bif::Int,
	optionsCont::ContinuationPar;
	kwargs...)
```

### Examples

An example of use is provided in [2d Bratu–Gelfand problem](@ref gelfand). A much simpler example is given now. It is a bit artificial because the vector field is its own normal form at the bifurcation point located at 0.

```@example TUT2_ABS_EQ_EQ
using BifurcationKit, Plots

function FbpD6(x, p)
    return [p.μ * x[1] + (p.a * x[2] * x[3] - p.b * x[1]^3 - p.c*(x[2]^2 + x[3]^2) * x[1]),
           p.μ * x[2] + (p.a * x[1] * x[3] - p.b * x[2]^3 - p.c*(x[3]^2 + x[1]^2) * x[2]),
           p.μ * x[3] + (p.a * x[1] * x[2] - p.b * x[3]^3 - p.c*(x[2]^2 + x[1]^2) * x[3])]
end

# model parameters
pard6 = (μ = -0.2, a = 0.3, b = 1.5, c = 2.9)

# problem
prob = BifurcationProblem(FbpD6, zeros(3), pard6, (@optic _.μ);
	record_from_solution = (x, p) -> (n = norminf(x)))

# continuation options
opts_br = ContinuationPar(dsmin = 0.001, dsmax = 0.02, ds = 0.01, 
	# parameter interval
	p_max = 0.4, p_min = -0.2, 
	nev = 3, 
	newton_options = NewtonPar(tol = 1e-10, max_iterations = 20), 
	max_steps = 1000, 
	n_inversion = 6)

br = continuation(prob, PALC(), opts_br)
```

You can now branch from the `nd` point

```@example TUT2_ABS_EQ_EQ
br2 = continuation(br, 1, opts_br; δp = 0.02)

plot(br, br2...)
```

## Assisted branching from non-simple bifurcation point

It may happen that the general procedure fails. We thus expose the procedure `multicontinuation` in order to let the user tune it to its need.

The first step is to compute the reduced equation, say of the first bifurcation point in `br`.

```@example TUT2_ABS_EQ_EQ
bp = get_normal_form(br, 1)
```

Next, we want to find the zeros of the reduced equation. This is usually achieved by calling the predictor

```@example TUT2_ABS_EQ_EQ
δp = 0.005
pred = predictor(bp, δp)
```

which returns zeros of `bp` before and after the bifurcation point. You could also use your preferred procedure from `Roots.jl` (or other) to find the zeros of the polynomials `bp(Val(:reducedForm), z, p)`.

We can use these zeros to form guesses to apply Newton for the full functional:

```@example TUT2_ABS_EQ_EQ
pts = BifurcationKit.get_first_points_on_branch(br, bp, pred, opts_br; δp)
```

We can then use this to continue the different branches

```@example TUT2_ABS_EQ_EQ
brbp = BifurcationKit.multicontinuation(br, bp, pts.before, pts.after, opts_br)

plot(br, brbp...)
```

Note that you can chose another predictor which uses all vertices of the cube as initial guesses

```@example TUT2_ABS_EQ_EQ
pred = predictor(bp, Val(:exhaustive), δp)
pts = BifurcationKit.get_first_points_on_branch(br, bp, pred, opts_br; δp)
```

```@example TUT2_ABS_EQ_EQ
brbp = BifurcationKit.multicontinuation(br, bp, pts.before, pts.after, opts_br)

plot(br, brbp...)
```

## predictors 

```@docs
BifurcationKit.predictor(bp::BifurcationKit.NdBranchPoint, δp::T; k...) where T
```

```@docs
BifurcationKit.predictor(bp::BifurcationKit.NdBranchPoint, algo::Val{:exhaustive}, δp::T;k...) where T
```