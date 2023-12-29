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
using BifurcationKit, Setfield, Plots

# vector field of transcritical bifurcation
F(x, p) = [x[1] * (p.μ - x[1])]

# parameters of the vector field
par = (μ = -0.2, )

# problem (automatic differentiation)
prob = BifurcationProblem(F, [0.], par, (@lens _.μ); record_from_solution = (x, p) -> x[1])

# compute branch of trivial equilibria and detect a bifurcation point
opts_br = ContinuationPar()
br = continuation(prob, PALC(), opts_br)
	
# perform branch switching on both sides of the bifurcation point
br1 = continuation(br, 1, setproperties(opts_br; max_steps = 14), bothside = true )

scene = plot(br, br1; branchlabel = ["br", "br1"], legend = :topleft)
```

### Simple example: pitchfork bifurcation

```@example TUT1b_ABS_EQ_EQ
using BifurcationKit, Setfield, Plots

# vector field of pitchfork bifurcation
F(x, p) = [x[1] * (p.μ - x[1]^2)]

# parameters of the vector field
par = (μ = -0.2, )

# problem (automatic differentiation)
prob = BifurcationProblem(F, [0.], par, (@lens _.μ); record_from_solution = (x, p) -> x[1])

# compute branch of trivial equilibria and 
# detect a bifurcation point with improved precision
opts_br = ContinuationPar(n_inversion = 6)
br = continuation(prob, PALC(), opts_br)
	
# perform branch switching on both sides of the bifurcation point
br1 = continuation(br, 1, setproperties(opts_br; max_steps = 14), bothside = true )

scene = plot(br, br1; branchlabel = ["br", "br1"], legend = :topleft)
```

### Algorithm
- the normal form is computed and non-trivial zeros are used to produce guesses for points on the bifurcated branch.


## [From non simple branch point to equilibria](@id abs-simple-eq)

We provide an automatic branch switching method in this case. The underlying method is to first compute the reduced equation (see [Non-simple branch point](@ref)) and its zeros. These zeros are then seeded as initial guess for [`continuation`](@ref). Hence, you can perform automatic branch switching by calling `continuation` with the following options:

```julia
continuation(br::ContResult, 
	ind_bif::Int,
	optionsCont::ContinuationPar;
	kwargs...)
```

An example of use is provided in [2d Bratu–Gelfand problem](@ref gelfand).	