# [🟠 2d Ginzburg-Landau equation (finite differences, codim 2, Hopf aBS)](@id cgl)

```@contents
Pages = ["tutorialsCGL.md"]
Depth = 3
```

> This example is also treated in the MATLAB library [pde2path](http://www.staff.uni-oldenburg.de/hannes.uecker/pde2path/).

We look at the Ginzburg-Landau equations in 2d. The code is very similar to the Brusselator example except that some special care has to be taken in order to cope with the "high" dimensionality of the problem.

Note that we try to be pedagogical here. Hence, we may write "bad" code that we improve later. Finally, we could use all sort of tricks to take advantage of the specificity of the problem. Rather, we stay quite close to the example in the MATLAB library [pde2path](http://www.staff.uni-oldenburg.de/hannes.uecker/pde2path/) (and discussed in **Hopf Bifurcation and Time Periodic Orbits with Pde2path – Algorithms and Applications.**, Uecker, Hannes, Communications in Computational Physics 25, no. 3 (2019)) for fair comparison.


!!! info "Goal"
    We do not use automatic branch switching here. The goal is to show how to use the internals of the package to squeeze most of the performances, use tailored options...


The equations are as follows

$$\partial_{t} u=\Delta u+(r+\mathrm{i} v) u-\left(c_{3}+\mathrm{i} \mu\right)|u|^{2} u-c_{5}|u|^{4} u+\gamma, \quad u=u(t, x) \in \mathbb{C},\quad \gamma\in\mathbb R$$

with Dirichlet boundary conditions. We discretize the square $\Omega = (0,L_x)\times(0,L_y)$ with $2N_xN_y$ points. We start by writing the Laplacian:

```@example CGL2d
using Revise, ForwardDiff
using BifurcationKit, LinearAlgebra, Plots, SparseArrays
const BK = BifurcationKit

function Laplacian2D(Nx, Ny, lx, ly)
    hx = 2lx/Nx
    hy = 2ly/Ny
    D2x = spdiagm(0 => -2ones(Nx), 1 => ones(Nx-1), -1 => ones(Nx-1) ) / hx^2
    D2y = spdiagm(0 => -2ones(Ny), 1 => ones(Ny-1), -1 => ones(Ny-1) ) / hy^2

    D2x[1,1] = -2/hx^2
    D2x[end,end] = -2/hx^2

    D2y[1,1] = -2/hy^2
    D2y[end,end] = -2/hy^2

    D2xsp = sparse(D2x)
    D2ysp = sparse(D2y)
    A = kron(sparse(I, Ny, Ny), D2xsp) + kron(D2ysp, sparse(I, Nx, Nx))
    return A
end
nothing #hide
```

It is then straightforward to write the vector field

```@example CGL2d
# this encodes the nonlinearity
function NL(u, p)
	(;r, μ, ν, c3, c5, γ) = p
	n = div(length(u), 2)
	u1 = @view u[1:n]
	u2 = @view u[n+1:2n]

	ua = u1.^2 .+ u2.^2

	f = similar(u, promote_type(eltype(u), typeof(r), typeof(γ)))
	f1 = @view f[1:n]
	f2 = @view f[n+1:2n]

	@. f1 .= r * u1 - ν * u2 - ua * (c3 * u1 - μ * u2) - c5 * ua^2 * u1 + γ
	@. f2 .= r * u2 + ν * u1 - ua * (c3 * u2 + μ * u1) - c5 * ua^2 * u2

	return f
end

function Fcgl(u, p)
	f = similar(u, promote_type(eltype(u), typeof(p.r), typeof(p.γ)))
	mul!(f, p.Δ, u)
	f .= f .+ NL(u, p)
end
nothing #hide
```

and its jacobian:

```@example CGL2d
function Jcgl(u, p)
	(;r, μ, ν, c3, c5, Δ) = p

	n = div(length(u), 2)
	u1 = @view u[1:n]
	u2 = @view u[n+1:2n]

	ua = u1.^2 .+ u2.^2

	f1u = zero(u1)
	f2u = zero(u1)
	f1v = zero(u1)
	f2v = zero(u1)

	@. f1u =  r - 2 * u1 * (c3 * u1 - μ * u2) - c3 * ua - 4 * c5 * ua * u1^2 - c5 * ua^2
	@. f1v = -ν - 2 * u2 * (c3 * u1 - μ * u2)  + μ * ua - 4 * c5 * ua * u1 * u2
	@. f2u =  ν - 2 * u1 * (c3 * u2 + μ * u1)  - μ * ua - 4 * c5 * ua * u1 * u2
	@. f2v =  r - 2 * u2 * (c3 * u2 + μ * u1) - c3 * ua - 4 * c5 * ua * u2 ^2 - c5 * ua^2

	jacdiag = vcat(f1u, f2v)

	Δ + spdiagm(0 => jacdiag, n => f1v, -n => f2u)
end
nothing #hide
```

We now define the parameters and the stationary solution:

```@example CGL2d
Nx = 41
Ny = 21
n = Nx * Ny
lx = pi
ly = pi/2

Δ = Laplacian2D(Nx, Ny, lx, ly)
par_cgl = (r = 0.5, μ = 0.1, ν = 1.0, c3 = -1.0, c5 = 1.0, Δ = blockdiag(Δ, Δ), γ = 0.)
sol0 = zeros(2Nx, Ny)

# we make a problem
prob = BK.BifurcationProblem(Fcgl, vec(sol0), par_cgl, (@optic _.r); J = Jcgl)
nothing #hide
```

and we continue it to find the Hopf bifurcation points. We use a Shift-Invert eigensolver.

```@example CGL2d
# Shift-Invert eigensolver
eigls = EigArpack(1.0, :LM) # shift = 1.0
opt_newton = NewtonPar(tol = 1e-10, eigsolver = eigls)
opts_br = ContinuationPar(dsmin = 0.001, dsmax = 0.005, ds = 0.001, p_max = 2., detect_bifurcation = 3, nev = 5, plot_every_step = 50, newton_options = opt_newton, max_steps = 1060)

br = continuation(prob, PALC(), opts_br)
```

which gives

```@example CGL2d
scene = plot(br, ylims=(-0.1,0.1))
```		

## Normal form computation

We compute the Hopf normal form of the first bifurcation point.

```@example CGL2d
hopfpt = get_normal_form(br, 1)
```

So the Hopf branch is subcritical.

## Codim 2 Hopf continuation

Having detected 2 hopf bifurcation points, we now continue them in the plane $(\gamma, r)$.

Before we start the codim 2 continuation, we tell `BifurcationKit.jl` to use the spectral information `start_with_eigen = true` because the left eigenvector of the Jacobian is simply not the conjugate of the right one.

```@example CGL2d
# we perform Hopf continuation of the first Hopf point in br
ind_hopf = 1
br_hopf = @time continuation(
	br, ind_hopf, (@optic _.γ),
	ContinuationPar(dsmin = 0.001, dsmax = 0.02, ds= 0.01, p_max = 6.5, p_min = -10., newton_options = opts_br.newton_options, detect_bifurcation = 1); plot = true,
	update_minaug_every_step = 1,
	normC = norminf,
	detect_codim2_bifurcation = 2,
	jacobian_ma = BK.MinAug(), # specific to high dimensions
	bdlinsolver = BorderingBLS(solver = DefaultLS(), check_precision = false))

plot(br_hopf, branchlabel = "Hopf curve", legend = :topleft)
```

with detailed information

```@example CGL2d
br_hopf
```

We can now construct the curve of Fold points branching off the Bogdanov-Takens bifurcation we just detected.

```@example CGL2d
# find the index of the BT point
indbt = findfirst(x -> x.type == :bt, br_hopf.specialpoint)
# branch from the BT point
brfold = continuation(br_hopf, indbt, setproperties(br_hopf.contparams; detect_bifurcation = 1, max_steps = 20);
	update_minaug_every_step = 1,
	detect_codim2_bifurcation = 2,
	callback_newton = BK.cbMaxNorm(1e5),
	bdlinsolver = BorderingBLS(solver = DefaultLS(), check_precision = false),
	jacobian_ma = BK.MinAug(), # !! specific to high dimensions
	bothside = true, normC = norminf)

plot(br_hopf, branchlabel = "Hopf"); plot!(brfold, legend = :topleft, branchlabel = "Fold")
```

## Periodic orbits continuation with stability
Having found two Hopf bifurcation points, we aim at computing the periodic orbits branching from them. Like for the Brusselator example, we need to find some educated guess for the periodic orbits in order to have a successful Newton call.

The following code is very close to the one explained in the tutorial [1d Brusselator (advanced user)](@ref) so we won't give too much details here.

We focus on the first Hopf bifurcation point. Note that, we do not improve the guess for the Hopf bifurcation point, *e.g.* by calling `newtonHopf`, as this is not really needed.

```julia
# index of the Hopf point we want to branch from
ind_hopf = 1

# number of time slices in the periodic orbit
M = 30

# periodic orbit initial guess from Hopf point
r_hopf, Th, orbitguess2, hopfpt, eigvec = guess_from_hopf(br, ind_hopf, opt_newton.eigsolver,
	# we pass the number of time slices M, the amplitude 22*sqrt(0.1) and phase
	M, 22*sqrt(0.1); phase = 0.25)

# flatten the initial guess
orbitguess_f2 = reduce(hcat, orbitguess2)
orbitguess_f = vcat(vec(orbitguess_f2), Th) |> vec
```

We create a problem to hold the functional and compute Periodic orbits based on Trapezoidal rule

```julia
poTrap = PeriodicOrbitTrapProblem(
# vector field and sparse Jacobian
	re_make(prob, params = (@set par_cgl.r = r_hopf - 0.01)),
# parameters for the phase condition
	real.(eigvec),
	hopfpt.u,
# number of time slices
	M,
# space dimension
	2n,
# jacobian of the periodic orbit functional
  jacobian = BK.MatrixFree())
```

We can use this (family) problem `poTrap` with `newton` on our periodic orbit guess to find a periodic orbit. Hence, one can be tempted to use


!!! danger "Don't run this!!"
    It uses too much memory

    ```julia
    opts_po_cont = ContinuationPar(dsmin = 0.0001, dsmax = 0.03, ds= 0.001, p_max = 2.5, max_steps = 250, plot_every_step = 3, newton_options = (@set opt_po.linsolver = DefaultLS()))
    br_po = @time continuation(poTrap, vec(sol0), PALC(), opts_po_cont)
    ```


**However, the linear system associated to the newton iterations will be solved by forming the sparse jacobian of size $(2N_xN_yM+1)^2$ and the use of `\` (based on LU decomposition). It takes way too much time and memory.**

Instead, we use a preconditioner. We build the jacobian once, compute its **incomplete LU decomposition** (ILU) and use it as a preconditioner.

```julia
using IncompleteLU

# Sparse matrix representation of the jacobian of the periodic orbit functional
Jpo = poTrap(Val(:JacFullSparse), orbitguess_f, @set par_cgl.r = r_hopf - 0.01)

# incomplete LU factorization with threshold
Precilu = @time ilu(Jpo, τ = 0.005);

# we define the linear solver with left preconditioner Precilu
ls = GMRESIterativeSolvers(verbose = false, reltol = 1e-3, N = size(Jpo,1), restart = 40, maxiter = 50, Pl = Precilu, log=true)

# we try the linear solver
ls(Jpo, rand(ls.N))
```

This converges in `7` iterations whereas, without the preconditioner, it does not converge after `100` iterations.

We set the parameters for the `newton` solve.

```julia
opt_po = @set opt_newton.verbose = true
outpo_f = @time newton(poTrap, orbitguess_f,  (@set opt_po.linsolver = ls); normN = norminf)
BK.converged(outpo_f) && printstyled(color=:red, "--> T = ", outpo_f.u[end], "\n")
BK.plot_periodic_potrap(outpo_f.u, M, Nx, Ny; ratio = 2);
```

which gives

```julia
┌─────────────────────────────────────────────────────┐
│ Newton step         residual      linear iterations │
├─────────────┬──────────────────────┬────────────────┤
│       0     │       6.5442e-03     │        0       │
│       1     │       1.4382e-03     │        7       │
│       2     │       3.7238e-04     │        8       │
│       3     │       6.4118e-05     │       10       │
│       4     │       4.2419e-06     │       10       │
│       5     │       5.6974e-08     │       11       │
│       6     │       3.1774e-10     │       12       │
│       7     │       3.1674e-13     │       14       │
└─────────────┴──────────────────────┴────────────────┘
  0.793448 seconds (143.31 k allocations: 1.242 GiB, 4.77% gc time)
--> T = 6.532023020978835, amplitude = 0.2684635643839235
```

and

![](cgl2d-po-newton.png)

At this point, we are still wasting a lot of resources, because the matrix-free version of the jacobian of the functional uses the jacobian of the vector field `x ->  Jcgl(x, p)`. Hence, it builds `M` sparse matrices for each evaluation!! Let us create a problem which is fully Matrix Free:

```julia
# computation of the first derivative using automatic differentiation
d1Fcgl(x, p, dx) = ForwardDiff.derivative(t -> Fcgl(x .+ t .* dx, p), 0.)

# linear solver for solving Jcgl*x = rhs. Needed for Floquet multipliers computation
ls0 = GMRESIterativeSolvers(N = 2Nx*Ny, reltol = 1e-9, Pl = lu(I + par_cgl.Δ))

poTrapMF = PeriodicOrbitTrapProblem(
# vector field and sparse Jacobian
	re_make(prob, params = (@set par_cgl.r = r_hopf - 0.01), J = (x, p) ->  (dx -> d1Fcgl(x, p, dx))),
# parameters for the phase condition
	real.(eigvec),
	hopfpt.u,
# number of time slices
	M,
# space dimension
	2n,
  ls0,
# jacobian of the periodic orbit functional
  jacobian = BK.MatrixFree())
```

We can now use newton

```julia
outpo_f = @time newton(poTrapMF, orbitguess_f, (@set opt_po.linsolver = ls); normN = norminf)
BK.converged(outpo_f) && printstyled(color=:red, "--> T = ", outpo_f.u[end], "\n")
```

which gives

```julia
┌─────────────────────────────────────────────────────┐
│ Newton step         residual      linear iterations │
├─────────────┬──────────────────────┬────────────────┤
│       0     │       6.5442e-03     │        0       │
│       1     │       1.4382e-03     │        7       │
│       2     │       3.7238e-04     │        8       │
│       3     │       6.4118e-05     │       10       │
│       4     │       4.2419e-06     │       10       │
│       5     │       5.6974e-08     │       11       │
│       6     │       3.1774e-10     │       12       │
│       7     │       3.1681e-13     │       14       │
└─────────────┴──────────────────────┴────────────────┘
  0.607167 seconds (46.11 k allocations: 461.511 MiB, 2.03% gc time)
```

The speedup will increase a lot for larger $N_x, N_y$. Also, for Floquet multipliers computation, the speedup will be substantial.

### Removing most allocations (Advanced and Experimental)

We show here how to remove most allocations and speed up the computations. This is an **experimental** feature as the Floquet multipliers computation is not yet readily available in this case. To this end, we rewrite the functional using *inplace* formulation and trying to avoid allocations. This can be done as follows:

```julia
# compute just the nonlinearity
function NL!(f, u, p, t = 0.)
	(;r, μ, ν, c3, c5) = p
	n = div(length(u), 2)
	u1v = @view u[1:n]
	u2v = @view u[n+1:2n]

	f1 = @view f[1:n]
	f2 = @view f[n+1:2n]

	@inbounds for ii = 1:n
		u1 = u1v[ii]
		u2 = u2v[ii]
		ua = u1^2+u2^2
		f1[ii] = r * u1 - ν * u2 - ua * (c3 * u1 - μ * u2) - c5 * ua^2 * u1
		f2[ii] = r * u2 + ν * u1 - ua * (c3 * u2 + μ * u1) - c5 * ua^2 * u2
	end
	return f
end

# derivative of the nonlinearity
function dNL!(f, u, p, du)
	(;r, μ, ν, c3, c5) = p
	n = div(length(u), 2)
	u1v = @view u[1:n]
	u2v = @view u[n+1:2n]

	du1v = @view du[1:n]
	du2v = @view du[n+1:2n]

	f1 = @view f[1:n]
	f2 = @view f[n+1:2n]

	@inbounds for ii = 1:n
		u1 = u1v[ii]
		u2 = u2v[ii]
		du1 = du1v[ii]
		du2 = du2v[ii]
		ua = u1^2+u2^2
		f1[ii] = (-5*c5*u1^4 + (-6*c5*u2^2 - 3*c3)*u1^2 + 2*μ*u1*u2 - c5*u2^4 - c3*u2^2 + r) * du1 +
		(-4*c5*u2*u1^3 + μ*u1^2 + (-4*c5*u2^3 - 2*c3*u2)*u1 + 3*u2^2*μ - ν) * du2

		f2[ii] = (-4*c5*u2*u1^3 - 3*μ*u1^2 + (-4*c5*u2^3 - 2*c3*u2)*u1 - u2^2*μ + ν) * du1 + (-c5*u1^4 + (-6*c5*u2^2 - c3)*u1^2 - 2*μ*u1*u2 - 5*c5*u2^4 - 3*c3*u2^2 + r) * du2
	end

	return f
end

# inplace vector field
function Fcgl!(f, u, p, t = 0.)
	NL!(f, u, p)
	mul!(f, p.Δ, u, 1., 1.)
end

# inplace derivative of the vector field
function dFcgl!(f, x, p, dx)
	dNL!(f, x, p, dx)
	mul!(f, p.Δ, dx, 1., 1.)
end

# we create a problem inplace
probInp = BifurcationProblem(Fcgl!, vec(sol0), (@set par_cgl.r = r_hopf - 0.01), (@optic _.r); J = dFcgl!, inplace = true)
```

We can now define an inplace functional

```julia
ls0 = GMRESIterativeSolvers(N = 2Nx*Ny, reltol = 1e-9)#, Pl = lu(I + par_cgl.Δ))

poTrapMFi = PeriodicOrbitTrapProblem(
# vector field and sparse Jacobian
	probInp,
# parameters for the phase condition
	real.(eigvec),
	hopfpt.u,
# number of time slices
	M,
# space dimension
	2n,
  ls0,
# jacobian of the periodic orbit functional
  jacobian = BK.MatrixFree())
```
and run the `newton` method:

```julia
outpo_f = @time newton(poTrapMFi, orbitguess_f, (@set opt_po.linsolver = ls); normN = norminf)
```
It gives

```julia
┌─────────────────────────────────────────────────────┐
│ Newton step         residual      linear iterations │
├─────────────┬──────────────────────┬────────────────┤
│       0     │       6.5442e-03     │        0       │
│       1     │       1.4382e-03     │        7       │
│       2     │       3.7238e-04     │        8       │
│       3     │       6.4118e-05     │       10       │
│       4     │       4.2419e-06     │       10       │
│       5     │       5.6974e-08     │       11       │
│       6     │       3.1774e-10     │       12       │
│       7     │       3.1674e-13     │       14       │
└─────────────┴──────────────────────┴────────────────┘
  0.583849 seconds (5.55 k allocations: 151.581 MiB)
```

Notice the small speed boost but the reduced allocations. At this stage, further improvements could target the use of `BlockBandedMatrices.jl` for the Laplacian operator, etc.


### Other linear formulation

We could use another way to "invert" jacobian of the functional based on bordered techniques. We try to use an ILU preconditioner on the cyclic matrix $J_c$ (see [Periodic orbits based on Trapezoidal rule](@ref)) which has a smaller memory footprint:

```julia
Jpo2 = poTrap(Val(:JacCyclicSparse), orbitguess_f, @set par_cgl.r = r_hopf - 0.1)
Precilu = @time ilu(Jpo2, τ = 0.005)
ls2 = GMRESIterativeSolvers(verbose = false, reltol = 1e-3, N = size(Jpo2, 1), restart = 30, maxiter = 50, Pl = Precilu, log=true)

opt_po = @set opt_newton.verbose = true
outpo_f = @time newton((@set poTrapMF.jacobian = BK.BorderedMatrixFree()),
	orbitguess_f,
	(@set opt_po.linsolver = ls2), 
	normN = norminf)
```

but it gives:

```julia
┌─────────────────────────────────────────────────────┐
│ Newton step         residual      linear iterations │
├─────────────┬──────────────────────┬────────────────┤
│       0     │       6.5442e-03     │        0       │
│       1     │       1.4262e-03     │       26       │
│       2     │       3.6549e-04     │       28       │
│       3     │       6.3752e-05     │       32       │
│       4     │       4.2248e-06     │       37       │
│       5     │       2.8050e-08     │       41       │
│       6     │       5.5126e-11     │       48       │
└─────────────┴──────────────────────┴────────────────┘
  1.581654 seconds (84.17 k allocations: 947.476 MiB, 2.92% gc time)
```

**Hence, it seems better to use the previous preconditioner.**

## Continuation of periodic solutions

We can now perform continuation of the newly found periodic orbit and compute the Floquet multipliers using Matrix-Free methods.

```julia
# set the eigensolver for the computation of the Floquet multipliers
opt_po = @set opt_po.eigsolver = EigKrylovKit(tol = 1e-3, x₀ = rand(2n), verbose = 2, dim = 25)

# parameters for the continuation
opts_po_cont = ContinuationPar(dsmin = 0.0001, dsmax = 0.02, ds = 0.001, p_max = 2.2, max_steps = 250, plot_every_step = 3, newton_options = (@set opt_po.linsolver = ls),
	nev = 5, tol_stability = 1e-7, detect_bifurcation = 3)

br_po = @time continuation(poTrapMF, outpo_f.u, PALC(),	opts_po_cont;
	verbosity = 2,	plot = true,
	linear_algo = BorderingBLS(solver = ls, check_precision = false),
	plot_solution = (x, p; kwargs...) -> BK.plot_periodic_potrap(x, M, Nx, Ny; ratio = 2, kwargs...),
	record_from_solution = (u, p; k...) -> BK.getamplitude(poTrapMF, u, par_cgl; ratio = 2), normC = norminf)
```

This gives the following bifurcation diagram:

![](cgl2d-po-cont.png)

!!! tip "Improved performances"
    Although it would be "cheating" for fair comparisons with existing packages, there is a trick to compute the bifurcation diagram without using preconditionners. We will not detail it here but it allows to handle the case `Nx = 200; Ny = 110; M = 30` and above.

We did not change the preconditioner in the previous example as it does not seem needed. Let us show how to do this nevertheless:

```julia
# callback which will be sent to newton.
# `iteration` in the arguments refers to newton iterations
function callbackPO(state; linsolver = ls, prob = poTrap, p = par_cgl, kwargs...)
	@show ls.N keys(kwargs)
	# we update the preconditioner every 10 continuation steps
	if mod(kwargs[:iterationC], 10) == 9 && state.step == 1
		@info "update Preconditioner"
		Jpo = poTrap(Val(:JacCyclicSparse), state.x, (@set p.r = state.p))
		Precilu = @time ilu(Jpo, τ = 0.003)
		ls.Pl = Precilu
	end
	true
end

br_po2 = @time continuation(poTrapMF, outpo_f.u, PALC(), opts_po_cont;
	verbosity = 2,	plot = true,
	callback_newton = callbackPO,
	linear_algo = BorderingBLS(solver = ls, check_precision = false),
	plot_solution = (x, p; kwargs...) -> BK.plot_periodic_potrap(x, M, Nx, Ny; ratio = 2, kwargs...),
	normC = norminf)
```

## Continuation of Fold of periodic orbits

We continue the Fold point of the first branch of the previous bifurcation diagram in the parameter plane $(r, c_5)$. To this end, we need to be able to compute the Hessian of the periodic orbit functional. This is not yet readily available so we turn to automatic differentiation.

```julia
using ForwardDiff

# computation of the second derivative of a function f
function d2Fcglpb(f, x, dx1, dx2)
   return ForwardDiff.derivative(t2 -> ForwardDiff.derivative( t1-> f(x .+ t1 .* dx1 .+ t2 .* dx2,), 0.), 0.)
end
```

We select the Fold point from the branch `br_po` and redefine our linear solver to get the ILU preconditioner tuned close to the Fold point.

```julia
indfold = 1
foldpt = foldpoint(br_po, indfold)

Jpo = poTrap(Val(:JacFullSparse), orbitguess_f, (@set par_cgl.r = r_hopf - 0.1))
Precilu = @time ilu(Jpo, τ = 0.005);
ls = GMRESIterativeSolvers(verbose = false, reltol = 1e-4, N = size(Jpo, 1), restart = 40, maxiter = 60, Pl = Precilu);
```

We can then use our functional to call `newtonFold` unlike for a regular function (see Tutorial 1). Indeed, we specify the change the parameters too much to rely on a generic algorithm.

```julia
# define a bifurcation problem for the fold
# this will be made automatic in the future
probFold = BifurcationProblem((x, p) -> poTrap(x, p), foldpt, getparams(br), getlens(br);
			J = (x, p) -> poTrap(Val(:JacFullSparse), x, p),
			d2F = (x, p, dx1, dx2) -> d2Fcglpb(z -> poTrap(z, p), x, dx1, dx2))

outfold = @time BK.newton_fold(
  br_po, indfold; #index of the fold point
  prob = probFold,
	# we change the linear solver for the one we
	# defined above
	options = (@set opt_po.linsolver = ls),
	bdlinsolver = BorderingBLS(solver = ls, check_precision = false))
BK.converged(outfold) && printstyled(color=:red, "--> We found a Fold Point at α = ", outfold.u.p," from ", br_po.specialpoint[indfold].param,"\n")
```

and this gives

```julia
┌─────────────────────────────────────────────────────┐
│ Newton step         residual     linear iterations  │
├─────────────┬──────────────────────┬────────────────┤
│       0     │       4.5937e-01     │        0       │
│       1     │       5.6013e-01     │       20       │
│       2     │       3.1385e-02     │       23       │
│       3     │       6.0620e-05     │       29       │
│       4     │       2.7839e-08     │       39       │
│       5     │       8.1593e-12     │       45       │
└─────────────┴──────-───────────────┴────────────────┘
 27.289005 seconds (1.07 M allocations: 24.444 GiB, 10.12% gc time)
--> We found a Fold Point at α = 0.9470569704262517 from 0.9481896723164748
```

Finally, one can perform continuation of the Fold bifurcation point as follows

```julia
optcontfold = ContinuationPar(dsmin = 0.001, dsmax = 0.05, ds= 0.01, p_max = 40.1, p_min = -10., newton_options = (@set opt_po.linsolver = ls), max_steps = 20, detect_bifurcation = 0)

outfoldco = BK.continuation_fold(probFold,
	br_po, indfold, (@optic _.c5),
	optcontfold;
	jacobian_ma = BK.MinAug(),
	bdlinsolver = BorderingBLS(solver = ls, check_precision = false),
	plot = true, verbosity = 2)
```

which yields:

![](cgl2d-po-foldcont.png)

There is still room for a lot of improvements here. Basically, the idea would be to use full Matrix-Free the jacobian functional and its transpose.

## Continuation of periodic orbits on the GPU (Advanced)

!!! tip ""
    This is a very neat example **all done** on the GPU using the following ingredients: Matrix-Free computation of periodic orbits using preconditioners.

We now take advantage of the computing power of GPUs. The section is run on an NVIDIA Tesla V100. Given the small number of unknowns, we can (only) expect significant speedup in the application of the **big** preconditioner.

> Note that we use the parameters `Nx = 82; Ny = 42; M=30`.

```julia
# computation of the first derivative
d1Fcgl(x, p, dx) = ForwardDiff.derivative(t -> Fcgl(x .+ t .* dx, p), 0.)

d1NL(x, p, dx) = ForwardDiff.derivative(t -> NL(x .+ t .* dx, p), 0.)

function dFcgl(x, p, dx)
	f = similar(dx)
	mul!(f, p.Δ, dx)
	nl = d1NL(x, p, dx)
	f .= f .+ nl
end
```

We first load `CuArrays`

```julia
using CUDA
CUDA.allowscalar(false)
import LinearAlgebra: mul!, axpby!
mul!(x::CuArray, y::CuArray, α::T) where {T <: Number} = (x .= α .* y)
mul!(x::CuArray, α::T, y::CuArray) where {T <: Number} = (x .= α .* y)
axpby!(a::T, X::CuArray, b::T, Y::CuArray) where {T <: Number} = (Y .= a .* X .+ b .* Y)
```

and update the parameters

```julia
par_cgl_gpu = @set par_cgl.Δ = CUDA.CUSPARSE.CuSparseMatrixCSC(par_cgl.Δ);
```

Then, we precompute the preconditioner on the CPU:

```julia
Jpo = poTrap(Val(:JacFullSparse), orbitguess_f, @set par_cgl.r = r_hopf - 0.01)
Precilu = @time ilu(Jpo, τ = 0.003)
```

To invert `Precilu` on the GPU, we need to define a few functions which are not in `CuArrays` and which are related to LU decomposition:

```julia
struct LUperso
	L
	Ut	# transpose of U in LU decomposition
end

import Base: ldiv!
function LinearAlgebra.ldiv!(_lu::LUperso, rhs::CUDA.CuArray)
	_x = UpperTriangular(_lu.Ut) \ (LowerTriangular(_lu.L) \ rhs)
	rhs .= vec(_x)
	CUDA.unsafe_free!(_x)
	rhs
end
```

Finally, for the methods in `PeriodicOrbitTrapProblem` to work, we need to redefine the following method. Indeed, we disable the use of scalar on the GPU to increase the speed.

```julia
import BifurcationKit: extractPeriodFDTrap
extractPeriodFDTrap(x::CuArray) = x[end:end]
```

We can now define our functional:

```julia
# we create a problem for the vector field on the GUP
probGPU = BifurcationProblem(Fcgl,
    vec(CuArrays(zeros((x, p) ->  (dx -> dFcgl(x, p, dx))))),
    (@set par_cgl_gpu.r = r_hopf - 0.01),
    (@optic _.r);
    J = (x, p) ->  (dx -> dFcgl(x, p, dx)))

# matrix-free problem on the gpu
ls0gpu = GMRESKrylovKit(rtol = 1e-9)
poTrapMFGPU = PeriodicOrbitTrapProblem(
	probGPU,
	CuArray(real.(eigvec)),
	CuArray(hopfpt.u),
	M, 2n, ls0gpu;
  jacobian = BK.MatrixFree(),
	ongpu = true) # this is required to alter the way the constraint is handled
```

Let us have a look at the linear solvers and compare the speed on CPU and GPU:

```julia
ls = GMRESKrylovKit(verbose = 2, Pl = Precilu, rtol = 1e-3, dim  = 20)
   # runs in 	2.990495 seconds (785 allocations: 31.564 MiB, 0.98% gc time)
	outh, = @time ls((Jpo), orbitguess_f)

Precilu_gpu = LUperso(LowerTriangular(CUDA.CUSPARSE.CuSparseMatrixCSR(I+Precilu.L)), UpperTriangular(CUDA.CUSPARSE.CuSparseMatrixCSR(sparse(Precilu.U'))));
lsgpu = GMRESKrylovKit(verbose = 2, Pl = Precilu_gpu, rtol = 1e-3, dim  = 20)
	Jpo_gpu = CUDA.CUSPARSE.CuSparseMatrixCSR(Jpo);
	orbitguess_cu = CuArray(orbitguess_f)
	# runs in 1.751230 seconds (6.54 k allocations: 188.500 KiB, 0.43% gc time)
	outd, = @time lsgpu(Jpo_gpu, orbitguess_cu)
```	 

So we can expect a pretty descent x2 speed up in computing the periodic orbits. We can thus call newton:

```julia
opt_po = @set opt_newton.verbose = true
outpo_f = @time newton(poTrapMFGPU, orbitguess_cu,
		(@set opt_po.linsolver = lsgpu);
		normN = x->maximum(abs.(x)))
```
The computing time is `6.914367 seconds (2.94 M allocations: 130.348 MiB, 1.10% gc time)`. The same computation on the CPU, runs in `13.972836 seconds (551.41 k allocations: 1.300 GiB, 1.05% gc time)`.

You can also perform continuation, here is a simple example:

```julia
opts_po_cont = ContinuationPar(dsmin = 0.0001, dsmax = 0.02, ds= 0.001, p_max = 2.2, max_steps = 250, plot_every_step = 3, newton_options = (@set opt_po.linsolver = lsgpu), detect_bifurcation = 0)
br_po = @time continuation(poTrapMFGPU, orbitguess_cu, PALC(),
   opts_po_cont;
   verbosity = 2,
   record_from_solution = (u, p; k...) -> getAmplitude(poTrapMFGPU, u, par_cgl_gpu), normC = x->maximum(abs.(x)))
```

!!! info "Preconditioner update"
    For now, the preconditioner has been precomputed on the CPU which forbids its (efficient) update during continuation of a branch of periodic orbits. This could be improved using `ilu0!` and friends in `CuArrays`.
