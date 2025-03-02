# [🟢 2d Ginzburg-Landau equation (Shooting)](@id cglshoot)

```@contents
Pages = ["tutorialsCGLShoot.md"]
Depth = 3
```

In this tutorial, we re-visit the example [2d Ginzburg-Landau equation (finite differences, codim 2, Hopf aBS)](@ref cgl) using a Standard Simple Shooting method. In the tutorial [1d Brusselator (advanced user)](@ref), we used the implicit solver `Rodas4P` for the shooting. We will use the exponential-RK scheme `ETDRK2` ODE solver to compute the solution of cGL equations. This method is convenient for solving semilinear problems of the form

$$\dot x = Ax+g(x)$$

where $A$ is the infinitesimal generator of a $C_0$-semigroup. We use the same beginning as in [2d Ginzburg-Landau equation (finite differences, codim 2, Hopf aBS)](@ref cgl):

```julia
using Revise, DifferentialEquations
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
    return A, D2x
end
```

We then encode the PDE:

```julia
function NL!(f, u, p, t = 0.)
	(;r, μ, ν, c3, c5) = p
	n = div(length(u), 2)
	u1 = @view u[1:n]
	u2 = @view u[n+1:2n]

	ua = u1.^2 .+ u2.^2

	f1 = @view f[1:n]
	f2 = @view f[n+1:2n]

	@. f1 .= r * u1 - ν * u2 - ua * (c3 * u1 - μ * u2) - c5 * ua^2 * u1
	@. f2 .= r * u2 + ν * u1 - ua * (c3 * u2 + μ * u1) - c5 * ua^2 * u2

	return f
end

function Fcgl!(f, u, p, t = 0.)
	mul!(f, p.Δ, u)
	f .= f .+ NL(u, p)
end

NL(u, p) = NL!(similar(u), u, p)
Fcgl(u, p, t = 0.) = Fcgl!(similar(u), u, p, t)

function Jcgl(u, p, t = 0.)
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
```

with parameters

```julia
Nx = 41
Ny = 21
n = Nx*Ny
lx = pi
ly = pi/2

Δ = Laplacian2D(Nx, Ny, lx, ly)[1]
par_cgl = (r = 0.5, μ = 0.1, ν = 1.0, c3 = -1.0, c5 = 1.0, Δ = blockdiag(Δ, Δ))
sol0 = 0.1rand(2Nx, Ny)
sol0_f = vec(sol0)

prob = BK.BifurcationProblem(Fcgl, sol0_f, par_cgl, (@optic _.r); J = Jcgl)
```

and the ODE problem

```julia
f1 = DiffEqArrayOperator(par_cgl.Δ)
f2 = NL!
prob_sp = SplitODEProblem(f1, f2, sol0_f, (0.0, 120.0), (@set par_cgl.r = 1.2), dt = 0.1)
# we solve the PDE!!!
sol = @time DifferentialEquations.solve(prob_sp, ETDRK2(krylov=true); abstol=1e-14, reltol=1e-14)
```

## Automatic branch switching from the Hopf points

We show how to use automatic branch switching from the Hopf points computed in the previous section. To compute the periodic orbits, we use a Standard Shooting method.

We first recompute the Hopf points as in the previous tutorial:

```julia
eigls = EigArpack(1.0, :LM)
opt_newton = NewtonPar(tol = 1e-9, verbose = true, eigsolver = eigls, max_iterations = 20)
opts_br = ContinuationPar(dsmax = 0.02, ds = 0.01, p_max = 2., nev = 15, newton_options = (@set opt_newton.verbose = false))

br = @time continuation(prob, PALC(), opts_br, verbosity = 0)
```

We then compute the differentials of the vector field, this is needed by the branch switching method because it first computes the Hopf normal form. Thankfully, this is little work using Automatic Differentiation.

We define the linear solvers to be use by the (Matrix-Free) shooting method

```julia
ls = GMRESIterativeSolvers(reltol = 1e-4, maxiter = 50, verbose = false)
eig = EigKrylovKit(tol = 1e-7, x₀ = rand(2Nx*Ny), verbose = 2, dim = 40)
optn = NewtonPar(verbose = true, tol = 1e-9,  max_iterations = 25, linsolver = ls, eigsolver = eig)
opts_po_cont = ContinuationPar(dsmin = 0.001, dsmax = 0.02, ds= 0.01, p_max = 2.5, max_steps = 32, newton_options = optn, nev = 7, tol_stability = 1e-3, plot_every_step = 1)
```

as

```julia
Mt = 1 # number of time sections
br_po = continuation(
	# we want to compute the bifurcated branch from
	# the first Hopf point
	br, 1,
	# arguments for continuation
	opts_po_cont,
	# this is how to pass the method to compute the periodic
	# orbits. We shall use 1 section and the ODE solver ETDRK2
	ShootingProblem(Mt, prob_sp, ETDRK2(krylov = true); abstol = 1e-10, reltol = 1e-8, jacobian = BK.FiniteDifferencesMF()) ;
	# linear solver for bordered linear system
	# we combine the 2 solves. It is here faster than BorderingBLS()
	linear_algo = MatrixFreeBLS(@set ls.N = Mt*2n+2),
	# regular parameters for the continuation
	verbosity = 3, plot = true,
	# plotting of a section
	plot_solution = (x, p; k...) -> heatmap!(reshape(x[1:Nx*Ny], Nx, Ny); color=:viridis, k...),
	# print the Floquet exponent
	finalise_solution = (z, tau, step, contResult; k...) ->
		(Base.display(contResult.eig[end].eigenvals) ;true),
	normC = norminf)
```

![](cgl-sh-br.png)

## Manual branch switching from the Hopf points

The goal of this section is to show how to use the package in case automatic branch switching fails. This can happen for tedious PDEs and "one has to get his hands dirty".

We decide to use Standard Shooting and thus define a Shooting functional

```julia
probSh = ShootingProblem(
	# we pass the ODEProblem encoding the flow and the time stepper
	remake(prob_sp, p = (@set par_cgl.r = 1.2)), ETDRK2(krylov = true),

	# this is the phase condition
	[sol[:, end]];

	# parameter axis
	lens = (@optic _.r),

	# jacobian of the periodic orbit functional
	jacobian = BK.FiniteDifferencesMF(),

	# these are options passed to the ODE time stepper
	abstol = 1e-14, reltol = 1e-14)
```

We use the solution from the ODE solver as a starting guess for the shooting method.

```julia
# initial guess with period 6.9 using solution at time t = 116
initpo = vcat(sol(116.), 6.9) |> vec

# linear solver for shooting functional
ls = GMRESIterativeSolvers(reltol = 1e-4, N = 2Nx * Ny + 1, maxiter = 50, verbose = false)

# newton parameters
optn = NewtonPar(verbose = true, tol = 1e-9,  max_iterations = 20, linsolver = ls)

# continuation parameters
eig = EigKrylovKit(tol=1e-7, x₀ = rand(2Nx*Ny), verbose = 2, dim = 40)
opts_po_cont = ContinuationPar(dsmin = 0.001,
				dsmax = 0.01,
				ds = -0.01,
				p_max = 1.5,
				max_steps = 60,
				newton_options = (@set optn.eigsolver = eig),
				nev = 5,
				tol_stability = 1e-3,
				detect_bifurcation = 0)

br_po = @time continuation(probSh,
	initpo, PALC(), opts_po_cont;
	verbosity = 3, plot = true,
	linear_algo = MatrixFreeBLS(@set ls.N = probSh.M*2n+2),
	plot_solution = (x, p; kwargs...) -> heatmap!(reshape(x[1:Nx*Ny], Nx, Ny); color=:viridis, kwargs...),
	normC = norminf)
```
