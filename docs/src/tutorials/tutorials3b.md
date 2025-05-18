# [🟠 1d Brusselator](@id bruss)

```@contents
Pages = ["tutorials3b.md"]
Depth = 3
```

!!! info "Goal"
    The goal of this tutorial is to show similar computations as in the previous tutorial but without using the automatic branch switching tools. This is for the experienced used who wants to dive more in the internals of the package.    

We look at the Brusselator in 1d (see [^Lust]). The equations are:

$$\begin{aligned} \frac { \partial X } { \partial t } & = \frac { D _ { 1 } } { l ^ { 2 } } \frac { \partial ^ { 2 } X } { \partial z ^ { 2 } } + X ^ { 2 } Y - ( β + 1 ) X + α \\ \frac { \partial Y } { \partial t } & = \frac { D _ { 2 } } { l ^ { 2 } } \frac { \partial ^ { 2 } Y } { \partial z ^ { 2 } } + β X - X ^ { 2 } Y \end{aligned}$$

with Dirichlet boundary conditions

$$\begin{array} { l } { X ( t , z = 0 ) = X ( t , z = 1 ) = α } \\ { Y ( t , z = 0 ) = Y ( t , z = 1 ) = β / α. } \end{array}$$

These equations have been introduced to reproduce an oscillating chemical reaction. There is an obvious equilibrium $(α, β / α)$. Here, we consider bifurcations with respect to the parameter $l$.

We start by writing the PDE

```@example TUTBRUmanual
using Revise
using BifurcationKit, LinearAlgebra, Plots, SparseArrays
const BK = BifurcationKit

f1(u, v) = u * u * v

function Fbru!(f, x, p, t = 0)
	(;α, β, D1, D2, l) = p
	n = div(length(x), 2)
	h2 = 1.0 / n^2
	c1 = D1 / l^2 / h2
	c2 = D2 / l^2 / h2

	u = @view x[1:n]
	v = @view x[n+1:2n]

	# Dirichlet boundary conditions
	f[1]   = c1 * (α	  - 2u[1] + u[2] ) + α - (β + 1) * u[1] + f1(u[1], v[1])
	f[end] = c2 * (v[n-1] - 2v[n] + β / α)			 + β * u[n] - f1(u[n], v[n])

	f[n]   = c1 * (u[n-1] - 2u[n] +  α   ) + α - (β + 1) * u[n] + f1(u[n], v[n])
	f[n+1] = c2 * (β / α  - 2v[1] + v[2])			 + β * u[1] - f1(u[1], v[1])

	for i=2:n-1
		  f[i] = c1 * (u[i-1] - 2u[i] + u[i+1]) + α - (β + 1) * u[i] + f1(u[i], v[i])
		f[n+i] = c2 * (v[i-1] - 2v[i] + v[i+1])			  + β * u[i] - f1(u[i], v[i])
	end
	return f
end

Fbru(x, p, t = 0) = Fbru!(similar(x), x, p, t)
nothing #hide
```

We use a sparse representation of the Jacobian:

```@example TUTBRUmanual
function Jbru_sp(x, p)
	(;α, β, D1, D2, l) = p
	# compute the Jacobian using a sparse representation
	n = div(length(x), 2)
	𝒯 = eltype(x)
	h = 1.0 / n; h2 = h*h

	c1 = D1 / p.l^2 / h2
	c2 = D2 / p.l^2 / h2

	u = @view x[1:n]
	v = @view x[n+1:2n]

	diag   = zeros(𝒯, 2n)
	diagp1 = zeros(𝒯, 2n-1)
	diagm1 = zeros(𝒯, 2n-1)

	diagpn = zeros(𝒯, n)
	diagmn = zeros(𝒯, n)

	@. diagmn = β - 2 * u * v
	@. diagm1[1:n-1] = c1
	@. diagm1[n+1:end] = c2

	@. diag[1:n]    = -2c1 - (β + 1) + 2 * u * v
	@. diag[n+1:2n] = -2c2 - u * u

	@. diagp1[1:n-1] = c1
	@. diagp1[n+1:end] = c2

	@. diagpn = u * u
	return spdiagm(0 => diag, 1 => diagp1, -1 => diagm1, n => diagpn, -n => diagmn)
end

nothing #hide
```

We shall now compute the equilibria and their stability.

```@example TUTBRUmanual
n = 500

# parameters of the Brusselator model and guess for the stationary solution
par_bru = (α = 2., β = 5.45, D1 = 0.008, D2 = 0.004, l = 0.3)
sol0 = vcat(par_bru.α * ones(n), par_bru.β / par_bru.α * ones(n))

# bifurcation problem
probBif = BK.BifurcationProblem(Fbru!, sol0, par_bru, (@optic _.l);
  J = Jbru_sp,
  plot_solution = (x, p; kwargs...) -> (plotsol(x; label="", kwargs... )),
  record_from_solution = (x, p; k...) -> x[div(n,2)])
nothing #hide
```

For the eigensolver, we use a Shift-Invert algorithm (see [Eigen solvers (Eig)](@ref))

```@example TUTBRUmanual
eigls = EigArpack(1.1, :LM)
nothing #hide
```

We continue the trivial equilibrium to find the Hopf points

```@example TUTBRUmanual
opt_newton = NewtonPar(eigsolver = eigls, tol = 1e-9)
opts_br_eq = ContinuationPar(dsmin = 0.001, dsmax = 0.1,
	p_max = 1.9, nev = 21,
	newton_options = opt_newton, max_steps = 1000,
	# specific options for precise localization of Hopf points
	n_inversion = 6)

br = continuation(probBif, PALC(),opts_br_eq, normC = norminf)
```

We obtain the following bifurcation diagram with 3 Hopf bifurcation points

```@example TUTBRUmanual
scene = plot(br)
```

## Normal form computation

We can compute the normal form of the Hopf points as follows

```@example TUTBRUmanual
hopfpt = get_normal_form(br, 1)
```

## Continuation of Hopf points

We use the bifurcation points guesses located in `br.specialpoint` to turn them into precise bifurcation points. For the second one, we have

```@example TUTBRUmanual
# index of the Hopf point in br.specialpoint
ind_hopf = 2

# newton iterations to compute the Hopf point
hopfpoint = newton(br, ind_hopf; normN = norminf)
BK.converged(hopfpoint) && printstyled(color=:red, "--> We found a Hopf Point at l = ", hopfpoint.u.p[1], ", ω = ", hopfpoint.u.p[2], ", from l = ", br.specialpoint[ind_hopf].param, "\n")
```

We now perform a Hopf continuation with respect to the parameters `l, β`

!!! tip "Tip"
    You don't need to call `newton` first in order to use `continuation`.

```@example TUTBRUmanual
optcdim2 = ContinuationPar(dsmin = 0.001, dsmax = 0.05, ds= 0.01, p_max = 6.5, p_min = 0.0, newton_options = opt_newton, detect_bifurcation = 0)
br_hopf = continuation(br, ind_hopf, (@optic _.β), optcdim2, normC = norminf, jacobian_ma = BK.MinAug())
scene = plot(br_hopf)
```

## Continuation of periodic orbits (Finite differences)

Here, we perform continuation of periodic orbits branching from the Hopf bifurcation points.We need an educated guess for the periodic orbit which is given by `guess_from_hopf`:

```@example TUTBRUmanual
# number of time slices
M = 51

l_hopf, Th, orbitguess2, hopfpt, vec_hopf = BK.guess_from_hopf(br, ind_hopf,
	opts_br_eq.newton_options.eigsolver,
	M, 2.7; phase = 0.25)

nothing #hide
```
We wish to make two remarks at this point. The first is that an initial guess is composed of a space time solution and of the guess for the period `Th` of the solution. Note that the argument `2.7` is a guess for the amplitude of the orbit.

```@example TUTBRUmanual
# orbit initial guess from guess_from_hopf, is not a vector, so we reshape it
orbitguess_f2 = reduce(vcat, orbitguess2)
orbitguess_f = vcat(vec(orbitguess_f2), Th) |> vec

nothing #hide
```

The second remark concerns the phase `0.25` written above. To account for the additional unknown (*i.e.* the period), periodic orbit localization using Finite Differences requires an additional constraint (see [Periodic orbits based on Trapezoidal rule](@ref) for more details). In the present case, this constraint is

$$< u(0) - u_{hopf}, \phi> = 0$$

where `u_{hopf}` is the equilibrium at the Hopf bifurcation and $\phi$ is `real.(vec_hopf)` where `vec_hopf` is the eigenvector. This is akin to a Poincaré section. We do not put any constraint on $u(t)$ albeit this is possible (see [Periodic orbits based on Trapezoidal rule](@ref).

The phase of the periodic orbit is set so that the above constraint is satisfied. We shall now use Newton iterations to find a periodic orbit.

Given our initial guess, we create a (family of) problem which encodes the functional associated to finding Periodic orbits based on Trapezoidal rule (see [Periodic orbits based on Trapezoidal rule](@ref) for more information):

```@example TUTBRUmanual
poTrap = PeriodicOrbitTrapProblem(
	probBif,    		   # pass the bifurcation problem
	real.(vec_hopf),	 # used to set ϕ, see the phase constraint
	hopfpt.u,          # used to set uhopf, see the phase constraint
	M, 2n;             # number of time slices
	jacobian = :FullSparseInplace) # jacobian of PO functional			      

nothing #hide 	
```

To evaluate the functional at `x`, you call it like a function: `poTrap(x, par)` for the parameter `par`.

!!! note "Using the functional for deflation, Fold of limit cycles..."
    The functional `poTrap` gives you access to the underlying methods to call a regular `newton`. For example the functional is `x -> poTrap(x, par)` at parameter `par`. The (sparse) Jacobian at `(x,p)` is computed like this `poTrap(Val(:JacFullSparse), x, p)` while the Matrix Free version is `dx -> poTrap((x, p, dx)`. This also allows you to call the newton deflated method (see [Deflated problems](@ref)) or to locate Fold point of limit cycles see [`PeriodicOrbitTrapProblem`](@ref). You can also use preconditioners. In the case of more computationally intense problems (like the 2d Brusselator), this might be mandatory as using LU decomposition for the linear solve will use too much memory. See also the example [2d Ginzburg-Landau equation (finite differences, codim 2, Hopf aBS)](@ref cgl)

For convenience, we provide a simplified newton / continuation methods for periodic orbits. One has just to pass a [`PeriodicOrbitTrapProblem`](@ref).

```@example TUTBRUmanual
# we use the linear solver LSFromBLS to speed up the computations
opt_po = NewtonPar(tol = 1e-10, verbose = true, max_iterations = 14, linsolver = BK.LSFromBLS())

# we set the parameter values
poTrap = @set poTrap.prob_vf.params = (@set par_bru.l = l_hopf + 0.01)

outpo_f = newton(poTrap, orbitguess_f, (@set opt_po.verbose = false), normN = norminf) # hide
outpo_f = @time newton(poTrap, orbitguess_f, opt_po, normN = norminf)
BK.converged(outpo_f) && printstyled(color=:red, "--> T = ", outpo_f.u[end], "\n")
# plot of the periodic orbit
BK.plot_periodic_potrap(outpo_f.u, n, M; ratio = 2)
```

Finally, we can perform continuation of this periodic orbit using the specialized call `continuationPOTrap`

```@example TUTBRUmanual
opt_po = @set opt_po.eigsolver = EigArpack(; tol = 1e-5, v0 = rand(2n))
opts_po_cont = ContinuationPar(dsmin = 0.001, dsmax = 0.03, ds= 0.01,
	p_max = 3.0, max_steps = 20,
	newton_options = opt_po, nev = 5, tol_stability = 1e-8, detect_bifurcation = 0)

br_po = continuation(poTrap,
	outpo_f.u, PALC(),
	opts_po_cont;
	verbosity = 2,	plot = true,
	plot_solution = (x, p;kwargs...) -> heatmap!(reshape(x[1:end-1], 2*n, M)'; ylabel="time", color=:viridis, kwargs...),
	normC = norminf)

Scene = title!("")
```

## Deflation for periodic orbit problems
Looking for periodic orbits branching of bifurcation points, it is very useful to use `newton` algorithm with deflation. We thus define a deflation operator (see previous example)

```julia
deflationOp = DeflationOperator(2, (x,y) -> dot(x[1:end-1], y[1:end-1]), 1.0, [zero(orbitguess_f)])
```

which allows to find periodic orbits different from `orbitguess_f `. Note that the `dot` product removes the last component, *i.e.* the period of the cycle is not considered during this particular deflation. We can now use

```julia
outpo_f = @time newton(poTrap, orbitguess_f, deflationOp, opt_po; normN = norminf)
```

## Floquet coefficients

A basic method for computing Floquet coefficients based on the eigenvalues of the monodromy operator is available (see [`FloquetQaD`](@ref)). It is precise enough to locate bifurcations. Their computation is triggered like in the case of a regular call to `continuation`:

```julia
opt_po = @set opt_po.eigsolver = DefaultEig()
opts_po_cont = ContinuationPar(dsmin = 0.001, dsmax = 0.04, ds= -0.01, p_max = 3.0, max_steps = 200, newton_options = opt_po, nev = 5, tol_stability = 1e-6)
br_po = @time continuation(poTrap, outpo_f.u, PALC(),
	opts_po_cont; verbosity = 3, plot = true,
	plot_solution = (x, p;kwargs...) -> heatmap!(reshape(x[1:end-1], 2*n, M)'; ylabel="time", color=:viridis, kwargs...), normC = norminf)
```

A more complete diagram can be obtained combining the methods (essentially deflation and Floquet) described above. It shows the period of the periodic orbits as function of `l`. See `example/brusselator.jl` for more information.

![](bru-po-cont-3br.png)

!!! danger "Floquet multipliers computation"
    The computation of Floquet multipliers is necessary for the detection of bifurcations of periodic orbits (which is done by analyzing the Floquet exponents obtained from the Floquet multipliers). Hence, the eigensolver needs to compute the eigenvalues with largest modulus (and not with largest real part which is their default behavior). This can be done by changing the option `which = :LM` of the eigensolver. Nevertheless, note that for most implemented eigensolvers in the current Package, the proper option is set when the computation of Floquet multipliers is requested.

!!! tip "Performances"
    This example is clearly not optimized because we wanted to keep it simple. We can use a Matrix-Free version of the functional and preconditioners to speed this up. Floquet multipliers could also be computed in a Matrix-Free manner. See `examples/brusselator.jl` for more efficient methods. See also [2d Ginzburg-Landau equation (finite differences, codim 2, Hopf aBS)](@ref cgl) for a more advanced example where we introduce those methods.

## Continuation of periodic orbits (Standard Shooting)

> Note that what follows is not really optimized on the `DifferentialEquations.jl` side. Indeed, we do not use automatic differentiation, we do not pass the sparsity pattern,...

We now turn to a different method based on the flow of the Brusselator. To compute this flow (time stepper), we need to be able to solve the differential equation (actually a PDE) associated to the vector field `Fbru`. We will show how to do this with an implicit method `Rodas4P` from `DifferentialEquations.jl`. Note that the user can pass its own time stepper but for convenience, we use the ones in `DifferentialEquations.jl`. More information regarding the shooting method is contained in [Periodic orbits based on the shooting method](@ref).

```julia
n = 100

# different parameters to define the Brusselator model and guess for the stationary solution
par_bru = (α = 2., β = 5.45, D1 = 0.008, D2 = 0.004, l = 0.3)
sol0 = vcat(par_bru.α * ones(n), par_bru.β/par_bru.α * ones(n))
probBif = re_make(probBif, u0 = sol0)

eigls = EigArpack(1.1, :LM)
opts_br_eq = ContinuationPar(dsmin = 0.001, dsmax = 0.00615, ds = 0.0061, p_max = 1.9,
	detect_bifurcation = 3, nev = 21, plot_every_step = 50,
	newton_options = NewtonPar(eigsolver = eigls, tol = 1e-9), max_steps = 1060)

br = @time continuation(probBif, PALC(), opts_br_eq, normC = norminf)
```

We need to create a guess for the periodic orbit. We proceed as previously:

```julia
# number of time slices
M = 5

# index of the Hopf point in the branch br
ind_hopf = 1

l_hopf, Th, orbitguess2, hopfpt, vec_hopf = BK.guess_from_hopf(br, ind_hopf,
	opts_br_eq.newton_options.eigsolver, M, 22*0.075)
#
orbitguess_f2 = reduce(hcat, orbitguess2)
orbitguess_f = vcat(vec(orbitguess_f2), Th) |> vec
```

Let us now initiate the Standard Shooting method. To this aim, we need to provide a guess of the periodic orbit at times $T/M_{sh}$ where $T$ is the period of the cycle and $M_{sh}$ is the number of slices along the periodic orbits. If $M_{sh} = 1$, this the Standard Simple Shooting and the Standard Multiple one otherwise. See [`ShootingProblem`](@ref) for more information.

```julia
dM = 2
orbitsection = Array(orbitguess_f2[:, 1:dM:M])

# the last component is an estimate of the period of the cycle.
initpo = vcat(vec(orbitsection), 3.0)
```

Finally, we need to build a problem which encodes the Shooting functional. This done as follows where we first create the time stepper. For performance reasons, we rely on `SparseDiffTools `

```julia
using DifferentialEquations, SparseDiffTools, SparseArrays

FOde(f, x, p, t) = Fbru!(f, x, p)

u0 = sol0 .+ 0.01 .* rand(2n)

# parameter close to the Hopf bifurcation point
par_hopf = (@set par_bru.l = l_hopf + 0.01)

jac_prototype = Jbru_sp(ones(2n), @set par_bru.β = 0)
jac_prototype.nzval .= ones(length(jac_prototype.nzval))
_colors = matrix_colors(jac_prototype)
vf = ODEFunction(FOde; jac_prototype = jac_prototype, colorvec = _colors)
prob = ODEProblem(vf,  sol0, (0.0, 520.), par_bru)
```

We create the parallel standard shooting problem:

```julia
# this encodes the functional for the Shooting problem
probSh = ShootingProblem(
	# we pass the ODEProblem encoding the flow and the time stepper
	prob, Rodas4P(),

	# this is for the phase condition, you can pass your own section as well
	[orbitguess_f2[:,ii] for ii=1:dM:M];

	# enable threading
	parallel = true,

	# these are options passed to the ODE time stepper
	abstol = 1e-10, reltol = 1e-8,

  # parameter axis
  lens = (@optic _.l),

  # parameters
  par = par_hopf,

  # jacobian of the periodic orbit functional
  jacobian = BK.FiniteDifferencesMF())
```

We are now ready to call `newton`

```julia
ls = GMRESIterativeSolvers(reltol = 1e-7, N = length(initpo), maxiter = 100)
optn_po = NewtonPar(verbose = true, tol = 1e-9,  max_iterations = 20, linsolver = ls)
outpo = @time newton(probSh, initpo, optn_po; normN = norminf)
plot(initpo[1:end-1], label = "Init guess")
plot!(outpo.u[1:end-1], label = "sol")
```

which gives (note that we did not have a really nice guess...)

```julia
┌─────────────────────────────────────────────────────┐
│ Newton step         residual      linear iterations │
├─────────────┬──────────────────────┬────────────────┤
│       0     │       2.6239e-01     │        0       │
│       1     │       8.8979e-03     │       27       │
│       2     │       9.6221e-04     │       28       │
│       3     │       5.0540e-02     │       31       │
│       4     │       3.4149e-03     │       28       │
│       5     │       1.1336e+00     │       31       │
│       6     │       9.5418e-02     │       25       │
│       7     │       1.1996e-02     │       28       │
│       8     │       1.2655e-02     │       28       │
│       9     │       3.8820e-03     │       30       │
│      10     │       4.1901e-04     │       31       │
│      11     │       5.4309e-06     │       31       │
│      12     │       1.1430e-09     │       33       │
│      13     │       5.2525e-14     │       34       │
└─────────────┴──────────────────────┴────────────────┘
  3.396345 seconds (5.38 M allocations: 17.675 GiB, 1.66% gc time)
```

and

![](brus-sh-new.png)

Note that using Simple Shooting, the convergence is much faster. Indeed, running the code above with `dM = 10` gives:

```julia
┌─────────────────────────────────────────────────────┐
│ Newton step         residual     linear iterations  │
├─────────────┬──────────────────────┬────────────────┤
│       0     │       6.1712e-03     │        0       │
│       1     │       3.4465e-03     │        6       │
│       2     │       1.0516e-01     │        8       │
│       3     │       7.4614e-03     │        6       │
│       4     │       1.6620e-03     │        7       │
│       5     │       3.9589e-04     │        7       │
│       6     │       4.3043e-05     │        8       │
│       7     │       1.7232e-06     │        8       │
│       8     │       8.0455e-09     │        8       │
│       9     │       3.9453e-11     │        8       │
└─────────────┴──────────────────────┴────────────────┘
  0.612070 seconds (217.21 k allocations: 523.069 MiB, 4.83% gc time)
```

!!! info "Convergence and speedup"
    The convergence is much worse for the multiple shooting than for the simple one. This is reflected above in the number of linear iterations made during the newton solve. The reason for this is because of the cyclic structure of the jacobian which impedes GMRES from converging fast. This can only be resolved with an improved GMRES which we'll provide in the future.


Finally, we can perform continuation of this periodic orbit using a specialized version of `continuation`:

```julia
# note the eigensolver computes the eigenvalues of the monodromy matrix. Hence
# the dimension of the state space for the eigensolver is 2n
opts_po_cont = ContinuationPar(dsmin = 0.001, dsmax = 0.05, ds= 0.01, p_max = 1.5,
	max_steps = 500, newton_options = (@set optn_po.tol = 1e-7), nev = 25,
	tol_stability = 1e-8, detect_bifurcation = 0)

br_po = @time continuation(probSh, outpo.u, PALC(),
	opts_po_cont; verbosity = 2,
	# specific bordered linear solver
	linear_algo = MatrixFreeBLS(@set ls.N = ls.N+1),
	plot = true,
	plot_solution = (x, p; kwargs...) -> BK.plot_periodic_shooting!(x[1:end-1], length(1:dM:M); kwargs...),
	record_from_solution = (u, p; k...) -> u[end], 
	normC = norminf)
```

We can observe that simple shooting is faster but the Floquet multipliers are less accurate than for multiple shooting. Also, when the solution is very unstable, simple shooting can have spurious branch switching. Finally, note the $0=\log 1$ eigenvalue of the monodromy matrix in the graph below.

![](brus-sh-cont.png)

## Continuation of periodic orbits (Poincaré Shooting)

We now turn to another Shooting method, namely the Poincaré one. We can provide this method thanks to the unique functionalities of `DifferentialEquations.jl`. More information is provided at [`PoincareShootingProblem`](@ref) and [Periodic orbits based on the shooting method](@ref) but basically, it is a shooting method between Poincaré sections $\Sigma_i$ (along the orbit) defined by hyperplanes. As a consequence, the dimension of the unknowns is $M_{sh}\cdot(N-1)$ where $N$ is the dimension of the phase space. Indeed, each time slice lives in an hyperplane $\Sigma_i$. Additionally, the period $T$ is not an unknown of the method but rather a by-product. However, the method requires the time stepper to find when the flow hits an hyperplane $\Sigma_i$, something called **event detection**.


We show how to use this method, the code is very similar to the case of the Standard Shooting. We first define the functional for Poincaré Shooting Problem

```julia
# sub-sampling factor of a initial guess for the periodic orbit
dM = 5

# vectors to define the hyperplanes Sigma_i
normals = [Fbru(orbitguess_f2[:,ii], par_hopf)/(norm(Fbru(orbitguess_f2[:,ii], par_hopf))) for ii = 1:dM:M]
centers = [orbitguess_f2[:,ii] for ii = 1:dM:M]

# functional to hold the Poincare Shooting Problem
probHPsh = PoincareShootingProblem(
	# ODEProblem, ODE solver used to compute the flow
	prob, Rodas4P(),

	# parameters for the Poincaré sections
	normals, centers;

	# enable threading
	parallel = true,

	# Parameters passed to the ODE solver
	abstol = 1e-10, reltol = 1e-8,

  # parameter axis
  lens = (@optic _.l),

  # parameters
  par = par_hopf,

  # jacobian of the periodic orbit functional
  jacobian = BK.FiniteDifferencesMF())
```

Let us now compute an initial guess for the periodic orbit, it must live in the hyperplanes $\Sigma_i$. Fortunately, we provide projections on these hyperplanes.

```julia
# projection of the initial guess on the hyperplanes. We assume that the centers[ii]
# form the periodic orbit initial guess.
initpo_bar = reduce(vcat, BK.projection(probHPsh, centers))
```

We can now call `continuation` to get the first branch.

```julia
# eigen / linear solver
eig = EigKrylovKit(tol= 1e-12, x₀ = rand(2n-1), dim = 40)
ls = GMRESIterativeSolvers(reltol = 1e-11, N = length(vec(initpo_bar)), maxiter = 500)

# newton options
optn = NewtonPar(verbose = true, tol = 1e-9,  max_iterations = 140, linsolver = ls)

# continuation options
opts_po_cont_floquet = ContinuationPar(dsmin = 0.0001, dsmax = 0.05, ds= 0.001,
	p_max = 2.5, max_steps = 500, nev = 10,
	tol_stability = 1e-5, detect_bifurcation = 3, plot_every_step = 3)
opts_po_cont_floquet = @set opts_po_cont_floquet.newton_options =
	NewtonPar(linsolver = ls, eigsolver = eig, tol = 1e-9, verbose = true)

# continuation run
br_po = @time continuation(probHPsh, vec(initpo_bar), PALC(),
	opts_po_cont_floquet; verbosity = 3,
	linear_algo = MatrixFreeBLS(@set ls.N = ls.N+1),
	plot = true,
	plot_solution = (x, p; kwargs...) -> BK.plot!(x; label="", kwargs...),
	normC = norminf)		
```

![](brus-psh-cont.png)

We also obtain the following information:

```julia
julia> br_po
Branch number of points: 41
Bifurcation points:
- #  1,    bp at p ≈ 1.20987963 ∈ (1.20128196, 1.20987963), |δp|=9e-03, [converged], δ = ( 1,  0), step =  21, eigenelements in eig[ 22], ind_ev =   1
- #  2,    ns at p ≈ 1.78687615 ∈ (1.77831727, 1.78687615), |δp|=9e-03, [converged], δ = ( 2,  2), step =  30, eigenelements in eig[ 31], ind_ev =   3
- #  3,    pd at p ≈ 1.85103701 ∈ (1.84676466, 1.85103701), |δp|=4e-03, [converged], δ = ( 1,  1), step =  31, eigenelements in eig[ 32], ind_ev =   4
- #  4,    ns at p ≈ 1.87667870 ∈ (1.86813520, 1.87667870), |δp|=9e-03, [converged], δ = ( 2,  2), step =  32, eigenelements in eig[ 33], ind_ev =   6
```

## References

[^Lust]:> **Numerical Bifurcation Analysis of Periodic Solutions of Partial Differential Equations,** Lust, 1997.


