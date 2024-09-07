# 🟡 Period doubling in the Barrio-Varea-Aragon-Maini model

```@contents
Pages = ["tutorialsPD.md"]
Depth = 3
```

The purpose of this example is to show how to handle period doubling bifurcations of periodic orbits. Note that we do not use automatic branch switching here although this is possible (but not yet for branching from period doubling points).

!!! info "Method and performance"
    We focus on the Shooting method but we could have based the computation of periodic orbits on finite differences instead. Performances of the current tutorial are directly linked to the ones of `DifferentialEquations.jl`.     

We focus on the following 1D model (see [^Aragon]):

$$\tag{E}\begin{aligned}
&\frac{\partial u}{\partial t}=D \nabla^{2} u+\eta\left(u+a v-C u v-u v^{2}\right)\\
&\frac{\partial v}{\partial t}=\nabla^{2} v+\eta\left(b v+H u+C u v+u v^{2}\right)
\end{aligned}$$
with Neumann boundary conditions. We start by encoding the model

```@example PDPDE
using Revise, ForwardDiff, DifferentialEquations, SparseArrays
using BifurcationKit, LinearAlgebra, Plots
const BK = BifurcationKit

f(u, v, p) = p.η * (      u + p.a * v - p.C * u * v - u * v^2)
g(u, v, p) = p.η * (p.H * u + p.b * v + p.C * u * v + u * v^2)

function NL!(dest, u, p, t = 0.)
	N = div(length(u), 2)
	u1 =  @view (u[1:N])
	u2 =  @view (u[N+1:end])
	dest[1:N]     .= f.(u1, u2, Ref(p))
	dest[N+1:end] .= g.(u1, u2, Ref(p))
	return dest
end

function Fbr!(f, u, p)
	NL!(f, u, p)
	mul!(f, p.Δ, u,1,1)
	f
end

NL(u, p) = NL!(similar(u), u, p)
Fbr(x, p, t = 0.) = Fbr!(similar(x), x, p)

# this is not very efficient but simple enough ;)
Jbr(x,p) = sparse(ForwardDiff.jacobian(x -> Fbr(x, p), x))
```

We can now perform bifurcation of the following Turing solution:

```@example PDPDE
N = 100
n = 2N
lx = 3pi/2
h = 2lx/N
X = LinRange(-lx,lx, N)

Δ = spdiagm(0 => -2ones(N), 1 => ones(N-1), -1 => ones(N-1) ) / h^2; Δ[1,1]=Δ[end,end]=-1/h^2
D = 0.08
par_br = (η = 1.0, a = -1., b = -3/2., H = 3.0, D = D, C = -0.6, Δ = blockdiag(D*Δ, Δ))

u0 = 1.0 * cos.(2X)
solc0 = vcat(u0, u0)

probBif = BK.BifurcationProblem(Fbr, solc0, par_br, (@optic _.C) ;J = Jbr,
		record_from_solution = (x, p; k...) -> norminf(x),
		plot_solution = (x, p; kwargs...) -> plot!(x[1:end÷2]; label="",ylabel ="u", kwargs...))

# parameters for continuation
eigls = EigArpack(0.5, :LM)
opt_newton = NewtonPar(eigsolver = eigls, tol=1e-9)
opts_br = ContinuationPar(dsmax = 0.04, ds = -0.01, p_min = -1.8,
	nev = 21, plot_every_step = 50, newton_options = opt_newton, max_steps = 400)

br = continuation(re_make(probBif, params = (@set par_br.C = -0.2)), PALC(), opts_br;
	plot = true, verbosity = 3)

plot(br)
```

## Periodic orbits from the Hopf point (Standard Shooting)

We continue the periodic orbit form the first Hopf point around $C\approx -0.8598$ using a Standard Simple Shooting method (see [Periodic orbits based on the shooting method](@ref)). To this end, we define a `SplitODEProblem` from `DifferentialEquations.jl` which is convenient for solving semilinear problems of the form

$$\dot x = Ax+g(x)$$

where $A$ is the infinitesimal generator of a $C_0$-semigroup. We use the exponential-RK scheme `ETDRK2` ODE solver to compute the solution of (E) just after the Hopf point.

```julia
# parameters close to the Hopf bifurcation
par_br_hopf = @set par_br.C = -0.86
# parameters for the ODEProblem
f1 = DiffEqArrayOperator(par_br.Δ)
f2 = NL!
prob_sp = SplitODEProblem(f1, f2, solc0, (0.0, 280.0), @set par_br.C = -0.86)

sol = @time solve(prob_sp, ETDRK2(krylov=true); abstol=1e-14, reltol=1e-14, dt = 0.1)
```

We estimate the period of the limit cycle to be around $T\approx 3$. We then use this as a guess for the shooting method:

```julia
# compute the guess for the shooting method
orbitsection = Array(sol[:, end])
initpo = vcat(vec(orbitsection), 3.)

# define the functional for the standard simple shooting based on the
# ODE solver ETDRK2. SectionShooting implements an appropriate phase condition
probSh = ShootingProblem(prob_sp, ETDRK2(krylov=true),
	[sol(280.0)]; abstol=1e-14, reltol=1e-12, dt = 0.1,
  optic = (@optic _.C),
  jacobian = BK.FiniteDifferencesMF())

# parameters for the Newton-Krylov solver
ls = GMRESIterativeSolvers(reltol = 1e-7, N = length(initpo), maxiter = 50, verbose = false)
optn = NewtonPar(verbose = true, tol = 1e-9,  max_iterations = 120, linsolver = ls)
# Newton-Krylov solver
out_po_sh = @time newton(probSh, initpo, optn; normN = norminf)
BK.converged(out_po_sh) && printstyled(color=:red, "--> T = ", out_po_sh.u[end], ", amplitude = ", BK.getamplitude(probSh, out_po_sh.u, par_br_hopf; ratio = 2),"\n")
```

which gives

```julia
--> T = 2.94557883943451, amplitude = 0.05791350025709674
```

We can now continue this periodic orbit:

```julia
eig = DefaultEig()
opts_po_cont = ContinuationPar(dsmin = 0.0001, dsmax = 0.01, ds= 0.005, p_min = -1.8, max_steps = 170, newton_options = (@set optn.eigsolver = eig),
	nev = 10, tol_stability = 1e-2)
br_po_sh = @time continuation(probSh, out_po_sh.u, PALC(), opts_po_cont; verbosity = 3,
	plot = true,
	linear_algo = MatrixFreeBLS(@set ls.N = probSh.M*n+2),
  plot_solution = (x, p; kwargs...) -> BK.plot_periodic_shooting!(x[1:end-1], 1; kwargs...),
  record_from_solution = (u, p; k...) -> BK.getmaximum(probSh, u, (@set par_br_hopf.C = p.p); ratio = 2), normC = norminf)
```

We plot the result using `plot(br_po_sh, br, label = "")`:

![](br_pd2.png)

!!! tip "Numerical precision for stability"
    The Floquet multipliers are not very precisely computed here using the Shooting method. We know that `1=exp(0)` should be a Floquet multiplier but this is only true here at precision ~1e-3. In order to prevent spurious bifurcation detection, there is a threshold `tol_stability` in `ContinuationPar` for declaring an unstable eigenvalue. Another way would be to use Poincaré Shooting so that this issue does not show up.

## Periodic orbits from the PD point (Standard Shooting)

We now compute the periodic orbits branching of the first Period-Doubling bifurcation point. It is straightforward to obtain an initial guess using the flow around the bifurcation point:

```julia
par_br_pd = @set par_br.C = -1.32
prob_sp = SplitODEProblem(f1, f2, solc0, (0.0, 300.0), par_br_pd)
# solution close to the PD point.
solpd = @time solve(prob_sp, ETDRK2(krylov=true); abstol=1e-14, reltol=1e-12, dt = 0.1)
```
The estimated period is $T_{pd}=6.2$:

```julia
orbitsectionpd = Array(solpd[:,end-100])
initpo_pd = vcat(vec(orbitsectionpd), 6.2)
```

For educational purposes, we show the newton outputs:

```julia
out_po_sh_pd = newton(BK.set_params_po(probSh, par_br_pd), initpo_pd , optn; normN = norminf)
BK.converged(out_po_sh_pd) && printstyled(color=:red, "--> T = ", out_po_sh_pd.u[end], ", amplitude = ", BK.getamplitude(probSh, out_po_sh_pd.u, (@set par_br.C = -0.86); ratio = 2),"\n")
```
which gives

```julia
┌─────────────────────────────────────────────────────┐
│ Newton step         residual     linear iterations  │
├─────────────┬──────────────────────┬────────────────┤
│       0     │       4.3061e+00     │        0       │
│       1     │       5.3254e-02     │       10       │
│       2     │       5.5053e-03     │       12       │
│       3     │       1.6203e-03     │       13       │
│       4     │       5.6236e-05     │       12       │
│       5     │       7.6290e-07     │       12       │
│       6     │       3.0634e-08     │       13       │
│       7     │       1.2282e-09     │       12       │
│       8     │       4.9233e-11     │       13       │
└─────────────┴──────────────────────┴────────────────┘
  4.249150 seconds (497.85 k allocations: 2.030 GiB, 11.13% gc time)
--> T = 6.126399996979465, amplitude = 1.410164896740365
```

We also compute the branch of periodic orbits using the following command:

```julia
opts_po_cont = ContinuationPar(dsmin = 0.0001, dsmax = 0.005, ds= 0.001, p_min = -1.8, max_steps = 100, newton_options = (@set optn.eigsolver = eig), nev = 5, tol_stability = 1e-3)
br_po_sh_pd = @time continuation(BK.set_params_po(probSh,par_br_pd), out_po_sh_pd.u, PALC(),	opts_po_cont;
  verbosity = 2, plot = true,
  linear_algo = MatrixFreeBLS(@set ls.N = probSh.M*n+2),
  plot_solution = (x, p; kwargs...) -> (BK.plot_periodic_shooting!(x[1:end-1], 1; kwargs...); plot!(br_po_sh; subplot=1, legend=false)),
  record_from_solution = (u, p; k...) -> BK.getmaximum(probSh, u, (@set par_br_pd.C = p.p); ratio = 2), normC = norminf)
```

and plot it using `plot(br_po_sh, br, br_po_sh_pd, label = "")`:

![](br_pd3.png)

# References

[^Aragon]:> Aragón, J. L., R. A. Barrio, T. E. Woolley, R. E. Baker, and P. K. Maini. “Nonlinear Effects on Turing Patterns: Time Oscillations and Chaos.” Physical Review E 86, no. 2 (2012)

