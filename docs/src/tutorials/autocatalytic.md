# 🟤 Modulated fronts in 1d autocatalytic model (Manual)

```@contents
Pages = ["autocatalytic.md"]
Depth = 3
```

We consider the model [^Balmforth][^Malham] which is also treated in [^Beyn]

$$\begin{array}{l}
u_{t}=a u_{x x}-u f(v), \quad a>0, u, v: \mathbb{R} \rightarrow \mathbb{R} \\
v_{t}=v_{x x}+u f(v)
\end{array}$$

where $f(u) = u^m 1_{u\geq 0}$. We chose the boundary conditions

$$\left(u_{-}, v_{-}\right)=(0,1),\quad \left(u_{+}, v_{+}\right)=(1,0)\tag{BC}.$$

It is straightforward to implement this problem as follows:

```@example TUTAUTOCAT
using Revise
using ForwardDiff, SparseArrays
using BifurcationKit, LinearAlgebra, Plots
const BK = BifurcationKit

# supremum norm
f(u) = u^9 # solutions are positive, so remove the heaviside

# helper function to plot solutions
function plotsol!(x; k...)
	u = @view x[1:end÷2]
	v = @view x[end÷2:end]
	plot!(u; label="u", k...)
	plot!(v; label="v", k...)
end
plotsol(x; k...) = (plot();plotsol!(x; k...))

# encode the nonlinearity
@views function NL!(dest, U, p, t = 0.)
	N = p.N
	u = U[1:N]
	v = U[N+1:2N]
	dest[1:N]    .= -u .* f.(v)
	dest[N+1:2N] .=  -dest[1:N]#u .* f.(v)
	return dest
end

# function for the differential with specific boundary conditions
# for fronts
@views function applyD_add!(f, U, p, a)
	uL = 0; uR = 1; vL = 1; vR = 0
	n = p.N
	u = U[1:n]
	v = U[n+1:2n]

	c1 = 1 / (2p.h)
	f[1]   += a * (uL      - u[2] ) * c1
	f[end] += a * (v[n-1]  - vR   ) * c1

	f[n]   += a * (u[n-1] - uR  ) * c1
	f[n+1] += a * (    vL - v[2] ) * c1

	@inbounds for i=2:n-1
		  f[i] += a * (u[i-1] - u[i+1]) * c1
		f[n+i] += a * (v[i-1] - v[i+1]) * c1
	end
	return f
end

# function which encodes the full PDE
@views function Fcat!(f, U, p, t = 0)
	uL = 0; uR = 1; vL = 1; vR = 0
	n = p.N
	# nonlinearity
	NL!(f, U, p)

	# Dirichlet boundary conditions
	h2 = p.h * p.h
	c1 = 1 / h2

	u = U[1:n]
	v = U[n+1:2n]

	f[1]   += p.a * (uL      - 2u[1] + u[2] ) * c1
	f[end] +=       (v[n-1]  - 2v[n] + vR   ) * c1

	f[n]   += p.a * (u[n-1] - 2u[n] +  uR  ) * c1
	f[n+1] +=       (vL - 2v[1] + v[2] ) * c1

	@inbounds for i=2:n-1
		  f[i] += p.a * (u[i-1] - 2u[i] + u[i+1]) * c1
		f[n+i] +=       (v[i-1] - 2v[i] + v[i+1]) * c1
	end
	return f
end
Jcat(x,p) = sparse(ForwardDiff.jacobian(x -> Fcat!(similar(x), x, p), x))
nothing #hide
```

We chose the following parameters:

```@example TUTAUTOCAT
N = 200
lx = 25.
X = LinRange(-lx,lx, N)
par_cat = (N = N, a = 0.18, h = 2lx/N)

u0 = @. (tanh(2X)+1)/2
U0 = vcat(u0, 1 .- u0)

# we define a problem to hold the vector field
prob = BifurcationProblem(Fcat!, u0, par_cat, (@optic _.a); J = Jcat)
nothing #hide
```

## Freezing method

The problem may feature fronts, solutions of the form $u(x,t) = \tilde u(x-st)$ (same for $v$) for a fixed value of the profile $\tilde u$ and the speed $s$. The equation for the front profile is, up to an abuse of notations (we removed the tildes)

$$\begin{array}{l}
0=a u_{\xi\xi}+s\cdot u_{\xi}-u f(v)\\
0=v_{\xi\xi}+s\cdot v_{\xi}+u f(v)
\end{array}$$

with unknowns $u,v,s$. The front is solution of these equations but it is not uniquely determined because of the phase invariance. Hence, we add the phase condition (see [^Beyn])

$$0 = \left\langle (u,v), \partial_\xi (u_0,v_0) \right\rangle$$

where $U_0:=(u_0,v_0)$ is some fixed profile. This is easily coded in the following functional

```@example TUTAUTOCAT
@views function FcatWave!(out, x, p)
	N = p.N
	U = x[1:end-1]
	Fcat!(out[1:2N], U, p)
	applyD_add!(out[1:2N], U, p, x[end])
	# phase condition
	out[2N+1] = dot(U, p.Du0)
	return out
end
FcatWave(x, p, t = 0) = FcatWave!(similar(x), x, p)
JcatWave(u, p) = sparse(ForwardDiff.jacobian(z -> FcatWave!(similar(z),z,p), u))
nothing #hide
```

We now define the $U_0$ profile

```@example TUTAUTOCAT
uold = vcat(u0, (1 .- u0))
Duold = zero(uold); applyD_add!(Duold, uold, par_cat,1)

# update problem parameters for front problem
par_cat_wave = (par_cat..., u0Du0 = dot(uold, Duold), Du0 = Duold, uold = uold)
nothing #hide
```

Let us find the front using `newton`

```@example TUTAUTOCAT
# we define a problem for solving for the wave
probtw = BifurcationProblem(FcatWave!, vcat(U0, -1.), par_cat_wave, (@optic _.a);
	J = JcatWave,
	record_from_solution = (x,p;k...) -> (s = x[end], nrm = norm(x[1:end-1])),
	plot_solution = (x, p; k...) -> plotsol!(x[1:end-1];k...))

front = BK.solve(probtw, Newton(), NewtonPar())
println("front speed s = ", front.u[end], ", norm = ", front.u[1:end-1] |> norminf)
```

```@example TUTAUTOCAT
plotsol(front.u[1:end-1], title="front solution")
```

## Continuation of front solutions

Following [^Malham], the modulated fronts are solutions of the following DAE

$$\begin{array}{l}\tag{DAE}
u_{t}=a u_{x x}+s\cdot u_x-u f(v)\\
v_{t}=v_{x x}+s\cdot v_x+u f(v)\\
0 = \left\langle U, \partial_\xi U_0	\right\rangle
\end{array}$$

which can be written with a PDE $M_aU_t = G(u)$ with mass matrix $M_a = (Id, Id, 0)$. We have already written the vector field of (MF) in the function `FcatWave`.

Having found a front $U^f$, we can continue it as function of the parameter $a$ and detect instabilities. The stability of the front is linked to the eigenelements $(\lambda, V)$ solution of the generalized eigenvalue problem:

$$\lambda M_a\cdot V = dG(U^f)\cdot V.$$


However `BifurcationKit` does not provide a generalized eigenvalue solver for now, so we devise one:

```@example TUTAUTOCAT
# we need  a specific eigensolver
struct EigenWave <: BK.AbstractEigenSolver end

# implementation of the solver for the generalized Eigen problem
function (eig::EigenWave)(Jac, nev; k...)
	N = size(Jac,1)
	B = diagm(vcat(ones(N-1),0))
	F = eigen(Array(Jac), B)
	I = sortperm(F.values, by = real, rev = true)
	nev2 = min(nev, length(I))
	J = findall( abs.(F.values[I]) .< 100000)
	return Complex.(F.values[I[J[1:nev2]]]), Complex.(F.vectors[:, I[J[1:nev2]]]), true, 1
end

optn = NewtonPar(tol = 1e-8, eigsolver = EigenWave())
opt_cont_br = ContinuationPar(p_min = 0.05, p_max = 1., newton_options = optn, ds= -0.001, plot_every_step = 2, detect_bifurcation = 3, nev = 10, n_inversion = 6)
br = continuation(probtw, PALC(), opt_cont_br)
plot(br)
```

We have detected a Hopf instability in front dynamics, this will give rise of modulated fronts. Let us try to compute them.

## Branch of modulated fronts

To branch from the Hopf bifurcation point, we just have to pass the mass matrix as follows:

```@example TUTAUTOCAT
# we compute the periodic solutions using Mt time steps and a Trapezoidal time stepper
# note that we pass the parameter massmatrix which
# allows to solver the DAE
Mt = 30
probTP = PeriodicOrbitTrapProblem(M = Mt ;
		massmatrix = spdiagm(0 => vcat(ones(2N),0.)),
		update_section_every_step = 1,
		# linear solver for the periodic orbit problem
		# OPTIONAL, one could use the default
		jacobian = BK.BorderedLU())

opts_po_cont = ContinuationPar(dsmin = 0.0001, dsmax = 0.01, ds= -0.001, p_min = 0.05, max_steps = 130, newton_options = optn, nev = 7, tol_stability = 1e-3, detect_bifurcation = 0, plot_every_step = 1)
opts_po_cont = @set opts_po_cont.newton_options.max_iterations = 10
opts_po_cont = @set opts_po_cont.newton_options.tol = 1e-6

br_po = continuation(
	# we want to compute the bifurcated branch from
	# the first Hopf point
	br, 1,
	# arguments for continuation
	opts_po_cont,
	# this is how we pass the method to compute the periodic
	probTP ;
	# OPTIONAL parameters
	# we want to jump on the new branch at phopf + δp
	δp = -0.0025,
	# tangent predictor
	alg = PALC(tangent = Secant(),
			# linear solver specific to PALC
			bls = BorderingBLS(solver = DefaultLS(), check_precision = false)),
	# regular parameters for the continuation
	# a few parameters saved during run
	record_from_solution = (u, p; k...) -> begin
		outt = BK.get_periodic_orbit(p.prob, u, (@set  par_cat_wave.a=p))
		m = maximum(outt.u[end,:])
		return (s = m, period = u[end])
	end,
	# plotting of a section
	plot_solution = (x, p; k...) -> begin
		outt = BK.get_periodic_orbit(p.prob, x, (@set  par_cat_wave.a=p.p))
		plot!(outt.t, outt.u[end, :]; label = "", subplot=3)
		plot!(br, subplot=1)
	end,
	# print the Floquet exponent
	finalise_solution = (z, tau, step, contResult; k...) -> begin
		true
	end,
	plot = true,
	normC = norminf)

plot(br);plot!(br_po, label = "modulated fronts")
```

Let us plot one modulated front:

```@example TUTAUTOCAT
modfront = get_periodic_orbit(br_po, length(br_po))
plot(plot(modfront.t, modfront.u[end,:], xlabel = "t", ylabel = "s", label = ""),
	contour(modfront.t, X, modfront.u[1:N,:], color = :viridis, xlabel = "t", title = "u for a = $(round(br_po.sol[length(br_po)].p,digits=4))", fill = true, ylims=(-10,10)))
```


## References

[^Balmforth]:> N. J. Balmforth, R. V. Craster, and S. J. A. Malham. Unsteady fronts in an autocatalytic system. R. Soc. Lond. Proc. Ser. A Math. Phys. Eng. Sci., 455(1984):1401–1433, 1999.

[^Malham]:> S. J. A. Malham and M. Oliver. Accelerating fronts in autocatalysis. R. Soc. Lond. Proc. Ser. A Math. Phys. Eng. Sci., 456(1999):1609–1624, 2000.

[^Beyn]:> Beyn, Wolf-Jürgen, and Vera Thümmler. “Phase Conditions, Symmetries and PDE Continuation.” In Numerical Continuation Methods for Dynamical Systems: Path Following and Boundary Value Problems Springer Netherlands, 2007. https://doi.org/10.1007/978-1-4020-6356-5_10.
