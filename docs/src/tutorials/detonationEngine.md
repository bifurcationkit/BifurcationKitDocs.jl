# [🟠 Detonation engine](@id detonation)

```@contents
Pages = ["detonationEngine.md"]
Depth = 3
```

This is a model of a detonation engine, a new kind of reactor developed for planes.
The model[^Koch] quantifies the spatio-temporal evolution of a property analogous to specific internal energy, $u(x,t)$, on a one dimensional (1D) periodic domain:


$$\left\{
\begin{gathered}
\frac{\partial u}{\partial t}=\nu_{1} \frac{\partial^{2} u}{\partial x^{2}}-u \frac{\partial u}{\partial x}+k q(1-\lambda) \exp \left(\frac{u-u_{c}}{\alpha}\right)-\epsilon u^{2} \\
\frac{\partial \lambda}{\partial t}=\nu_{2} \frac{\partial^{2} \lambda}{\partial x^{2}}+k(1-\lambda) \exp \left(\frac{u-u_{c}}{\alpha}\right) -\frac{s u_{p} \lambda}{1+\exp \left(r\left(u-u_{p}\right)\right)}
\end{gathered}\right.$$

## Problem discretization

We start by discretizing the above PDE based on finite differences.

```@example DETENGINE
using Revise, ForwardDiff
using DifferentialEquations, SparseArrays
using BifurcationKit, LinearAlgebra, Plots
const BK = BifurcationKit

ω(u, p) = p.k * exp((u - p.uc) / p.α)
β(u, p) = p.s * p.up /(1 + exp(p.r * (u - p.up)) )

# utilities for plotting solutions
function plotsol!(x; k...)
	n = length(x) ÷ 2
	u = @view x[1:n]
	v = @view x[n+1:2n]
	plot!(u; label="u", k...)
	plot!(v; label="λ", k...)
end
plotsol(x; k...) = (plot();plotsol!(x; k...))

# nonlinearity of the model
function NL!(dest, U, p, t = 0.)
	N = p.N
	u = @view U[1:N]
	λ = @view U[N+1:2N]
	dest[1:N]    .= p.q .* (1 .- λ) .* ω.(u, Ref(p)) .- p.ϵ .* u.^2
	dest[N+1:2N] .=        (1 .- λ) .* ω.(u, Ref(p)) .- λ .* β.(u, Ref(p))
	return dest
end

# function which encodes the right hand side of the PDE
@views function Fdet!(f, u, p, t = 0)
	N = p.N
	NL!(f, u, p) 						# nonlinearity
	mul!(f, p.Δ, u, p.ν1, 1) 			# put Laplacian
	f[1:p.N] .-= (p.D * u[1:N].^2) ./ 2	# add drift a la Burger
	return f
end
NL(U, p, t = 0.) = NL!(similar(U), U, p, t)
JNL(U, p, t = 0.) = ForwardDiff.jacobian(x -> NL(x, p), U)
Fdet(x, p, t = 0) = Fdet!(similar(x), x, p, t)
Jdet(x, p) = sparse(ForwardDiff.jacobian(x -> Fdet(x, p), x))
nothing #hide
```

We can now instantiate the model

```@example DETENGINE
N = 300
lx = 2pi
X = LinRange(0, lx, N)

h = lx/N
Δ = spdiagm(0 => -2ones(N), 1 => ones(N-1), -1 => ones(N-1) ) / h^2; Δ[1,end] = 1/h^2; Δ[end,1]=1/h^2
D = spdiagm(1 => ones(N-1), -1 => -ones(N-1)) / (2h); D[1,end] = -1/(2h);D[end, 1] = 1/(2h)

_ν1 = 0.0075
# model parameters
par_det = (N = N, q = 0.5, α = 0.3, up = 0.0, uc = 1.1, s = 3.5, k = 1., ϵ = 0.15, r = 5.0, ν1 = _ν1, ν2 = _ν1, Δ = blockdiag(Δ, Δ), D = D, Db = blockdiag(D, D))

# initial conditions
u0 = 0.5ones(N)
λ0 = 0.5ones(N)
U0 = vcat(u0, λ0)
U0cons = vcat(copy(U0), 1.)
nothing #hide
```

## Jacobian with sparsity detection

Writing the jacobian explicitly is cumbersome. We rely on automatic differentiation to get the sparse jacobian.

```@example DETENGINE
# improved jacobian with sparse coloring
using SparseArrays, SparseDiffTools, Test
const L1 = copy(Jdet(rand(2N), par_det))
const colors = matrix_colors(L1)
function JlgvfColorsAD(J, u, p, colors)
	SparseDiffTools.forwarddiff_color_jacobian!(J, (out, x) -> out .= Fdet(x,p), u, colorvec = colors)
	J
end
JdetAD(x,p) = JlgvfColorsAD(L1, x, p, colors)
nothing #hide
```

We are now ready to compute the bifurcation of the trivial (constant in space) solution:

```@example DETENGINE
# bifurcation problem
prob = BifurcationProblem(Fdet, U0, (par_det..., q = 0.5), (@optic _.up); 
	J = JdetAD,
	plot_solution = (x, p; k...) -> plotsol!(x; k...),
	record_from_solution = (x, p; k...) -> (u∞ = norminf(x[1:N]), n2 = norm(x)))

prob = re_make(prob, params = (@set par_det.up = 0.56))

# iterative eigen solver
eig = EigArpack(0.2, :LM, tol = 1e-13, v0 = rand(2N))
# eig = EigArnoldiMethod(sigma=0.2, which = BifurcationKit.LM(), x₀ = rand(2N ))

# newton options
optnew = NewtonPar(verbose = true, eigsolver = eig)
solhomo = BK.solve(prob, Newton(), optnew; normN = norminf)
optcont = ContinuationPar(newton_options = NewtonPar(optnew, verbose = false),
	detect_bifurcation = 3, nev = 50, n_inversion = 8, max_bisection_steps = 25,
	dsmax = 0.01, ds = 0.01, p_max = 1.4, max_steps = 1000, plot_every_step = 50)

br = continuation(
		re_make(prob, params = (@set par_det.q = 0.5), u0 = solhomo.u),
		PALC(), optcont; plot = true)
Scene = title!("")
```

```@example DETENGINE
br
```

We have detected 6 Hopf bifurcations. We now study the periodic orbits branching from them.

## Computing the branches of Traveling waves

The periodic orbits emanating from the Hopf points look like traveling waves. This is intuitive because the equation is mostly advective as the diffusion coefficients $\nu_i$ are small. We will thus seek for traveling waves instead of periodic orbits. The advantage is that the possible Neimark-Sacker bifurcation is transformed into a regular Hopf one which allows the study of modulated traveling waves.

As we will do the same thing 3 times, we bundle the procedure in functions. We first use the regular Hopf normal form to create a guess for the traveling wave:

```@example DETENGINE
function getGuess(br, nb; δp = 0.005)
	nf = get_normal_form(br, nb; verbose  = false)
	pred = predictor(nf, δp)
	return pred.p, pred.orbit(0)
end
nothing #hide
```

Using this guess, we can continue the traveling wave as function of a parameter. Note that in the following code, a generalized eigensolver is automatically created during the call to `continuation` which properly computes the stability of the wave.

```@example DETENGINE
function computeBranch(br, nb; δp = 0.005, max_steps = 190)
	_p, sol = getGuess(br, nb)
	# traveling wave problem
	probTW = TWProblem(
		re_make(br.prob, params = (getparams(br)..., up = _p)),
		getparams(br).Db,
		copy(sol),
		jacobian = :AutoDiff)
	# newton parameters with iterative eigen solver
	# eig = EigArnoldiMethod(sigma=0.2, which = BifurcationKit.LM(),x₀ = rand(2N ))
	eig = EigArpack(nev = 10, which = :LM, sigma = 0.4)
	optn = NewtonPar(verbose = true, eigsolver = eig)
	# continuation parameters
	opt_cont_br = ContinuationPar(p_min = 0.1, p_max = 1.3, newton_options = optn, ds= -0.001, dsmax = 0.01, plot_every_step = 5, detect_bifurcation = 3, nev = 10, max_steps = max_steps)
	# we build a guess for the traveling wave with speed -0.9
	twguess = vcat(sol, -0.9)
	br_wave = continuation(probTW, twguess, PALC(), opt_cont_br;
		verbosity = 3, plot = true, bothside = true,
		record_from_solution = (x, p; k...) -> (u∞ = maximum(x[1:N]), s = x[end], amp = amplitude(x[1:N])),
		plot_solution = (x, p; k...) -> (plotsol!(x[1:end-1];k...);plot!(br,subplot=1, legend=false)),
		callback_newton = BK.cbMaxNorm(1e2),
		finalise_solution = (z, tau, step, contResult; k...) -> begin
			amplitude(z.u[N+1:2N]) > 0.01
		end,
		)
end
nothing #hide
```

We can try this continuation as follows

```@example DETENGINE
amplitude(x) = maximum(x) - minimum(x)
br_wave = computeBranch(br, 1; max_steps = 10)
Scene = title!("")
```

## Building the full diagram

```@example DETENGINE
branches = [computeBranch(br, i) for i in 1:3]
plot(br, branches..., legend=:topleft, xlims = (0.5, 1.25), ylims=(0.5, 2.3))
```


## References

[^Koch]:> Koch, James, Mitsuru Kurosaka, Carl Knowlen, and J. Nathan Kutz. “Multi-Scale Physics of Rotating Detonation Engines: Autosolitons and Modulational Instabilities.” ArXiv:2003.06655 [Nlin, Physics:Physics], March 14, 2020. http://arxiv.org/abs/2003.06655.
