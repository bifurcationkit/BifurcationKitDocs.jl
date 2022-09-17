# we use this hacky way because AsymptoticNumericalMethod is not registered
using Pkg
pkg"add https://github.com/bifurcationkit/AsymptoticNumericalMethod.jl"

using Documenter, BifurcationKit, Setfield, AsymptoticNumericalMethod
# using DocThemeIndigo
ENV["GKSwstype"] = "100"

makedocs(doctest = false,
	sitename = "Bifurcation Analysis in Julia",
	format = Documenter.HTML(collapselevel = 1,assets = ["assets/indigo.css"]),
	# format = DocumenterLaTeX.LaTeX(),
	authors = "Romain Veltz",
	pages = Any[
		"Home" => "index.md",
		"Educational introduction" => "educational.md",
		"Tutorials" => "tutorials/tutorials.md",
		"Basics" => [
			"Overview" => "guidelines.md",
			"Plotting" => "plotting.md",
			],
		"Problems" => [
			"Bifurcation Problem" => "BifProblem.md",
			"DiffEq wrapper" => "diffeq.md",
			"Periodic Orbits" => [
				"Introduction" => "periodicOrbit.md",
				"Trapezoid" => "periodicOrbitTrapeze.md",
				"Collocation" => "periodicOrbitCollocation.md",
				"Shooting" => "periodicOrbitShooting.md",
				],
			"Symmetries / Waves" => [
				"Introduction" => "intro_wave.md",
				"Eigen Solvers" => "waveEigen.md",
				"Modulated Travelling waves" => "ModulatedTW.md",
				],
		],
		"Continuation methods" => [
			"Introduction" => "IntroContinuation.md",
			"Predictors / correctors" => "Predictors.md",
			"PALC" => "PALC.md",
			"Moore-Penrose Continuation" => "MooreSpence.md",
			"ANM" => "ANM.md",
			"Deflated Continuation" => "DeflatedContinuation.md",
			"Event Handling and Callback" => "EventCallback.md",
				],
		"Functionalities" => [
			"Nonlinear Equations (Newton)" => ["newton.md", "deflatedproblem.md"],
			"Bifurcations" => [
				"Bifurcation detection (codim 1)" => "detectionBifurcation.md",
				"Fold / Hopf Continuation (codim 2)" => "codim2Continuation.md",
				],
			"Normal form" => [
				"Simple branch point" => "simplebp.md",
				"Non-simple branch point" => "nonsimplebp.md",
				"Simple Hopf point" => "simplehopf.md",
				"Cusp point" => "cusp.md",
				"Bogdanov-Takens point" => "bt.md",
				"Bautin point" => "bautin.md",
			],
			"Branch switching" => "branchswitching.md",
			"Bifurcation diagram" => "BifurcationDiagram.md",
			"Constrained problem" => "constrainedproblem.md",
			"Iterator Interface" => "iterator.md",
		],
		"Options" => [
			"Linear Solvers" => "linearsolver.md",
			"Bordered linear solvers" => "borderedlinearsolver.md",
			"Eigen Solvers" => "eigensolver.md",
			"Bordered arrays" => "Borderedarrays.md",
		],
		"Contributing" => [
			"Interfaces" => [
				"Vector" => "interfaceLS.md",
				"Linear / Eigen Solvers" => "interfaceLS.md",
				"Predictor / Corrector" => "interfacePred.md",
				"Flow" => "interfaceFlow.md",
				]
		],
		"Frequently Asked Questions" => "faq.md",
		"Library" => "library.md"
	]
	)

deploydocs(
	repo = "github.com/bifurcationkit/BifurcationKitDocs.jl.git",
	devbranch = "main"
)
