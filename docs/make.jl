using Pkg
cd(@__DIR__)
pkg" activate ."
pkg" dev AbstractTrees BandedMatrices"
pkg" add BifurcationKit AsymptoticNumericalMethod"


using Documenter, BifurcationKit, AsymptoticNumericalMethod
# using DocumenterMermaid
# using DocThemeIndigo
ENV["GKSwstype"] = "100"

# to display progress
ENV["JULIA_DEBUG"] = Documenter

makedocs(
	modules = [BifurcationKit],
	doctest = false,
	pagesonly = false, # this is on Documenter#master, do not compile what is not in pages =
	draft = false,
	warnonly = true,
	sitename = "Bifurcation Analysis in Julia",
	format = Documenter.HTML(;
		collapselevel = 1,
		size_threshold_warn=300 * 2^10, # raise slightly from 100 to 200 KiB
		size_threshold=400 * 2^10,      # raise slightly 200 to to 300 KiB
	),# assets = ["assets/indigo.css"]),
	authors = "Romain Veltz",
	pages = Any[
		"Home" => "index.md",
		"Overview of capabilities" => "capabilities.md",
		"Getting Started with BifurcationKit" => "gettingstarted.md",
		"Tutorials" => "tutorials/tutorials.md",
		"Basics" => [
			"Guidelines" => "guidelines.md",
			"Educational introduction" => "educational.md",
			"Overview" => "overview.md",
			"Plot functions" => "plotting.md",
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
				# "Modulated Travelling waves" => "ModulatedTW.md",
				],
		],
		"Continuation methods" => [
			"Introduction" => "IntroContinuation.md",
			"Predictors / correctors" => "Predictors.md",
			"PALC" => "PALC.md",
			"Moore-Penrose continuation" => "MooreSpence.md",
			"AutoSwitch" => "hybrid.md",
			"ANM" => "ANM.md",
			"Deflated continuation" => "DeflatedContinuation.md",
				],
		"Functionalities" => [
			"Nonlinear equations (Newton)" => ["newton.md", "deflatedproblem.md"],
			"Bifurcations (equilibria)" => [
									"Bifurcation detection (1 param)" => "detectionBifurcation.md",
									"Fold / Hopf Continuation (2 params)" => "codim2Continuation.md",
									"Bogdanov-Takens refinement (3 params)" => "codim3Continuation.md",
								],
			"Bifurcations (periodic orbits)" => [
									"Bifurcation detection (1 param)" => "detectionBifurcationPO.md",
									"Fold continuation (2 params)" => "FoldContinuationPO.md",
									"Period-Doubling continuation (2 params)" => "PDContinuationPO.md",
									"Neimark-Sacker continuation (2 params)" => "NSContinuationPO.md",
								],
			"Normal form (equilibria)" =>	[
						"Simple branch point" => "simplebp.md",
						"Non-simple branch point" => "nonsimplebp.md",
						"Simple Hopf" => "simplehopf.md",
						"Cusp" => "cusp.md",
						"Bogdanov-Takens" => "bt.md",
						"Bautin" => "bautin.md",
						"Zero-Hopf" => "zh.md",
						"Hopf-Hopf" => "hh.md",
								],
			"Normal form (periodic orbits)" => [
						"Simple branch point" => "bppo.md",
						"Period-doubling" => "pd.md",
						"Neimark-Sacker" => "ns.md",
									],
			"Branch switching" => [
						"Introduction" => "intro-abs.md",
						"From equilibria to equilibria" => "abs-from-eq.md",
						"From Hopf/PD/Branch to periodic orbits" => "abs-from-hopf.md",
						"From codim 2 to equilibria" => "abs-from-codim2-eq.md",
						"From codim 2 to periodic orbits" => "abs-from-codim2-po.md",
						],
			"Automatic Bifurcation diagram" => "BifurcationDiagram.md",
			"Event handling and Callback" => "EventCallback.md",
			# "Constrained problem" => "constrainedproblem.md",
			"Iterator Interface" => "iterator.md",
		],
		"Options" => [
			"Linear solvers" => "linearsolver.md",
			"Bordered linear solvers" => "borderedlinearsolver.md",
			"Eigen solvers" => "eigensolver.md",
			"Bordered arrays" => "Borderedarrays.md",
		],
		"Contributing" => [
			"Interfaces" => [
				"Vector" => "interfaceLS.md",
				"Linear / eigen Solvers" => "interfaceLS.md",
				"Predictor / corrector" => "interfacePred.md",
				"Flow" => "interfaceFlow.md",
				]
		],
		"Frequently asked questions" => "faq.md",
		"Debugging" => "debug.md",
		"Migration from old versions" => "migration.md",
		"Library" => "library.md"
	]
	)

deploydocs(;
	repo = "github.com/bifurcationkit/BifurcationKitDocs.jl.git",
	push_preview = true, 
	target = "build", 
	devbranch="main"
	)
