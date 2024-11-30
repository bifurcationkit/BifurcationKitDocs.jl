# Wrapper to the package `DifferentialEquations.jl`

!!! warning "Warning"
    This is work in progress.  

Several packages in the [SciML](https://sciml.ai) organization provide wrappers to `BifurcationKit`. On can mention

1. [ModelingToolkit](https://docs.sciml.ai/ModelingToolkit/stable/) and the [tutorials](https://docs.sciml.ai/ModelingToolkit/stable/tutorials/bifurcation_diagram_computation/)
2. [Catalyst](https://docs.sciml.ai/Catalyst/stable/) and the [tutorials](https://docs.sciml.ai/Catalyst/stable/steady_state_functionality/bifurcation_diagrams/)

# Work in progress

Use the [LinearSolve.jl](https://github.com/SciML/LinearSolve.jl) for handling the linear problems and also [NonlinearSolve.jl](https://github.com/SciML/NonlinearSolve.jl) whenever possible. The use of [LinearSolve.jl](https://github.com/SciML/LinearSolve.jl) allows to re-use Krylov spaces in between continuation steps.