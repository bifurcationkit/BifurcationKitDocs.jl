# Interface for Predictor / Corrector

Here is a description of the interface that is used to specify predictor / corrector for continuation methods. The type must be a subtype of the abstract type `AbstractTangentPredictor`.

In the table below, we assume `it::AbstractContinuationIterable` and `M = BorderedArray `

| Required methods|   Brief description                                                                     |
|:------------------------------ |:------------------------------------------- |
| `getPredictor!(z_pred::M, z_old::M, τ::M, ds, pred::AbstractTangentPredictor, nrm = false)`                | Write in `z_pred` a prediction for the new Newton guess given the current solution `z_old`, tangent `τ` and arc length `ds`.        |
|`getTangent!(τ::M, z_new::M, z_old::M, it, ds, θ, pred:: AbstractTangentPredictor, verbosity)`| Generate an estimate of the tangent to the branch of solutions at positions `z_new`, `z_old`|
| `corrector(it, z_old::M, τ::M, z_pred::M, ds, θ, pred::AbstractTangentPredictor, linearalgo; normC = norm, callback = cbDefault, kwargs...)` | Correct the guess `z_pred`. Must return a tuple `(sol, residuals, isconverged, it numbers, itlinear)`|
| **Optional methods** | **Brief description**                                                                 |
| `Base.empty!(pred::AbstractTangentPredictor)`       |  Reset the predictor         |
