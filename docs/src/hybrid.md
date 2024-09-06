# Automatic switching between Natural and PALC 

This is one of the continuation methods implemented in the package. It is set by the option `AutoSwitch()` in [`continuation`](@ref). See also [`AutoSwitch`](@ref) for more information.

This algorithm uses a `Natural` continuation when the branch is "horizontal" and switches to `PALC` otherwise.

## Example

We provide an example of use. We define a `BifurcationProblem` as usual and pass the continuation algorithm `AutoSwitch`.

```@example AutoSwitch
using Plots
using BifurcationKit
const BK = BifurcationKit

function F(x, p)
	(;α) = p
	f = similar(x)

	f[1] = (-2x[1]+x[2]) + α * exp(x[1])
	f[2] = ( x[1]-2x[2]) + α * exp(x[2])

	return f
end

sol0 = zeros(2)
par = (α = 0.0, )
prob = BifurcationProblem(F, sol0, par, (@optic _.α); record_from_solution = (x,p) -> norminf(x))
```

```@example AutoSwitch
opts = ContinuationPar(dsmax = 0.15, max_steps = 40, )

br = continuation(prob, AutoSwitch(), opts, normC = norminf)
```

You can plot the result as usual:

```@example AutoSwitch
plot(br)
```