# Results investigation

This section demonstrates how to diagnose and handle unexpected behavior in BifurcationKit. The primary strategy involves enabling verbose output by setting `verbose = true` in `NewtonPar(...)` and `verbosity = 3` in `continuation(...)` and making use of debug mode.

```julia
ENV["JULIA_DEBUG"] = BifurcationKit
```

Note that you can restore normal logging by using 

```julia
ENV["JULIA_DEBUG"] = Nothing
```