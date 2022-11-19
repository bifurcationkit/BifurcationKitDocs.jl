# Migration from v0.1.x to v0.2.x

In v0.2.x, we introduced problem based bifurcation diagram, meaning that you have now to wrap your vector field in a `BifurcationProblem`. You also need to pass your plot/record functions.

## Don't use AD yourself

There is nothing wrong with doing so but this is done in the constructor of `BifurcationPoblem`, so 

```
prob = BifurcationProblem(F, x, p, lens ; J = myJac) 
```

should be 

```
prob = BifurcationProblem(F, x, p, lens) 
```

## Error: no method matching iterate(::BifurcationKit.ContResult 

This is because you use the old syntax 

```julia
br, = continuation(...)
```

instead of (no comma)

```julia
br = continuation(...)
```

## Arguments to `continuation`

`recordFromSolution` and `plotFromSolution` should be passed to `BifurcationProblem` instead of `continuation`.