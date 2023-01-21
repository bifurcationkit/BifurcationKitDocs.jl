# Migration from v0.1.x to v0.2.x

## IMPORTANT NOTICE
**New version of the package with modified interface. You are now required to define a `BifurcationProblem` to perform continuation or bifurcation analysis. The previous interface is available under the tag 0.1.12 which can be installed by doing**

<span style="color:red">`] add BifurcationKit@0.1.12`</span>

**The new version provides many bugs fix though.
(Please note that the docs are up to date).**

## Introduction

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