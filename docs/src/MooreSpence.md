# Moore-Penrose continuation

```@contents
Pages = ["MooreSpence.md"]
Depth = 3
```

This is one of the various continuation methods implemented in `BifurcationKit.jl`. It is set by the option `alg = MoorePenrose()` in [`continuation`](@ref). See also [`MoorePenrose`](@ref) for more information.

For solving

$$\mathbb R^n\ni F(x,p) = 0 \quad\tag{E}$$

using a Newton algorithm, we miss an equation. Hence, we proceed as follows [^Meijer]. Starting from a predictor $(x_1,p_1)$, we look for the solution to (E) that is closest to $(x_1,p_1)$. Hence, we optimise

$$\min_{(x,p)} \{ \|(x,p)-(x_1,p_1)\| \text{ such that } F(x,p)=0\} \tag{MS}$$  

It can be interpreted as a PALC in which the hyperplane is adapted at every step.  

## Predictor

The possible predictors `tangent::AbstractTangentPredictor` are listed in [Predictors - Correctors](@ref). They can be used to create a Moore-Penrose algorithm  like `MoorePenrose(tangent = PALC())`

## Corrector

The corrector is the Gauss Newton algorithm applied to (MS).

## Linear Algebra


### Norm

The option `normC` [`continuation`](@ref) specifies the norm used to evaluate the distance in (MS). The dot product (resp. norm) used in the (iterative) linear solvers is `LinearAlgebra.dot` (resp. `LinearAlgebra.norm`). It can be changed by importing these functions and redefining it. Note that by default, the ``L^2`` norm is used.

### Linear problem

The linear solver for the linear problem associated to (MS) is set by the option `linearAlgo` in [`continuation`](@ref): it is one of [Bordered linear solvers (BLS)](@ref).

## Algorithm for solving (MS)

Let us write $y\equiv(x,p)\in\mathbb R^{N+1}$.
In order to solve for the argmin, we apply the newton algorithm with jacobian belonging to $\mathbb R^{N\times (N+1)}$:

$$y^{k+1} = y^k -d_yF(y^k)^+F(y^k)$$

where the superscript $^+$ indicates the Moore-Penrose pseudoinverse of rank $N$.

### Direct case
In this case, triggered by the option `MoorePenrose(method = BifurcationKit.direct)`, the pseudoinverse is computed with `\`.

the option `MoorePenrose(method = BifurcationKit.pInv)`, the pseudoinverse is computed with `pinv`.

### Iterative case
In this case, triggered by the option `MoorePenrose(method = BifurcationKit.iterative)`, the pseudoinverse is computed with an iterative method described in [^Meijer]:

$$\left\{\begin{array}{l}
y_{1}^{j+1}=y_{1}^{j}-\left(\begin{array}{c}
F_{y}\left(y_{1}^{j}\right) \\
\left(\phi_{1}^{j}\right)^{\top}
\end{array}\right)^{-1}\left(\begin{array}{c}
F\left(y_{1}^{j}\right) \\
0
\end{array}\right) \\
\phi_{1}^{j+1}=\left(\begin{array}{c}
F_{y}\left(y_{1}^{j+1}\right) \\
\left(\phi_{1}^{j}\right)^{\top}
\end{array}\right)^{-1}\left(\begin{array}{l}
0 \\
1
\end{array}\right), \quad j=0,1,2, \ldots
\end{array}\right.$$

We initialise $\phi_1^0$ with the tangent.


## Step size control

Each time the corrector fails, the step size ``ds`` is halved. This has the disadvantage of having lost Newton iterations (which costs time) and imposing small steps (which can be slow as well). To prevent this, the step size is controlled internally with the idea of having a constant number of Newton iterations per point. This is in part controlled by the aggressiveness factor `a` in `ContinuationPar`.


## References

[^Meijer]:> Meijer, Dercole, and Oldeman, “Numerical Bifurcation Analysis.”
