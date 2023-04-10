# Periodic orbits computation

Consider the Cauchy problem

$$\frac{du}{dt}=F(u,p).$$

A periodic solution with period $T$ satisfies

$$\begin{align}
\frac{du}{dt}&=F(u,p)\\
u(0)&=u(T).
\end{align}$$

We provide 4 methods for computing periodic orbits (PO):

1. one (Trapezoid) based on finite differences to discretize a Cauchy problem,
2. one (Collocation) based on orthogonal collocation to discretize a Cauchy problem,
3. two (Shooting) based on the flow associated to a Cauchy problem.

It is important to understand the pro and cons of each method to compute PO in large dimensions.


### Trapezoid method
The Trapezoid method (or the Collocation one) is usually faster than the ones based or Shooting but it requires more memory as it saves the whole orbit. However the main drawback of this method is that the associated linear solver is not "nice", being composed of a cyclic matrix for which no generic Matrix-free preconditioner is known. Hence, the Trapezoid method is **often used with an ILU preconditioner** which is severely constrained by memory. Also, when the period of the cycle is large, finer time discretization (or mesh adaptation which is not yet implemented) must be employed which is also a limiting factor both in term of memory and preconditioning.

### Collocation method

The Collocation method is (for now) the slowest of the three methods implemented for computing periodic orbits. However, it is by far the most precise one. Additionally, the mesh can be automatically adapted during the continuation. The implementation will be improved for large dimensional systems like the Trapezoid method one.

### Shooting method
The methods based on Shooting do not share the same drawbacks because the associated linear system is usually well conditioned, at least in the simple shooting case. There are thus often used **without preconditioner at all**. Even in the case of multiple shooting, this can be alleviated by a simple generic preconditioner based on deflation of eigenvalues (see [Linear solvers (LS)](@ref)). Also, the time stepper will automatically adapt to the stiffness of the problem, putting more time points where needed unlike the method based on finite differences which requires an adaptive (time) meshing to provide a similar property. The main drawback of the method is to find a fast time stepper, at least to compete with the method based on finite differences. The other drawback is the precision of the method which cannot compete with the collocation method.

## Important notes

We regroup here some important notes which are valid for all methods above. 

### 1. Accessing the periodic orbit

In `recordFromSolution`, `plotoSlution` or after the computation of a branch of periodic orbits, how do I obtain the periodic orbit, for plotting purposes for example? If `x` it the solution from `newton` for example for the parameter `p`, you can obtain the periodic orbit as follows

```
xtt = getPeriodicOrbit(x, p)
```
where `xtt.t` contains the time mesh and `xtt[:,:]` contains the different components. Note that for Trapezoid and collocation methods, calling `getPeriodicOrbit` is essentially free as it is a reshape of `x`. However, in the case of Shooting methods, this requires recomputing the periodic orbit which can be costly for large scale problems.

### 2. Finaliser
If you pass a `finaliseSolution` function to [`continuation`](@ref), the following occurs:

1. If the newton solve was successful, we update the phase condition every `updateSectionEveryStep`
2. we call the user defined finalizer `finaliseSolution`

### 3. recordFromSolution

You can pass your own function to [`continuation`](@ref). In the particular case of periodic orbits, the method is called like `recordFromSolution(x, opt; k...)` where `opt.p` is the current value of the continuation parameter and `opt.prob` is the current state of the continuation problem. You can then obtain the current periodic orbit using (see above)

```julia
xtt = getPeriodicOrbit(x, opt.p)
``` 

### 4. plotSolution

Similarly to `recordFromSolution`, the method is called like `plotSolution(x, opt; k...)` where `opt.p` is the current value of the continuation parameter and `opt.prob` is the current state of the continuation problem.

### 5. Most precise method for Floquet coefficients
The state of the art method is based on a Periodic Schur decomposition. It is available through the package [PeriodicSchurBifurcationKit.jl](https://github.com/bifurcationkit/PeriodicSchurBifurcationKit.jl). For more information, have a look at `FloquetPQZ`.