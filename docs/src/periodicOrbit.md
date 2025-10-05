# Periodic orbits computation

```@contents
Pages = ["periodicOrbit.md"]
Depth = 3
```

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

It is important to understand the pros and cons of each method to compute PO in large dimensions.

In small dimensions, the collocation method is the preferred choice. It is the most precise, the fastest with dedicated linear solvers, time adaptive and feature full.


### Trapezoid method
The Trapezoid method (or the Collocation one) is usually faster than the ones based on Shooting but it requires more memory as it saves the whole orbit. However the main drawback of this method is that the associated linear solver is not "nice", being composed of a cyclic matrix for which no generic Matrix-free preconditioner is known. Hence, the Trapezoid method is **often used with an ILU preconditioner** which is severely constrained by memory. Also, when the period of the cycle is large, finer time discretization (or mesh adaptation which is not yet implemented) must be employed which is also a limiting factor both in term of memory and preconditioning.

### Collocation method

The Collocation method is (for now) the slowest of the three methods implemented for computing periodic orbits in large dimensions. However, it is by far the most precise one. Additionally, the mesh can be automatically adapted during the continuation. The implementation will be improved for large dimensional systems like the Trapezoid method one.

### Shooting method
The methods based on Shooting do not share the same drawbacks because the associated linear system is usually well conditioned, at least in the simple shooting case. They are thus often used **without preconditioner at all**. Even in the case of multiple shooting, this can be alleviated by a simple generic preconditioner based on deflation of eigenvalues (see [Linear solvers (LS)](@ref)). Also, the time stepper will automatically adapt to the stiffness of the problem, putting more time points where needed unlike the method based on finite differences which requires an adaptive (time) meshing to provide a similar property. Finally, we can use parallel Shooting to greatly decrease the computation time.

The main drawback of the method is to find a fast time stepper, at least to compete with the method based on finite differences. The other drawback is the precision of the method which cannot compete with the collocation one.

## Important notes

We regroup here some important notes which are valid for all methods above. 

### 1. Accessing the periodic orbit

In `record_from_solution`, `plot_solution` or after the computation of a branch of periodic orbits, how do I obtain the periodic orbit, for plotting purposes for example? If `x` is the solution from `newton` for the parameter `p`, you can obtain the periodic orbit as follows

```
xtt = get_periodic_orbit(x, p)
```

where `xtt.t` contains the time mesh and `xtt[:,:]` contains the different components. Note that for Trapezoid and collocation methods, calling `get_periodic_orbit` is essentially free as it is a reshape of `x`. However, in the case of Shooting methods, this requires recomputing the periodic orbit which can be costly for large scale problems.

### 2. Finaliser
If you pass a `finalise_solution` function to [`continuation`](@ref), the following occurs:

1. If the newton solve was successful, we update the phase condition every `update_section_every_step`
2. we call the user defined finalizer `finalise_solution`

### 3. Record from solution

You can pass your own function to [`continuation`](@ref). In the particular case of periodic orbits, the method is called like `record_from_solution(x, opt; k...)` where `opt.p` is the current value of the continuation parameter and `opt.prob` is the current state of the continuation problem. You can then obtain the current periodic orbit using (see above)

```julia
xtt = get_periodic_orbit(x, opt.p)
``` 

### 4. Plot solution

Similarly to `record_from_solution`, the method is called like `plot_solution(x, opt; k...)` where `opt.p` is the current value of the continuation parameter and `opt.prob` is the current state of the continuation problem.

### 5. Most precise method for Floquet coefficients
The state of the art method is based on a Periodic Schur decomposition. It is available through the package [PeriodicSchurBifurcationKit.jl](https://github.com/bifurcationkit/PeriodicSchurBifurcationKit.jl). For more information, have a look at `FloquetPQZ`.

### 6. Misc

`set_params_po`
`generate_ci_problem`
`generate_solution`
