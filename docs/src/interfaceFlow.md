# Interface for Flow of Cauchy problem

Here is a description of the interface that is used to specify flows or semigroups of solutions of a Cauchy problem
$$\frac{du}{dt} = F(u, p),\quad u(0) = u_0.$$

More precisely, we call flow the mapping $Flow(x, p, t) = u(t)$.

The flow `fl` must be a subtype of the abstract type `AbstractFlow`. Note that in most cases, we only need `u(t)`. However, for plotting, we need **optionally** the full trajectory with time stamps in $[0,T]$.

Another point is the need for implementing serial versions of multithreaded flows.

| Required methods |   Brief description                                                                     |
|:------------------------------ |:------------------------------------------- |
| `vf(fl, x, p)`                | The vector field `F(x, p)` associated to a Cauchy problem. Used for the differential of the shooting problem. Must return `F(x, p)`       |
| `evolve(fl, x, par, t; k...)`        | the function implements the flow (or semigroup) `(x, p, t) -> flow(x, p, t)` associated to an autonomous Cauchy problem. Only the last time point must be returned in the form Named Tuple (u = ..., t = t). In the case of Poincaré Shooting, one must be able to call the flow like `evolve(fl, x, par, Inf)`. |
| **Optional methods** | **Brief description**                                                                 |
| `evolve(fl, x, par, dx, t; k...)`       |  The differential `dflow` of the flow *w.r.t.* `x`, `(x, p, dx, t) -> dflow(x, p, dx, t)`. One important thing is that we require `dflow(x, dx, t)` to return a Named Tuple: `(t = t, u = flow(x, p, t), du = dflow(x, p, dx, t))`, the last component being the value of the derivative of the flow.          |
| `evolve(fl, ::Val{:Full}, x, par, t; k...)`       |  The function implements the flow (or semigroup) associated to an autonomous Cauchy problem `(x, p, t) -> flow(x, p, t)`. The whole solution on the time interval [0,t] must be returned. It is not strictly necessary to provide this, it is mainly used for plotting on the user side. In the case of Poincaré Shooting, one must be able to call the flow like `evolve(fl, Val(:Full), x, par, Inf)`.          |
| `evolve(fl, ::Val{:SerialTimeSol}, x, par, t; k...)`       |  Serial version of the flow. Used for Matrix based jacobian (Shooting and Poincaré Shooting) and diffPoincareMap         |
| `evolve(fl, ::Val{:TimeSol}, x, par, t = Inf; k...)`       |  Flow which returns the tuple `(t, u(t))`. Optional, mainly used for plotting on the user side.       |
| `evolve(fl, ::Val{:SerialdFlow}, x, par, dx, tΣ; k...)`       |  Serial version of `dflow`. Used internally when using parallel multiple shooting. Named Tuple `(u = ..., du = ..., t = t)`. |
