# Neimark-Sacker point

```@contents
Pages = ["ns.md"]
Depth = 2
```

At a Neimark-Sacker (NS) bifurcation of a periodic orbit $\gamma$ (with period $T$) for parameter value $p_0$ for the Cauchy problem 

$$\frac{du}{dt}=F(u,p),\tag{E}$$

the eigenvalues (Floquet coefficients) of the monodromy operator $\mathcal M=Y(T)$ solution to

$$\frac{dY}{dt}=A(t)Y(t), Y(0)=I_n$$

contain a pair of complex conjugate eigenvalues $e^{\pm i \theta}$ on the unit circle with $\theta\in(0, \pi)$, which satisfy the non-resonance conditions

$$e^{i q \theta}-1 \neq 0, \quad q=1,2,3,4 \text { (no strong resonances). }$$

At such a bifurcation point, an invariant torus of the system (E) branches off from the periodic orbit.

There are two ways to compute the normal form of this bifurcation:

1. using the Poincaré return map [^Kuznetsov],
2. using the method of [^Iooss], see also [^Kuz2].

You can obtain the normal form of a NS bifurcation using 

```
ns = get_normal_form(br, ind; prm = Val(false))
```

where `prm` indicates whether you want to use the method based on Poincaré return map (`Val(true)`) or the one based on Iooss method (`Val(false)`). The call returns a point of type `NeimarkSackerPO` from which you can access

- `ns.p` the parameter value at the bifurcation point,
- `ns.T` the period of the periodic orbit,
- `ns.ω` the frequency $\theta$,
- `ns.ζ`, `ns.ζ★` the right / left eigenvectors,
- `ns.nf.nf` the coefficients of the normal form.

You can pass `verbose = true` to get more information during the computation. Note that, for the [periodic orbits based on orthogonal collocation](@ref), you can skip the (costly) computation of the normal form and only get the bifurcation point characteristics by using `detailed = Val(false)`.

## Which method to use?

Depending on the method used for computing the periodic orbits, you have several possibilities:

- For shooting, you can only use the PRM method. Shooting is the preferred way for large scale systems. Note that the PRM method is not very precise numerically.
- For collocation, you can use PRM and Iooss methods. The Iooss method (`prm = Val(false)`, the default) is **the most precise**.
- For Trapezoid method, NS normal form is **not yet implemented**.

!!! note "Branch switching"
    The bifurcating object at a NS point is an invariant torus which cannot be represented by the periodic orbit methods. Hence, automatic branch switching (`continuation(br, ind)`) is **not available** from a NS point, contrary to the case of [Period-doubling point](@ref).

## Normal form based on Poincaré return map

Given a transversal section $\Sigma$ to $\gamma$ at $\gamma(0)$, the Poincaré return map $\mathcal P$ associates to each point $x\in\Sigma$ close to $\gamma(0)$ the first point $\mathcal P(x,p)\in\Sigma$ where the orbit of (E) with initial condition $x$ intersects again $\Sigma$. Hence, the discrete map $x_{n+1}=\mathcal P(x_n,p)$ has normal form

$$z_{n+1} = z_ne^{i\theta}(1+d|z_n|^2)$$

where[^Kuz2]

$$d=\frac{1}{2} e^{-i \theta}\left\langle \zeta^*, \mathcal{C}(\zeta, \zeta, \bar{\zeta})+2 \mathcal{B}\left(\zeta,\left(I_{n-1}-\mathcal{A}\right)^{-1} \mathcal{B}(\zeta, \bar{\zeta})\right)+\mathcal{B}\left(\bar{\zeta},\left(e^{2 i \theta} I_{n-1}-\mathcal{A}\right)^{-1} \mathcal{B}(\zeta, \zeta)\right)\right\rangle$$

where $\mathcal C=d_1^3\mathcal P(\gamma(0),p_0)$, $\mathcal B = d_1^2\mathcal P(\gamma(0),p_0)$ and $\mathcal A = d_1\mathcal P(\gamma(0),p_0)$. Also:

$$\mathcal{A} \zeta=e^{i \theta} \zeta, \mathcal{A}^{\mathrm{T}} \zeta^*=e^{-i \theta} \zeta^*, \text { and }\left\langle \zeta^*, \zeta\right\rangle=1$$

The coefficients `a` and `b` are returned in `ns.nf.nf`. The bifurcation is **supercritical** if $\Re(b) < 0$ and **subcritical** if $\Re(b) > 0$.

!!! danger "Large scale problems"
    The computation of the normal form is not optimized for Matrix-Free problems (e.g. Monodromy) yet.

!!! info "Collocation case"
    The monodromy matrix and other flow differentials are computed using finite differences.

## Normal form based on Iooss method

This is based on [^Iooss],[^Kuz2]. Suppose that the $T$ periodic orbit $\gamma(\tau)$ has a Neimark-Sacker bifurcation for a parameter value $p_0$. We also assume that there are no strong resonances.
Locally, the orbits can be represented by 

$$x(\tau) = \gamma(\tau)+Q_0(\tau)\xi+\Phi(\tau, \xi)$$

where 

$$\left\{\begin{aligned}
\frac{d \tau}{d t} & =1+a|\xi|^2+\cdots \\
\frac{d \xi}{d t} & =\frac{i \theta}{T} \xi+d \xi|\xi|^2+\cdots
\end{aligned}\right.$$

with center manifold correction $\Phi(\tau, \xi)$ being $T$ periodic in $\tau$ and $Q_0(\tau)$ built from the Floquet eigenvectors.

The coefficients `a` and `d` are returned in `ns.nf.nf`. The bifurcation is **supercritical** if $\Re(d) < 0$ and **subcritical** if $\Re(d) > 0$.

## See also

- [Detection of bifurcation points of periodic orbits](@ref)
- [Continuation of Neimark-Sacker (NS) bifurcations of periodic orbits](@ref)
- [Steinmetz-Larter model](@ref steinmetz): continuation of a curve of NS points of periodic orbits and computation of the NS normal form
- [Lorenz-84 model, take 2](@ref lorenz98-take2): computation of a curve of NS points of periodic orbits from a Hopf-Hopf bifurcation

## References

[^Kuznetsov]: > Yu. A. Kuznetsov, "Elements of Applied Bifurcation Theory", 2nd ed., 1998.

[^Kuz2]: > Kuznetsov et al., “Numerical Periodic Normalization for Codim 1 Bifurcations of Limit Cycles.”

[^Iooss]: > Iooss, "Global Characterization of the Normal Form for a Vector Field near a Closed Orbit.", 1988
