# Neimark-Sacker point

```@contents
Pages = ["ns.md"]
Depth = 2
```

At a Neimark-Sacker (NS) bifurcation of a periodic orbit $\gamma$ (with period $T$) for parameter value $p_0$ for the Cauchy problem 

$$\frac{du}{dt}=F(u,p),\tag{E}$$

the eigenvalues (Floquet coefficients) of the monodromy operator $\mathcal M=Y(T)$ solution to

$$\frac{dY}{dt}=A(t)Y(t), Y(0)=I_n$$

contain the eigenvalues $e^{\pm i \theta}$ with $\theta$ and 

$$e^{i q \theta}-1 \neq 0, \quad q=1,2,3,4 \text { (no strong resonances). }$$

There are two ways to compute the normal form of this bifurcation

1. using the Poincaré return map [^Kuznetsov]
2. using the method of [^Iooss], see also [^Kuz2]

You can obtain the normal form of a NS bifurcation using 

```
pd = get_normal_form(br, ind; prm = false)
```

where `prm` indicates whether you want to use the method based on Poincaré return map (PRM) or the one based on Iooss method.

## Which method to use?

Depending on the method used for computing the periodic orbits, you have several possibilities:

- For shooting, you can only the PRM method. Shooting is the preferred way for large scale systems. Note that the PRM method is not very precise numerically.
- For collocation, you can use PRM and Iooss methods. Note that the Iooss method is **the most precise**.
- For Trapezoid method, NS normal form is **not yet implemented**.

## Normal form based on Poincaré return map

Given a transversal section $\Sigma$ to $\gamma$ at $\gamma(0)$, the Poincaré return map $\mathcal P$ associates to each point $x\in\Sigma$ close to $\gamma(0)$ the first point $\mathcal P(x,p)\in\Sigma$ where the orbit of (E) with initial condition $x$ intersects again $\Sigma$ at $\mathcal P(x,p)$. Hence, the discrete map $x_{n+1}=\mathcal P(x_n,p)$ has normal form

$$z_{n+1} = z_ne^{i\theta}(1+d|z_n|^2)$$

where[^Kuz2]

$$d=\frac{1}{2} e^{-i \theta}\left\langle v^*, \mathcal{C}(v, v, \bar{v})+2 \mathcal{B}\left(v,\left(I_{n-1}-\mathcal{A}\right)^{-1} \mathcal{B}(v, \bar{v})\right)+\mathcal{B}\left(\bar{v},\left(e^{2 i \theta} I_{n-1}-\mathcal{A}\right)^{-1} \mathcal{B}(v, v)\right)\right\rangle$$

where $\mathcal C=d^3\mathcal P(\gamma(0))$, $\mathcal B = d^2\mathcal P(\gamma(0))$ and $\mathcal A = d\mathcal P(\gamma(0))$. Also:

$$\mathcal{A} v=e^{i \theta} v, \mathcal{A}^{\mathrm{T}} v^*=e^{-i \theta} v^*, \text { and }\left\langle v^*, v\right\rangle=1$$

## Normal form based on Iooss method

This is based on [^Iooss],[^Kuz2]. Suppose that the $T$ periodic orbit $\gamma(\tau)$ has a Neimark-Sacker bifurcation for a parameter value $p_0$. We also assume that there are no strong resonances.
Locally, the orbits can be represented by 

$$x(\tau) = \gamma(\tau)+Q_0(\tau)\xi+\Phi(\tau, \xi)$$

where 

$$\left\{\begin{aligned}
\frac{d \tau}{d t} & =1+a|\xi|^2+\cdots \\
\frac{d \xi}{d t} & =\frac{i \theta}{T} \xi+d \xi|\xi|^2+\cdots
\end{aligned}\right.$$

with center manifold correction $\Phi(\tau, \xi)$ being $T$ periodic in $\tau$ and $Q_0(\tau)$ is the Floquet operator.


## References

[^Kuznetsov]: > Yu. A. Kuznetsov, "Elements of Applied Bifurcation Theory", 2nd ed., 1998.

[^Kuz2]: > Kuznetsov et al., “Numerical Periodic Normalization for Codim 1 Bifurcations of Limit Cycles.”

[^Iooss]: > Iooss, "Global Characterization of the Normal Form for a Vector Field near a Closed Orbit.", 1988