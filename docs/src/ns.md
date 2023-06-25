# Neimark-Sacker point

At a Neimark-Sacker (NS) bifurcation of a periodic orbit $\gamma$ (with period $T$) for parameter value $p_0$ for the Cauchy problem 

$$\frac{du}{dt}=F(u,p),\tag{E}$$

the eigenvalues (Floquet coefficients) of the monodromy operator $\mathcal M=Y(T)$ solution to

$$\frac{dY}{dt}=A(t)Y(t), Y(0)=I_n$$

contains the eigenvalues $e^{\pm i \theta}$ with $\theta$ and 

$$e^{i q \theta}-1 \neq 0, \quad q=1,2,3,4 \text { (no strong resonances). }$$

There are two ways to compute the normal form of this bifurcation

1. using the Poincaré return map [^Kuznetsov]
2. using the method of [^Iooss], see also [^Kuz2]

You can obtain the normal form of a NS bifurcation using 

```
pd = getNormalForm(br, ind)
```


## Normal form based on Poincaré return map

Given a transversal section $\Sigma$ to $\gamma$ at $\gamma(0)$, the Poincaré return map $\mathcal P$ associates to each point $x\in\Sigma$ close to $\gamma(0)$ the point $\mathcal P(x,p)\in\Sigma$ where the orbit of (E) with initial condition $x$ intersects again $\Sigma$ at $\mathcal P(x,p)$. Hence, the discrete map $x_{n+1}=\mathcal P(x_n,p)$ has normal form

$$z_{n+1} = z_ne^{i\theta}(1+d|z_n|^2)$$

where[^Kuz2]

$$d=\frac{1}{2} e^{-i \theta}\left\langle v^*, \mathcal{C}(v, v, \bar{v})+2 \mathcal{B}\left(v,\left(I_{n-1}-\mathcal{A}\right)^{-1} \mathcal{B}(v, \bar{v})\right)+\mathcal{B}\left(\bar{v},\left(e^{2 i \theta} I_{n-1}-\mathcal{A}\right)^{-1} \mathcal{B}(v, v)\right)\right\rangle$$

where $\mathcal C=d^3\mathcal P(\gamma(0))$, $\mathcal B = d^2\mathcal P(\gamma(0))$ and $\mathcal A = d\mathcal P(\gamma(0))$. Also:

$$\mathcal{A} v=e^{i \theta} v, \mathcal{A}^{\mathrm{T}} v^*=e^{-i \theta} v^*, \text { and }\left\langle v^*, v\right\rangle=1$$

## Normal form based on Iooss method

Not implemented.

## References

[^Kuznetsov] :> Yu. A. Kuznetsov, "Elements of Applied Bifurcation Theory", 2nd ed., 1998.

[^Kuz2] :> Kuznetsov et al., “Numerical Periodic Normalization for Codim 1 Bifurcations of Limit Cycles.”

[^Iooss] :> Iooss, "Global Characterization of the Normal Form for a Vector Field near a Closed Orbit.", 1988