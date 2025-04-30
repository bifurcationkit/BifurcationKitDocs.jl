# Period-doubling point

```@contents
Pages = ["pd.md"]
Depth = 2
```

At a period-doubling (PD) bifurcation of a periodic orbit $\gamma$ (with period $T$) for parameter value $p_0$ for the Cauchy problem 

$$\frac{du}{dt}=F(u,p),\tag{E}$$

the eigenvalues (Floquet coefficients) of the monodromy operator $\mathcal M=Y(T)$ solution to

$$\frac{dY}{dt}=A(t)Y(t), Y(0)=I_n$$

contain the simple eigenvalue $\mu=-1$.

There are two ways to compute the normal form of this bifurcation

1. using the Poincaré return map [^Kuznetsov]
2. using the method of [^Iooss] see also [^Kuz2]

You can obtain the normal form of a PD bifurcation using 

```
pd = get_normal_form(br, ind; prm = false)
```

where `prm` indicates whether you want to use the method based on Poincaré return map (PRM) or the one based on Iooss method.

## Which method to use?

Depending on the method used for computing the periodic orbits, you have several possibilities:

- For shooting, you can only the PRM method. Shooting is the preferred way for large scale systems. Note that the PRM method is not very precise numerically.
- For collocation, you can use PRM and Iooss methods. Note that the Iooss method is the most precise.
- For Trapezoid method, PD normal form is not yet implemented.

## Predictor

The predictor for a non trivial guess at distance $\delta p$ from the bifurcation point is provided by the method

```@docs
predictor(hp::BifurcationKit.PeriodDoublingPO{ <: ShootingProblem }, δp, 
                    ampfactor; 
                    override = false)
```

If `override = true`, then the predictor is simply `x0 .+ ampfactor .* e` for the parameter `p0 + δp` and where `e` is a bifurcating eigenvector.

## Normal form based on Poincaré return map

Given a transversal section $\Sigma$ to $x_0$ at $x_0(0)$, the Poincaré return map $\mathcal P$ associates to each point $x\in\Sigma$ close to $x_0(0)$ the first point $\mathcal P(x,p)\in\Sigma$ where the orbit of (E) with initial condition $x$ intersects again $\Sigma$ at $\mathcal P(x,p)$. Hence, the discrete map $x_{n+1}=\mathcal P(x_n,p)$ has normal form

$$x_{n+1} = -x_n+cx_n^3+...$$

where [^Kuz2]

$$c =\frac{1}{6}\left\langle p^*, \mathcal{C}(p, p, p)+3 \mathcal{B}\left(p,\left(I_{n-1}-\mathcal{A}\right)^{-1} \mathcal{B}(p, p)\right)\right\rangle$$

where $\mathcal C=d^3\mathcal P(x_0(0))$, $\mathcal B = d^2\mathcal P(x_0(0))$ and $\mathcal A = d\mathcal P(x_0(0))$. Also:

$$\mathcal{A} p=-p, \mathcal{A}^{\mathrm{T}} p^*=-p^*$$

## Normal form based on Iooss method

This is based on [^Iooss],[^Kuz2]. Suppose that the $T$ periodic orbit $x_0(\tau)$ has a Period-Doubling bifurcation for a parameter value $p_0$.
Locally, the orbits can be represented by $p-p_0:=\mu$ and

$$x(\tau) = x_0(\tau)+\xi v(\tau)+H(\tau, \xi, \mu)$$

where 

$$\left\{\begin{array}{l}
\frac{d \tau}{d t}=1+a_{01}\cdot(p-p_0)+a_2 \xi^2+\cdots \\
\frac{d \xi}{d \tau}=c_{11}\cdot(p-p_0)\xi+c_3 \xi^3+\cdots
\end{array}\right.$$

with center manifold correction $H(\tau, \xi, \mu)$ being $2T$ periodic in $\tau$ and $v(\tau)$ is a Floquet eigenvector for the eigenvalue -1.


## References

[^Kuznetsov]: > Yu. A. Kuznetsov, "Elements of Applied Bifurcation Theory", 2nd ed., 1998.

[^Kuz2]: > Kuznetsov et al., “Numerical Periodic Normalization for Codim 1 Bifurcations of Limit Cycles.”

[^Iooss]: > Iooss, "Global Characterization of the Normal Form for a Vector Field near a Closed Orbit.", 1988