# Branch point of periodic orbit

```@contents
Pages = ["bppo.md"]
Depth = 3
```

At a Branch point (BP) of a periodic orbit $\gamma$ (with period $T$)  for the Cauchy problem at parameter value $p_0$:

$$\frac{du}{dt}=F(u,p),\tag{E}$$

the eigenvalues (Floquet coefficients) of the monodromy operator $\mathcal M=Y(T)$ solution to

$$\frac{dY}{dt}=A(t)Y(t), Y(0)=I_n$$

contain the eigenvalue $1$ with algebraic multiplicity 2.
There are two ways to compute the normal form of this bifurcation

1. using the Poincaré return map [^Kuznetsov]
2. using the method of [^Iooss], see also [^Kuz2]

You can obtain the normal form of a BP bifurcation using 

```
pd = get_normal_form(br, ind; prm = false)
```

where `prm` indicates whether you want to use the method based on Poincaré return map (PRM) or the one based on Iooss method.

## Predictor

The predictor for a non trivial guess at distance $\delta p$ from the bifurcation point is provided by the method

```@docs
predictor(hp::BifurcationKit.BranchPointPO{ <: ShootingProblem }, δp, 
                    ampfactor; 
                    override = false)
```

## Which method to use?

Depending on the method used for computing the periodic orbits, you have several possibilities:

- For shooting, you can only the PRM method. Shooting is the preferred way for large scale systems. Note that the PRM method is not very precise numerically.
- For collocation, you can use PRM and Iooss methods. Note that the Iooss method is the most precise. This is not yet implemented.
- For Trapezoid method, the normal form is not yet implemented.

## Normal form based on Poincaré return map

Given a transversal section $\Sigma$ to $\gamma$ at $\gamma(0)$, the Poincaré return map $\mathcal P$ associates to each point $x\in\Sigma$ close to $\gamma(0)$ the first point $\mathcal P(x,p)\in\Sigma$ where the orbit of (E) with initial condition $x$ intersects again $\Sigma$ at $\mathcal P(x,p)$. Hence, the discrete map $x_{n+1}=\mathcal P(x_n,p)$ has normal form

$$x_{n+1} = x_n+a_{10}\delta p + a_{11}\delta p x+a_{02}x^2+a_{03}x^3$$
## References

[^Kuznetsov]: > Yu. A. Kuznetsov, "Elements of Applied Bifurcation Theory", 2nd ed., 1998.

[^Kuz2]: > Kuznetsov et al., “Numerical Periodic Normalization for Codim 1 Bifurcations of Limit Cycles.”

[^Iooss]: > Iooss, "Global Characterization of the Normal Form for a Vector Field near a Closed Orbit.", 1988

