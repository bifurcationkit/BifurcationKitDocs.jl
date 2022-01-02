# Interface for Linear Solvers

The linear solver `ls` must be a subtype of the abstract type `AbstractLinearSolver`. It is then called as follows

| Required methods               |                        | Brief description                                                                     |
|:------------------------------ |:---------------------- |:------------------------------------------------------------------------------------- |
| `ls(J, rhs; a₀ = 0, a₁ = 1, kwargs...)`                |                        | Compute the solution `sol` of the linear problem `(a₀ * I + a₁ * J) * sol = rhs` where `J` is the jacobian. Returns `(sol, success::Bool, itnumber)` where `itnumber` is the number of iterations for solving the problem.|


!!! tip "Shifts"
    Two methods `_axpy` and `_axpy_op` are provided to help implementing this shift `a₀ * I + a₁ * J`

# Interface for Eigen Solvers

The linear solver `eig` must be a subtype of the abstract type `AbstractEigenSolver`. It is then called as follows

| Required methods               |                        | Brief description                                                                     |
|:------------------------------ |:---------------------- |:------------------------------------------------------------------------------------- |
| `eig(J, nev::Int; kwargs...)`                |                        | Compute the `nev` eigen-elements with largest real part. Returns `(eigenvalues, eigenvectors, success:Bool, itnumber)` where `itnumber` is the number of iterations for completing the computation. `kwargs` can be used to send information about the algorithm (perform bisection,...).|
| `geteigenvector(eig, eigenvectors, i::Int)`         |                        | Returns the ith eigenvectors from the set of `eigenvectors`.|

# Interface for Bordered Linear Solvers

The bordered linear solver `bls` must be a subtype of the abstract type `AbstractBorderedLinearSolver` (which is itself a subtype of `AbstractLinearSolver`). It is used to solve

$$\tag{BLS}\left[\begin{array}{ll}
{\text{shift}\cdot I+J} & {dR} \\
{(\xi_u \cdot dz_u)^T} & {\xi_p \cdot dz_p}
\end{array}\right] \cdot\left[\begin{array}{l}
{dX} \\
{dl}
\end{array}\right] = \left[\begin{array}{l}
R \\
n
\end{array}\right]$$

where $\xi_u,\xi_p\in\mathbb C$ and $dR,\xi_u\in\mathbb C^N$.

| Required methods               |                        | Brief description                                                                     |
|:------------------------------ |:---------------------- |:------------------------------------------------------------------------------------- |
| `bls(J, dR, dzu, dzp, R, n, ξu::Number, ξp::Number; shift = nothing, kwargs...)`                |                        | Compute the solution `dX, dl` of the linear problem (BLS) where `J` is the jacobian and `dR, dzu` are vectors (not necessarily subtypes of `AbstractVector`). `shift = nothing` is used in place of saying `shift=0`. Returns `(dX, dl, success::Bool, itnumber)` where `itnumber` is the number of iterations for solving the problem.|
