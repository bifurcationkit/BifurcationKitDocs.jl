# Interface for Linear Solvers

The linear solver `ls` must be a subtype of the abstract type `AbstractLinearSolver`. It is then called as follows

| Required methods               |                        | Brief description                                                                     |
|:------------------------------ |:---------------------- |:------------------------------------------------------------------------------------- |
| `ls(J, rhs)`                |                        | Compute the solution `sol` of the linear problem `J * sol = rhs` where `J` is the jacobian. Returns `(sol, success::Bool, itnumber)` where `itnumber` is the number of iterations for solving the problem.|

# Interface for Eigen Solvers

The linear solver `eig` must be a subtype of the abstract type `AbstractEigenSolver`. It is then called as follows

| Required methods               |                        | Brief description                                                                     |
|:------------------------------ |:---------------------- |:------------------------------------------------------------------------------------- |
| `eig(J, nev::Int; kwargs...)`                |                        | Compute the `nev` eigen-elements with largest real part. Returns `(eigenvalues, eigenvectors, success:Bool, itnumber)` where `itnumber` is the number of iterations for completing the computation. `kwargs` can be used to send information about the algorithm (perform bisection,...).|
| `geteigenvector(eig, eigenvectors, i::Int)`         |                        | Returns the ith eigenvectors from the set of `eigenvectors`.|
