## OGS Lambda2 Criterion

The **OGS Compute Lambda2 Criterion** computes the lambda2 turbulence criterion from a _vtkRectilinearGrid_ (structured grid) and generates the lambda2 value. The filter uses the same approach as the **OGS Derivatives** for the computation of the gradient.

The lambda2 criterion (Jeong & Hussain, 1995) is based on computing the value of the second eigenvalue (lambda2) of the tensor S2 + O2 where:

<center><img src="https://github.com/inogs/OGSParaviewSuite/blob/master/OGSPlugins/OGSLambda2Criterion/doc/eqs.png" alt="" width="200"/></center>

Properly tuned, this criterion yelds the same result as the Q criterion.

The code makes use of the _matrixMN class_ to deal with matrix operations and the computation of the eigenvalues. For the latter, a Jacobi iteration algorithm is used.