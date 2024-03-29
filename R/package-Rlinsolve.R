#' A Collection of Iterative Solvers for (Sparse) Linear System of Equations
#'
#' Solving a system of linear equations is one of the most fundamental
#' computational problems for many fields of mathematical studies, such as
#' regression from statistics or numerical partial differential equations.
#' We provide a list of both stationary and nonstationary solvers. Sparse
#' matrix class from \pkg{Matrix} is also supported for large sparse system.
#'
#' @section Non-square matrix:
#' For a matrix \eqn{A} of size \code{(m-by-n)}, we say the system is
#' \strong{overdetermined} if \code{m>n}, \strong{underdetermined} if \code{m<n}, or
#' \strong{squared} if \code{m=n}. In the current version, underdetermined system is
#' not supported; it will later appear with sparse constraints. For an overdetermined system,
#' it automatically transforms the problem into \emph{normal equation}, i.e.,
#' \deqn{Ax=b \rightarrow A^T Ax = A^T b}
#' even though if suffers from worse condition number while having desirable property
#' of a system to be symmetric and positive definite.
#'
#' @section Sparsity:
#' \pkg{RcppArmadillo} is extensively used in the package. In order for
#' bullet-proof transition between dense and sparse matrix, only 3 of
#' 12 RcppArmadillo-supported sparse matrix formats have access to
#' our algorithms; \code{"dgCMatrix"},\code{"dtCMatrix"} and \code{"dsCMatrix"}.
#' Please see \href{https://CRAN.R-project.org/package=RcppArmadillo/vignettes/RcppArmadillo-sparseMatrix.pdf}{the vignette}
#' on sparse matrix support from RcppArmadillo. If either of two inputs \code{A} or \code{b} is
#' sparse, all matrices involved are automatically transformed into sparse matrices.
#'
#'
#' @section Composition of the Package:
#' Following is a list of stationary methods,
#' \describe{
#'   \item{\code{\link{lsolve.jacobi}}}{Jacobi method}
#'   \item{\code{\link{lsolve.gs}}}{Gauss-Seidel method}
#'   \item{\code{\link{lsolve.sor}}}{Successive Over-Relaxation method}
#'   \item{\code{\link{lsolve.ssor}}}{Symmetric Successive Over-Relaxation method}
#' } as well as nonstationary (or, Krylov subspace) methods,
#' \describe{
#'   \item{\code{\link{lsolve.bicg}}}{Bi-Conjugate Gradient method}
#'   \item{\code{\link{lsolve.bicgstab}}}{Bi-Conjugate Gradient Stabilized method}
#'   \item{\code{\link{lsolve.cg}}}{Conjugate Gradient method}
#'   \item{\code{\link{lsolve.cgs}}}{Conjugate Gradient Squared method}
#'   \item{\code{\link{lsolve.cheby}}}{Chebyshev method}
#'   \item{\code{\link{lsolve.gmres}}}{Generalized Minimal Residual method}
#'   \item{\code{\link{lsolve.qmr}}}{Quasi-Minimal Residual method}
#' }
#' Also, \code{\link{aux.fisch}} is provided to generate a sparse system of
#' discrete Poisson matrix from finite difference approximation scheme of Poisson equation
#' on 2-dimensional square domain.
#'
#' @references Demmel, J.W. (1997) \emph{Applied Numerical Linear Algebra, 1st ed.}, SIAM.
#' @references Barrett, R., Berry, M., Chan, T.F., Demmel, J., Donato, J., Dongarra, J.,
#' Eijkhout, V., Pozo, R., Romine, C., and van der Vorst, H. (1994) \emph{Templates for the Solution
#' of Linear Systems: Building Blocks for Iterative Methods, 2nd ed.} Philadelphia, SIAM.
#' 
#' @docType package
#' @name Rlinsolve
#' @import Matrix
#' @import Rdpack
#' @importFrom utils packageVersion
#' @importFrom stats rnorm
#' @importFrom Rcpp evalCpp
#' @useDynLib Rlinsolve, .registration=TRUE
NULL


