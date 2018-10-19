#' Check Matrix is k-banded
#'
#' \code{check_k_banded} returns a Boolean variable indicating whether or
#' not a matrix is k-banded.
#'
#' @param A The input matrix
#' @param k bandwidth parameter
#' @export
check_k_banded <- function(A, k) {
  if(!((k%%1==0) && (k>=0)) ) stop("k should be a nonnegative integer")
  size = dim(A)[1]
  if(is.null(size)){
    A = as.matrix(A)
    size = dim(A)[1]
  }
  if(all.equal(A,t(A))!=TRUE) stop("A is not symmetric")
  boolean = TRUE
  if(k<(size-1)){
    for(i in 1:(size-k-1)){
      for(j in (i+k+1):size){
        if(A[i,j]!=0){
          boolean = FALSE
          break
        }
      }
    }
  }
  return(boolean)
}


#' Banded Cholesky
#'
#' \code{chol_banded} computes the Cholesky decomposition of a banded
#' positive definite matrix.
#'
#' @param A The input matrix
#' @param k bandwidth parameter
#' @param checks_off Boolean variable to turn on / off checks
#' @import Matrix
#' @export
chol_banded <- function(A, k, checks_off=FALSE) {
  if(checks_off==FALSE){
    if(!check_k_banded(A,k)) stop("Matrix is not k banded")
  }
  size = dim(A)[1]
  if(is.null(size)){
    A = as.matrix(A)
    size = dim(A)[1]
  }
  A[lower.tri(A)] <- 0
  A = Matrix(A,sparse = TRUE)
  m = 1
  while(m<size){
    if(A[m,m]<=0) stop("Matrix not positive definite")
    rho = sqrt(A[m,m])
    A[m,m] = rho
    for(j in (m+1):min(size,m+k)){
      A[m,j] = A[m,j]/rho
      for(i in (m+1):j){
        A[i,j] = A[i,j] - A[m,i]*A[m,j]
      }
    }
    m = m+1
  }
  if(A[m,m]<=0) stop("Matrix not positive definite")
  A[m,m] = sqrt(A[m,m])
  return(A)
}

