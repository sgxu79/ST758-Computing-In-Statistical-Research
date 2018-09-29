

#' Sweep k
#'
#' \code{sweep_k} applies the sweep operator to a symmetric matrix
#' on the kth diagonal entry if it is possible.
#'
#' @param A The input symmetric matrix
#' @param k Diagonal index on which to sweep
#' @export
sweep_k <- function(A, k) {
  if(all.equal(A,t(A))!=TRUE) stop("Matrix not symmetric")
  if(A[k,k]<=0) stop("kth diagonal entry not positive")
  size = dim(A)[1]
  akk = A[k,k]
  A[lower.tri(A)] = 0
  if(k>1){
    sw_col = A[1:(k-1),k]
    sw_col_n = sw_col/akk
    A[1:(k-1),k] = sw_col/akk
  }
  if(k<size){
    sw_row = A[k,(k+1):size]
    sw_row_n = sw_row/akk
    A[k,(k+1):size] = sw_row/akk
  }
  A[k,k] = -1/akk
  r_r = (1:size)[-k]
  r_c = r_r
  for(i in r_r){
    for(j in r_c){
      if(i>k){aik = sw_row_n[i-k]}else{aik = sw_col_n[i] }
      if(j<k){akj = sw_col_n[j]}else{akj = sw_row_n[j-k]}
      A[i,j] = A[i,j] - (aik*akj*akk)
    }
    r_c = r_c[!r_c %in% i]
  }
  d = diag(A)
  return(A+t(A)-diag(d))
}



#' Inverse Sweep k
#'
#' \code{isweep_k} applies the inverse sweep operator to a symmetric matrix
#' on the kth diagonal entry if it is possible.
#'
#' @param A The input symmetric matrix
#' @param k Diagonal index on which to sweep
#' @export
isweep_k <- function(A, k) {
  if(all.equal(A,t(A))!=TRUE) stop("Matrix not symmetric")
  if(A[k,k]>=0) stop("kth diagonal entry not negative")
  A = -1*A
  return(-1*sweep_k(A,k))
}



#' Sweep
#'
#' @param A The input symmetric matrix
#' @param k Diagonal index entry set on which to sweep
#' @export
sweep <- function(A, k=NULL) {
  if(all.equal(A,t(A))!=TRUE) stop("Matrix not symmetric")
  if(is.null(k)){
    if(!all(diag(A)>0)) stop("Matrix is not positive definite")
    size = dim(A)[1]
    for(i in 1:size){
      A = sweep_k(A,i)
    }
    return(A)
  }else{
    for(i in k){
      A = sweep_k(A,i)
    }
    return(A)
  }
}



#' Inverse Sweep
#'
#' @param A The input symmetric matrix
#' @param k Diagonal index entry set on which to sweep
#' @export
isweep <- function(A, k=NULL) {
  if(all.equal(A,t(A))!=TRUE) stop("Matrix not symmetric")
  if(is.null(k)){
    if(!all(diag(A)<0)) stop("Matrix is not negative definite")
    size = dim(A)[1]
    for(i in 1:size){
      A = isweep_k(A,i)
    }
    return(A)
  }else{
    for(i in k){
      A = isweep_k(A,i)
    }
    return(A)
  }
}



