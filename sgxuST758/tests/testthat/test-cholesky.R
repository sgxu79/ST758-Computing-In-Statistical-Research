test_that("check_k_banded throws error when k is not a nonnegative integer",{
  m = 5
  k = sample(c(-1,1.5,-3.5),1)
  u <- matrix(rnorm(m), ncol=1)
  A <- tcrossprod(u)
  diag(A) <- diag(A) + 1
  A <- round(10*A)/10
  for (i in 1:m) {
    for (j in 1:m) {
      if (abs(i - j) > k) A[i,j] <- 0
    }
  }
  expect_error(check_k_banded(A,k),"k should be a nonnegative integer")
})

test_that("check_k_banded throws error when A is not a symmetric matrix",{
  m = 5
  k = 2
  A = matrix(rep(1:25),nrow=5)
  for (i in 1:m) {
    for (j in 1:m) {
      if (abs(i - j) > k) A[i,j] <- 0
    }
  }
  expect_error(check_k_banded(A,k),"A is not symmetric")
})

test_that("chol_banded throws error when A is not k-banded",{
  m = 5
  k = 2
  A = matrix(rep(2,25),nrow=5)
  for (i in 1:m) {
    for (j in 1:m) {
      if (abs(i - j) > k) A[i,j] <- 0
    }
  }
  expect_error(chol_banded(A,k-1),"Matrix is not k banded")
})

test_that("chol_banded correctly factorizes k-banded matrix",{
  m <- sample(1:50,1)
  k <- sample(1:round(m/2),1)
  u <- matrix(rgamma(m,5,5)+1, ncol=1)
  A <- tcrossprod(u)
  A[lower.tri(A)] <- 0
  for (i in 1:m) {
    for (j in 1:m) {
      if (abs(i - j) > k) A[i,j] <- 0
    }
  }
  A <- round(10*A)/10
  A = t(A)%*%A
  Rtest = chol_banded(A,k,TRUE)
  Rtrue = as(chol(A),"dgCMatrix")
  expect_equal(Rtest,Rtrue,tolerance = m*1e-14,check.attributes = FALSE)
})

test_that("check_k_banded checks correctly whether a matrix is k-banded or not",{
  m <- sample(1:100,1)
  k <- sample(1:m,1)
  l <- sample(1:(k-1),1)
  A <- matrix(rnorm(m**2), m, m)
  A <- A + t(A)
  A <- round(10*A)/10
  for (i in 1:m) {
    for (j in 1:m) {
      if (abs(i - j) > k) A[i,j] <- 0
    }
  }
  expect_equal(check_k_banded(A,k),TRUE)
  expect_equal(check_k_banded(A,l),FALSE)
})
