
test_that("sweep_k and isweep_k undo each other",{
  n <- sample(10:200,1)
  u <- matrix(rnorm(n), ncol=1)
  A <- tcrossprod(u)
  diag(A) <- diag(A) + 1
  k = sample(1:n,1)
  expect_equal(A,isweep_k(sweep_k(A,k),k),tolerance = n*1e-14)
})

test_that("sweep and isweep undo each other for all k",{
  n <- sample(10:200,1)
  u <- matrix(rnorm(n), ncol=1)
  A <- tcrossprod(u)
  diag(A) <- diag(A) + 1
  expect_equal(A,isweep(sweep(A)),tolerance = n*1e-14)
})

test_that("sweep and isweep undo each other for some k",{
  n <- sample(10:200,1)
  u <- matrix(rnorm(n), ncol=1)
  A <- tcrossprod(u)
  diag(A) <- diag(A) + 1
  k = (1:n)[sample(1:n,sample(1:n,1))]
  expect_equal(A,isweep(sweep(A)),tolerance = n*1e-14)
})


test_that("sweep works correctly",{
  n <- sample(10:200,1)
  u <- matrix(rnorm(n), ncol=1)
  A <- tcrossprod(u)
  diag(A) <- diag(A) + 1
  B = sweep(A)
  I = diag(n)
  expect_equal(norm(A%*%B+I),0,tolerance = n*1e-14)
})


