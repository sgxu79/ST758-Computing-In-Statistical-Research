test_that("ridge_regression throws error when k is not a nonnegative number",{
  lambda = c(1,-1,0)
  X = matrix(rnorm(4),ncol=2)
  y = rnorm(2)
  expect_error(ridge_regression(y,X,lambda),"negative")
})

test_that("ridge_regression throws error when X and y are not conformable",{
  lambda = c(1,1,0)
  X = matrix(rnorm(6),ncol=3)
  y = rnorm(2)
  expect_error(ridge_regression(y,X,lambda),"non-conformable")
})

test_that("ridge_regression correctly calculates the ridge regression coefficients",{
  lambda = sample(1:10,1)
  X = matrix(rnorm(500),ncol=20)
  y = matrix(rnorm(25),ncol=1)
  I = diag(20)
  b_lam = ridge_regression(y,X,lambda)
  Z = (t(X)%*%X+lambda*I)%*%b_lam
  W = t(X)%*%y
  expect_equal(Z,W)
})
