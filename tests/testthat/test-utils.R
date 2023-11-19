test_that("tr_estimate() works", {
  X <- matrix(rnorm(100), 20, 5)
  expect_equal(tr_estimate_R_impl(X), tr_estimate(X))
})
