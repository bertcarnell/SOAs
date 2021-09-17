test_that("phi_optimize", {
  temp <- phi_optimize(DoE.base::L16.4.5)
  expect_equal(dim(temp), c(16, 5))

  ## check balance
  tab <- as.table(rep(20,4))
  names(tab) <- 0:3
  names(dimnames(tab)) <- "temp"
  expect_equal(table(temp), tab)

  ## check symmetry error
  expect_error(phi_optimize(DoE.base::L18.3.6.6.1), regexp="same number of levels", fixed=TRUE)
})
