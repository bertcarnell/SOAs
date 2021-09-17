test_that("MDLE", {
  temp <- MDLEs(DoE.base::L16.4.5, 2, noptim.rounds = 1)

  expect_s3_class(temp, "MDLE")
  expect_equal(c("phi_p", "type", "optimized", "call") %in% names(attributes(temp)), rep(TRUE, 4))
  expect_equal(dim(temp), c(16, 5))
  expect_error(MDLEs(DoE.base::L16.4.5, 3, noptim.rounds = 1), regexp="(n/(s * ell))%%1 == 0", fixed=TRUE)
  expect_error(MDLEs(DoE.base::L72.2.6.3.7.4.1.6.6, el=3))
})
