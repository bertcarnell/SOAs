test_that("MDLE", {
  suppressMessages(temp <- MDLEs(DoE.base::L16.4.5, 2, noptim.rounds = 1))

  expect_s3_class(temp, "MDLE")
  expect_equal(c("phi_p", "type", "optimized", "call") %in% names(attributes(temp)), rep(TRUE, 4))
  expect_equal(dim(temp), c(16, 5))
  expect_error(suppressMessages(MDLEs(DoE.base::L16.4.5, 3, noptim.rounds = 1), regexp="(n/(s * ell))%%1 == 0", fixed=TRUE))
  expect_error(suppressMessages(MDLEs(DoE.base::L72.2.6.3.7.4.1.6.6, el=3)))

  # operational checks for code coverage
  suppressMessages(temp <- MDLEs(DoE.base::L16.4.5, 2, noptim.rounds = 1,
                                 noptim.oa = 2, storeperms = TRUE))
  expect_s3_class(temp, "MDLE")
  suppressMessages(temp <- MDLEs(DoE.base::L16.4.5, 2, optimize = FALSE))
  expect_s3_class(temp, "MDLE")
})
