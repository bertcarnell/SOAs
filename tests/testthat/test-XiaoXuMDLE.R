test_that("XiaoXuMDLE", {
  suppressMessages(temp <- XiaoXuMDLE(DoE.base::L16.4.5, 2, noptim.oa=5,
                                      nrounds = 5, nsteps=50))

  expect_true("phi_p" %in% names(attributes(temp)))
  expect_equal(dim(temp), c(16, 5))

  # operational checks for code coverage
  suppressMessages(temp <- XiaoXuMDLE(DoE.base::L16.4.5, 2, nrounds = 5, nsteps=50,
                                      dmethod="euclidean", p=20))

})
