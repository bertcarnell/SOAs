test_that("MDLE", {
  temp <- MDLEs(DoE.base::L16.4.5, 2, noptim.rounds = 1)

  expect_equal(class(temp), c("MDLE", "list"))
  expect_equal(names(temp), c("array", "phi_p", "optimized"))
  expect_equal(dim(temp$array), c(16, 5))
})
