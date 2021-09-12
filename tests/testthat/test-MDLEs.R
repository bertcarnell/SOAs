test_that("MDLE", {
  temp <- MDLEs(DoE.base::L16.4.5, 2, noptim.rounds = 1)

  expect_equal(class(temp)[1], "MDLE")
  expect_equal(c("phi_p", "type", "optimized", "call") %in% names(attributes(temp)), rep(TRUE, 4))
  expect_equal(dim(temp), c(16, 5))
})
