test_that("utilities", {
  temp <- ff(3,3,3)
  expect_equal(dim(temp), c(27,3))

  temp <- interleavecols(matrix(2,3,3), matrix(7,3,3))
  expect_equal(temp[1,], c(2,7,2,7,2,7))

})
