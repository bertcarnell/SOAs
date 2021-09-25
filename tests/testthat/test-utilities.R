test_that("utilities", {
  temp <- ff(3,3,3)
  expect_equal(dim(temp), c(27,3))

  temp <- interleavecols(matrix(2,3,3), matrix(7,3,3))
  expect_equal(temp[1,], c(2,7,2,7,2,7))

  expect_equal(mbound_LiuLiu(19, 3), 9)

  expect_error(ff("a", "b"),
               regexp="ff takes only integers as arguments",
               fixed=TRUE)
})
