test_that("createYcols", {
  temp <- createYcols(2)
  expect_true(all(vapply(temp, length, integer(1)) == 2^2-1))

  temp <- createYcols(5)
  expect_true(all(vapply(temp, length, integer(1)) == 2^5-1))
})

test_that("createABcols", {
  temp <- createABcols(4)
  expect_true(all(vapply(temp, length, integer(1)) == 5))

  temp <- createABcols(5)
  expect_true(all(vapply(temp, length, integer(1)) == 9))

  temp <- createABcols(6)
  expect_true(all(vapply(temp, length, integer(1)) == 20))

  temp <- createABcols(7)
  expect_true(all(vapply(temp, length, integer(1)) == 40))
})

test_that("create_ABC", {
  temp <- create_ABC(4)
  expect_equal(dim(temp$A), c(16, 4))
  expect_equal(dim(temp$B), c(16, 4))
  expect_equal(dim(temp$C), c(16, 4))
  expect_equal(length(temp$Yates.columns$A), 4)

  temp <- create_ABC(5)
  expect_equal(dim(temp$A), c(32, 8))

  temp <- create_ABC(5, m = 2)
  expect_equal(dim(temp$A), c(32, 2))

  temp <- create_ABC(5, m = 2, constr = "ShiTang_alpha")
  expect_equal(dim(temp$A), c(32, 2))
})

test_that("create_DfromABC", {
  temp <- create_DfromABC(create_ABC(4))
  expect_true(all(apply(temp, 2, sum) == sum(temp[,1])))

  temp <- create_DfromABC(create_ABC(7))
  expect_true(all(apply(temp, 2, sum) == sum(temp[,1])))

  temp <- create_DfromABC(create_ABC(9), random=TRUE)
  expect_equal(unname(c(table(temp))), rep(8192,8))

  temp <- create_DfromABC(create_ABC(9))
  expect_true(all(apply(temp, 2, sum) == sum(temp[,1])))

  temp <- create_DfromABC(create_ABC(4), random = FALSE)
  expect_true(all(apply(temp, 2, sum) == sum(temp[,1])))
})
