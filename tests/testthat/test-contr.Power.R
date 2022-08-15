test_that("contr.Power", {
  expect_error(contr.Power(15, 2),
               regexp = "contr.Power requires that the number of levels is a power of s")
  expect_error(contr.Power(36, 6),
               regexp = "s must be a prime or a prime power")

  temp <- contr.Power(16, 2)
  expect_equal(dim(temp), c(16, 15))
  expect_equal(unname(colSums(temp)), rep(0L, 15))

  ## a variant that uses a Galois field
  temp <- contr.Power(16, 4)
  expect_equal(dim(temp), c(16, 15))
  expect_equal(unname(colSums(temp)), rep(0L, 15))

  ## invalid contrasts argument
  expect_error(contr.Power(4,2,contrasts=FALSE))

  ## invalid n
  expect_error(contr.Power(1,2))

  ## vector-valued n
  temp <- contr.Power(1:4,2)
  expect_equal(dim(temp), c(4, 3))

  ## n=s
  temp <- contr.Power(9,9)
  expect_equal(dim(temp), c(9,8))

})
