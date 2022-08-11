test_that("contr.Power", {
  expect_error(contr.Power(15, 2),
               regexp = "contr.Power requires that the number of levels is a power of s")
  temp <- contr.Power(16, 2)
  expect_equal(dim(temp), c(16, 15))
  expect_true(all(colSums(round(temp,8))==0))

  ## a variant that uses a Galois field
  temp <- contr.Power(16, 4)
  expect_equal(dim(temp), c(16, 15))
  expect_true(all(colSums(round(temp,8))==0))

  expect_error(contr.Power(4,2,contrasts=FALSE))

})
