test_that("OSOAarbitrary", {
  #temp <- OSOAarbitrary(DoE.base::L16.4.5) # Errors.  Is this a bug?
  temp <- OSOAarbitrary(DoE.base::L16.4.5, m = 2)

  expect_equal(dim(temp), c(64, 2))
})
