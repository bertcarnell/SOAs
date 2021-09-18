test_that("OSOAarbitrary", {
  temp <- OSOAarbitrary(DoE.base::L16.4.5)
  expect_equal(dim(temp), c(64, 4))
  temp <- OSOAarbitrary(DoE.base::L16.4.5, el=2)
  expect_equal(dim(temp), c(64, 5))
  temp <- OSOAarbitrary(DoE.base::L16.4.5, m = 2)
  expect_equal(dim(temp), c(64, 2))
  expect_error(OSOAarbitrary(DoE.base::L16.4.5, m = 5), regexp="at most 4 columns", fixed=TRUE)

  # test the oa as a data.frame
  temp <- OSOAarbitrary(as.data.frame(DoE.base::L16.4.5))
  expect_equal(dim(temp), c(64, 4))
  expect_equal(dim(attr(temp, "A")), c(64, 4))
})
