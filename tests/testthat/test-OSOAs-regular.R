test_that("OSOAs_regular", {
  temp <- OSOAs_regular(s=3, k=3, el=3, optimize=FALSE)
  expect_s3_class(temp, "SOA")
  expect_equal(attr(temp, "type"), "OSOA")
  expect_equal(attr(temp, "strength"), "2*")
  expect_equal(dim(temp), c(27, 4))
  expect_equal(length(unique(c(temp))), 27)

  expect_snapshot_output((temp <- OSOAs_regular(s=3, k=3, el=2, optimize=FALSE)))
  expect_s3_class(temp, "SOA")
  expect_equal(attr(temp, "type"), "OSOA")
  expect_equal(attr(temp, "strength"), "2+")
  expect_equal(dim(temp), c(27, 4))
  expect_equal(length(unique(c(temp))), 9)

  expect_error(OSOAs_regular(s=3, k=3, el=2, m=5),
               regexp="m is too large", fixed=TRUE)

  expect_error(OSOAs_regular(s = 4, k = 3, el = 3, m = 5),
               regexp="m is too large in combination with el=3")

  temp <- OSOAs_regular(s = 3, k = 3, m = 3, el = 2, optimize = FALSE)
  expect_equal(dim(temp), c(27, 3))

  temp <- OSOAs_regular(s = 3, k = 3, m = 3, el = 3, optimize = FALSE)
  expect_equal(dim(temp), c(27, 3))
})
