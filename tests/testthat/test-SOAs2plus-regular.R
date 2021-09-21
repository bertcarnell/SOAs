test_that("SOAs2plus_regular", {
  ## s=2 follows a different approach than s>2
  temp <- SOAs2plus_regular(s=2, k=4, optimize=FALSE)
  expect_s3_class(temp, "SOA")
  expect_equal(dim(temp), c(16, 10))
  expect_equal(attr(temp, "type"), "SOA")
  expect_equal(attr(temp, "strength"), "2+")
  expect_equal(length(unique(c(temp))), 4)

  ## testing optimize
  set.seed(123)
  temp <- SOAs2plus_regular(s=3, k=3)
  expect_s3_class(temp, "SOA")
  expect_equal(dim(temp), c(27, 6))
  expect_equal(attr(temp, "type"), "OSOA")
  expect_equal(attr(temp, "strength"), "2+")
  expect_equal(length(unique(c(temp))), 9)

  temp <- SOAs2plus_regular(s=4, k=3, optimize=FALSE)
  expect_snapshot_output(temp)
  expect_s3_class(temp, "SOA")
  expect_equal(attr(temp, "type"), "OSOA")
  expect_equal(attr(temp, "strength"), "2+")
  expect_equal(dim(temp), c(64, (4^3-1)/(4-1) - ((4-1)^3-1)/(4-2)))
  expect_equal(length(unique(c(temp))), 16)

  expect_error(SOAs2plus_regular(s=4, k=3, m=9),
               regexp="m <= mbound is not TRUE", fixed=TRUE)
  expect_error(SOAs2plus_regular(s=2, k=3, optimize=FALSE),
               regexp="requires k >= 4", fixed=TRUE)
})
