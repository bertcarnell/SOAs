test_that("SOAs", {
  temp <- SOAs(DoE.base::L64.4.6, optimize=FALSE)
  expect_s3_class(temp, "SOA")
  expect_equal(attr(temp, "type"), "SOA")
  expect_equal(attr(temp, "strength"), "3")
  expect_equal(dim(temp), c(64, 5))
  expect_equal(length(unique(c(temp))), 64)

  suppressMessages(temp <- SOAs(DoE.base::L64.4.6, t=2, noptim.rounds=2, dmethod = "euclidean"))
  expect_s3_class(temp, "SOA")
  expect_equal(attr(temp, "type"), "SOA")
  expect_equal(attr(temp, "strength"), "2")
  expect_equal(dim(temp), c(64, 6))
  expect_equal(length(unique(c(temp))), 16)

  temp <- SOAs(DoE.base::L81.3.5, t=4, optimize=FALSE)
  expect_equal(attr(temp, "type"), "SOA")
  expect_equal(attr(temp, "strength"), "4")
  expect_equal(dim(temp), c(81, 2))
  expect_equal(length(unique(c(temp))), 81)

  expect_error(SOAs(DoE.base::L9.3.4, t=3), regexp="requires a strength 3", fixed=TRUE)
})
