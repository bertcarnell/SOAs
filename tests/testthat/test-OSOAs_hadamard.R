test_that("OSOAs_hadamard", {
  suppressMessages(temp <- OSOAs_hadamard(m = 6))
  expect_s3_class(temp, "SOA")
  expect_equal(dim(temp), c(16, 6))
  expect_equal(attr(temp, "type"), "OSOA")

  ## 9 8-level factors in 24 runs
  temp <- OSOAs_hadamard(9, optimize=FALSE)
  expect_snapshot_output(temp)
  expect_equal(dim(temp), c(24, 9))
  expect_equal(length(unique(c(temp))), 8)

  ## 6 8-level factors in 16 runs
  temp <- OSOAs_hadamard(n=16, optimize=FALSE)
  expect_equal(dim(temp), c(16, 6))
  expect_equal(length(unique(c(temp))), 8)

  ## 6 4-level factors in 24 runs
  temp <- OSOAs_hadamard(n=24, m=6, el=2, optimize=FALSE)
  expect_equal(dim(temp), c(24, 6))
  expect_equal(length(unique(c(temp))), 4)

  ## 35 8-level factors in 80 runs
  temp <- OSOAs_hadamard(m=35, optimize=FALSE)
  expect_equal(dim(temp), c(80, 35))
  expect_equal(length(unique(c(temp))), 8)

  # test error with both m and n null
  expect_error(OSOAs_hadamard(optimize = TRUE))

  # test null m and el = 2
  temp <- OSOAs_hadamard(m = NULL, n = 48, el = 2, optimize = FALSE)
  expect_equal(dim(temp), c(48, 48 / 2 - 1))

  # test with null n and el = 2
  temp <- OSOAs_hadamard(m = 7, n = NULL, el = 2, optimize = FALSE)
  expect_equal(dim(temp), c(16, 16 / 2 - 1))
})
