test_that("OSOAs_hadamard", {
  temp <- OSOAs_hadamard(m = 6)
  expect_equal(class(temp)[1], "SOA")
  expect_equal(dim(temp), c(16, 6))
  expect_equal(attr(temp, "type"), "OSOA")

  ## TODO:  verify this the number of columns that should be produced:  Issue #10
  ## 9 8-level factors in 24 runs
  temp <- OSOAs_hadamard(9, optimize=FALSE)
  expect_equal(dim(temp), c(24, 9))
  expect_equal(length(unique(c(temp))), 8)

  ## 6 8-level factors in 16 runs
  temp <- OSOAs_hadamard(n=16, optimize=FALSE)
  expect_equal(dim(temp), c(16, 6))
  expect_equal(length(unique(c(temp))), 8)

  temp <- OSOAs_hadamard(n=24, m=6, optimize=FALSE)
  expect_equal(dim(temp), c(24, 6))
  expect_equal(length(unique(c(temp))), 8)

  ## 35 8-level factors in 80 runs
  temp <- OSOAs_hadamard(m=35, optimize=FALSE)
  expect_equal(dim(temp), c(80, 35))
  expect_equal(length(unique(c(temp))), 8)
})
