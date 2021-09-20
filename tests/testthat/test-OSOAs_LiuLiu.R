test_that("OSOAs_LiuLiu", {
  temp <- OSOAs_LiuLiu(lhs::createBose(4, 4), optimize=FALSE)
  expect_equal(dim(temp), c(16, 4))

  suppressMessages(temp <- OSOAs_LiuLiu(DoE.base::L16.4.5, t = 2, m = 3))
  expect_equal(dim(temp), c(16, 3))
  expect_s3_class(temp, "SOA")
  expect_equal(attr(temp, "strength"), 2)
  expect_equal(attr(temp, "type"), "OSOA")
  expect_error(OSOAs_LiuLiu(DoE.base::L16.4.5, t = 4, m = 4))

  temp <- OSOAs_LiuLiu(FrF2::FrF2(32,6), t = 3, m = 2, optimize=FALSE)
  expect_equal(attr(temp, "strength"), 3)
  expect_equal(attr(temp, "type"), "OSOA")
  expect_equal(dim(temp), c(32, 2))

  temp <- OSOAs_LiuLiu(FrF2::FrF2(32, 6, randomize=FALSE), t = 4, m = 2, optimize=FALSE)
  expect_snapshot_output(temp)
  expect_equal(attr(temp, "strength"), 4)
  expect_equal(attr(temp, "type"), "OSOA")
  expect_equal(dim(temp), c(32, 2))

  expect_error(suppressMessages(OSOAs_LiuLiu(FrF2::FrF2(32,6), t = 3, m = 3), regexp="m <= boundm", fixed=TRUE))
  expect_error(OSOA_LiuLiut(DoE.base::L16.4.5, t = 4, m = 4))

  ## check 3-level
  suppressMessages(temp <- OSOAs_LiuLiu(lhs::createBusht(3, 4, 4)))
  expect_equal(attr(temp, "strength"), 4)
  expect_equal(attr(temp, "type"), "OSOA")
  expect_equal(dim(temp), c(81, 2))
})
