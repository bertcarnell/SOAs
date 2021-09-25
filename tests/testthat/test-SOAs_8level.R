test_that("SOAs_8level", {
  suppressMessages(temp <- SOAs_8level(64))
  expect_equal(nrow(temp), 64)
  expect_true(all(apply(temp, 2, sum) == sum(temp[,1])))
  expect_s3_class(temp, "SOA")
  expect_equal(attr(temp, "type"), "OSOA")
  expect_equal(attr(temp, "strength"), "3+")
  expect_equal(dim(temp), c(64, 15))

  suppressMessages(temp <- SOAs_8level(16, m = 4))
  expect_equal(nrow(temp), 16)
  expect_true(all(apply(temp, 2, sum) == sum(temp[,1])))
  expect_s3_class(temp, "SOA")
  expect_equal(attr(temp, "type"), "SOA")
  expect_equal(attr(temp, "strength"), "3")
  expect_equal(dim(temp), c(16, 4))

  ## m=NULL
  suppressMessages(temp <- SOAs_8level(16, constr="ShiTang_alpha"))
  expect_equal(dim(temp), c(16, 5))
  # m=5
  suppressMessages(temp <- SOAs_8level(16, m = 5, constr="ShiTang_alpha"))
  expect_equal(nrow(temp), 16)
  expect_true(all(apply(temp, 2, sum) == sum(temp[,1])))
  expect_s3_class(temp, "SOA")
  expect_equal(attr(temp, "type"), "SOA")
  expect_equal(attr(temp, "strength"), "3")
  expect_equal(dim(temp), c(16, 5))


  suppressMessages(temp <- SOAs_8level(64, m = 2))
  expect_true(all(apply(temp, 2, sum) == sum(temp[,1])))

  suppressMessages(temp <- SOAs_8level(16, m = 4, noptim.rounds = 3))
  expect_true(all(apply(temp, 2, sum) == sum(temp[,1])))

  suppressMessages(temp <- SOAs_8level(64, m = NULL, optimize = FALSE))
  expect_true(all(apply(temp, 2, sum) == sum(temp[,1])))
  expect_snapshot_output(temp)

  suppressMessages(temp <- SOAs_8level(64, dmethod = "euclidean"))
  expect_true(all(apply(temp, 2, sum) == sum(temp[,1])))

  expect_error(SOAs_8level(64, dmethod = "junk"))
  expect_error(SOAs_8level(64, constr = "junk"))
})
