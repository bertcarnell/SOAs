test_that("SOAs_8level", {
  suppressMessages(temp <- SOAs_8level(64))
  expect_equal(nrow(temp), 64)
  expect_true(all(apply(temp, 2, sum) == sum(temp[,1])))
  expect_equal(attr(temp, "type"), "SOA")

  suppressMessages(temp <- SOAs_8level(64, m = 2))
  expect_true(all(apply(temp, 2, sum) == sum(temp[,1])))

  suppressMessages(temp <- SOAs_8level(64, m = 2, constr = "ShiTang_alpha"))
  expect_true(all(apply(temp, 2, sum) == sum(temp[,1])))

  suppressMessages(temp <- SOAs_8level(64, m = 4, noptim.rounds = 3))
  expect_true(all(apply(temp, 2, sum) == sum(temp[,1])))

  suppressMessages(temp <- SOAs_8level(64, m = NULL, optimize = FALSE))
  expect_true(all(apply(temp, 2, sum) == sum(temp[,1])))

  suppressMessages(temp <- SOAs_8level(64, dmethod = "euclidean"))
  expect_true(all(apply(temp, 2, sum) == sum(temp[,1])))

  expect_error(SOAs_8level(64, dmethod = "junk"))
  expect_error(SOAs_8level(64, constr = "junk"))
})
