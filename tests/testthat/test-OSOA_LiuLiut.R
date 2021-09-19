test_that("OSOA_LiuLiut", {
  temp <- OSOA_LiuLiut(lhs::createBose(4, 4))
  expect_equal(dim(temp), c(16, 4))

  temp <- OSOA_LiuLiut(DoE.base::L16.4.5, t = 2, m = 4)
  expect_equal(dim(temp), c(16, 4))
  expect_error(OSOA_LiuLiut(DoE.base::L16.4.5, t = 4, m = 4))

  # test the oa as a data.frame
  oa <- as.data.frame(lhs::createBose(4, 4))
  temp <- OSOA_LiuLiut(oa)
  expect_equal(dim(temp), c(16, 4))

  # test a oa data.frame with a factor
  oa[,2] <- as.factor(oa[,2])
  temp <- OSOA_LiuLiut(oa)
  expect_equal(dim(temp), c(16, 4))

  # test a null t gets a t==3
  oa <- lhs::createBusht(8, 5, 3)
  temp <- OSOA_LiuLiut(oa)
  expect_equal(dim(temp), c(512, 2))

  # test a null t with strength 4 oa gets t==4
  oa <- lhs::createBusht(5, 4, 4)
  temp <- OSOA_LiuLiut(oa)
  expect_equal(dim(temp), c(625, 2))
})
