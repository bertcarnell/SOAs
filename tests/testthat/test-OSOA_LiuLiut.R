test_that("OSOA_LiuLiut", {
  temp <- OSOA_LiuLiut(lhs::createBose(4, 4))
  expect_equal(dim(temp), c(16, 4))

  temp <- OSOA_LiuLiut(DoE.base::L16.4.5, t = 2, m = 4)
  expect_equal(dim(temp), c(16, 4))
  expect_error(OSOA_LiuLiut(DoE.base::L16.4.5, t = 4, m = 4))
})
