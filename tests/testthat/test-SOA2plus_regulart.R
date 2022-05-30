test_that("SOA2plus_regulart", {
  set.seed(1234)
  temp <- SOA2plus_regulart(2, 4)
  expect_equal(dim(temp), c(16, 10))
  expect_false(ocheck(temp))
  set.seed(1234)
  temp <- SOA2plus_regulart(2, 4, m=5)
  expect_equal(dim(temp), c(16, 5))
  expect_true(ocheck(temp))
  set.seed(1234)
  temp <- SOA2plus_regulart(3, 3)
  expect_equal(dim(temp), c(27, 6))
  expect_true(ocheck(temp))
  set.seed(1234)
  temp <- SOA2plus_regulart(4, 3)
  expect_equal(dim(temp), c(64, 8))
  expect_true(ocheck(temp))

  temp <- SOA2plus_regulart(4, 3, m=3,
                            permlist = list(
                              list(0:3,0:3,0:3),
                              list(0:3,0:3,0:3),
                              list(0:3,0:3,0:3)))
  expect_equal(dim(temp), c(64,3))
  temp <- SOA2plus_regulart(4, 3,
                            permlist = list(
                              list(0:3,0:3,0:3),
                              list(0:3,0:3,0:3),
                              list(0:3,0:3,0:3),
                              list(0:3,0:3,0:3),
                              list(0:3,0:3,0:3),
                              list(0:3,0:3,0:3),
                              list(0:3,0:3,0:3),
                              list(0:3,0:3,0:3)))
  expect_equal(dim(temp), c(64,8))

  expect_error(SOA2plus_regulart(64, 3),
               regexp="must not be larger than s=2^5", fixed=TRUE)
  expect_error(SOA2plus_regulart(729, 3),
               regexp="must not be larger than s=3^3", fixed=TRUE)
})
