set.seed(123)
oa1 <- lhs::createBusht(3, 4, 3)
oa2 <- DoE.base:::L9.3.4-1
set.seed(12345)
oa3 <- lhs::createBose(8, ncol=3)
oa4 <- lhs::createBose(4, ncol=3)

nullcase <- matrix(0:7, nrow=8, ncol=4)
nullcase4 <- matrix(0:3, nrow=4, ncol=4)
## Shi and Tang strength 3+ construction in 7 8-level factors for 32 runs
st1 <- SOAs_8level(32, optimize=FALSE)

test_that("soacheck2D", {
  expect_false(soacheck2D(nullcase, s=2))
  capture_output(expect_false(soacheck2D(nullcase, s=2, verbose = TRUE)))
  expect_false(soacheck2D(nullcase, s=2, t=2))
  capture_output(expect_snapshot_output(soacheck2D(nullcase, s=2, t=2, verbose = TRUE)))
  expect_true(soacheck2D(st1, s=2, el=3, t=4))
  expect_true(soacheck2D(st1, s=2, el=3, t=3))

  # test when el = 2 and t = 4
  expect_message(expect_true(soacheck2D(oa4, s=2, el=2, t=4)))
  # test when el = 2 and t = 3
  expect_true(soacheck2D(oa4, s=2, el=2, t=3))
  # test when min(OA) == 1
  expect_true(soacheck2D(oa3 + 1, s=2, el=3, t=4))
  #
  expect_false(soacheck2D(nullcase, s=2, el=3, t=4))
  # test verbose output
  capture_output(expect_snapshot_output(expect_true(soacheck2D(st1, s=2, el=3, t=4, verbose=TRUE))))
  capture_output(expect_snapshot_output(expect_false(soacheck2D(nullcase, s=2, el=3, t=4, verbose=TRUE))))
  capture_output(expect_false(soacheck2D(nullcase, s=2, el=3, t=3, verbose=TRUE)))
  capture_output(expect_false(soacheck2D(nullcase4, s=2, el=2, t=3, verbose=TRUE)))

})

test_that("soacheck3D", {
  expect_false(soacheck3D(nullcase, s=2))
  capture_output(expect_false(soacheck3D(nullcase, s=2, verbose = TRUE)))
  expect_false(soacheck3D(nullcase, s=2, t=4))
  capture_output(expect_false(soacheck3D(nullcase, s=2, t=4, verbose = TRUE)))
  expect_true(soacheck3D(st1, s=2, el=3, t=4))
  expect_true(soacheck3D(st1, s=2, el=3, t=3))

  ## t=4, FALSE
  capture_output(expect_snapshot_output(expect_false(soacheck3D(oa3 + 1,
                                s=2, el=3, t=4, verbose = TRUE))))
  expect_false(soacheck3D(oa3 + 1, s=2, el=3, t=4))
  # test verbose output
  expect_true(soacheck3D(st1, s=2, el=3, t=4))
  # test error when min(OA) == 2 (not permitted)
  expect_error(soacheck2D(oa3 + 2, s=2, el=3, t=4))
  expect_error(soacheck3D(oa3 + 2, s=2, el=3, t=4))
  # test error when t is wrong
  expect_error(soacheck2D(oa3, s=2, el=3, t=5))

})


test_that("Spattern", {
  ## n not a power of s
  expect_error(Spattern(nullcase, s=4))

  ## mixed level D
  expect_error(Spattern(cbind(nullcase, rbind(nullcase4, nullcase4)), s=2))

  ## calculations correct ?
  temp <- Spattern(nullcase, s=2, maxdim=NULL, maxwt=NULL)
  attributes(temp) <- NULL
  expect_equal(temp, c(0, 6, 0, 13, 24, 36, 48, 128, 96, 96, 0, 64))

  ## with maxdim and maxwt specified larger than possible
  ##    is silently corrected
  temp <- Spattern(nullcase, s=2, maxdim=5, maxwt=20)
  attributes(temp) <- NULL
  expect_equal(temp, c(0, 6, 0, 13, 24, 36, 48, 128, 96, 96, 0, 64))

  ## same for data frame D
  temp <- Spattern(as.data.frame(nullcase), s=2, maxdim=NULL, maxwt=NULL)
  attributes(temp) <- NULL
  expect_equal(temp, c(0, 6, 0, 13, 24, 36, 48, 128, 96, 96, 0, 64))

  ## same for array starting with 1
  temp <- Spattern(nullcase+1, s=2, maxdim=NULL, maxwt=NULL)
  attributes(temp) <- NULL
  expect_equal(temp, c(0, 6, 0, 13, 24, 36, 48, 128, 96, 96, 0, 64))

  temp <- Spattern(nullcase, s=2, maxwt=NULL)
  ## unchanged, because the default maxdim 4 equals m
  attributes(temp) <- NULL
  expect_equal(temp, c(0, 6, 0, 13, 24, 36, 48, 128, 96, 96, 0, 64))

  temp <- Spattern(nullcase, s=2, maxdim=3, maxwt=NULL)
  ## reduced, missing the last three and earlier ones are reduced
  attributes(temp) <- NULL
  expect_equal(temp, c(0, 6, 0, 12, 24, 24, 48, 96, 0))

  temp <- Spattern(nullcase, s=2)
  attributes(temp) <- NULL
  expect_equal(temp, c(0, 6, 0, 13))

  temp <- Spattern(nullcase, s=2, maxdim=3)
  attributes(temp) <- NULL
  expect_equal(temp, c(0, 6, 0, 12))

  D1 <- cbind(c(18, 19, 11, 16, 20, 14,  4, 12, 10, 22,  2, 15,  1,  3, 23,  0,  8,  7,  9,  6, 24, 21, 13, 17,  5),
              c(16,  9,  1, 20, 22,  7, 17, 12, 15, 14, 21, 11,  5, 10,  3, 13, 23,  8,  4, 19, 18,  6, 24,  0,  2),
              c(14,  0,  2,  3, 12, 10,  1,  7, 24,  5, 22, 21,  4, 11,  8, 15,  6, 18, 13,  9, 19, 23, 17, 16, 20))

  ## s=5
  temp <- Spattern(D1, s=5, maxwt=NULL)
  attributes(temp) <- NULL
  expect_equal(temp[1:3], c(0, 0.64, 26.08))  ## Example 4 Tian and Xu
  expect_equal(sum(temp), (5^(3*2)/25)-1)     ## Theorem 4 Tian and Xu

  ## s = number of levels: works but uses the GWLP
  capture_output(expect_snapshot(Spattern(nullcase, s=8)))
  capture_output(expect_snapshot(Spattern(nullcase, s=8, maxdim=3, maxwt=3)))

  capture_output(expect_snapshot_output(Spattern(st1, s=2, maxdim=4, maxwt=4)))
  capture_output(expect_snapshot_output(dim_wt_tab(Spattern(st1, s=2))))
  capture_output(expect_snapshot_output(Spattern(nullcase4, s=2, maxwt=NULL)))
  capture_output(expect_snapshot_output(Spattern(nullcase4, s=2, maxdim=NULL)))
  capture_output(expect_snapshot_output(Spattern(nullcase4, s=2, maxdim=NULL, maxwt=NULL)))

})
