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

test_that("ocheck", {
  expect_true(ocheck(oa1))
  expect_true(ocheck(st1))
  expect_true(ocheck(as.data.frame(oa1)))
  capture_output(expect_snapshot_output(ocheck(nullcase, verbose = TRUE)))
})

test_that("ocheck3", {
  temp <- OSOAs_LiuLiu(DoE.base::L81.3.10, optimize=FALSE)
  expect_false(ocheck3(st1))
  expect_true(ocheck3(temp))
  expect_true(ocheck3(as.data.frame(temp)))

  expect_true(ocheck3(oa1))
  capture_output(expect_snapshot_output(ocheck3(oa1, verbose = TRUE)))

  expect_false(ocheck3(oa2))
  capture_output(expect_snapshot_output(ocheck3(oa2, verbose = TRUE)))
})

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
  expect_error(Spattern(nullcase, s=4))
  capture_output(expect_snapshot(Spattern(nullcase, s=2)))

  ## s = number of levels: works but uses the GWLP
  capture_output(expect_snapshot(Spattern(nullcase, s=8)))
  capture_output(expect_snapshot(Spattern(nullcase, s=8, maxdim=3, maxwt=3)))

  ## s = number of levels: does not resort to GWLP
  ##     because maxdim (3) and maxwt (default 4) are unequal
  expect_error(Spattern(nullcase, s=8, maxdim=3))

  capture_output(expect_snapshot_output(Spattern(st1, s=2, maxdim=4, maxwt=4)))
  capture_output(expect_snapshot_output(Spattern(st1, s=2, maxdim=3, maxwt=2, detailed=TRUE)))
  capture_output(expect_snapshot_output(Spattern(nullcase4, s=2, maxwt=NULL)))
  capture_output(expect_snapshot_output(Spattern(nullcase4, s=2, maxdim=NULL)))
  capture_output(expect_snapshot_output(Spattern(nullcase4, s=2, maxdim=NULL, maxwt=NULL)))

})


test_that("count_npairs", {
  temp <- count_npairs(oa1)
  expect_true(is.list(temp))
  expect_equal(names(temp), c("paircounts", "columnpaircounts"))
})

test_that("count_nallpairs", {
  expect_equal(count_nallpairs(c(2,3,4)), c(2*3, 2*4, 3*4))
})

