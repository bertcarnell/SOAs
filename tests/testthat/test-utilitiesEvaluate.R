oa1 <- lhs::createBusht(3, 4, 3)
nullcase <- matrix(0:7, nrow=8, ncol=4)
## Shi and Tang strength 3+ construction in 7 8-level factors for 32 runs
st1 <- SOAs_8level(32, optimize=FALSE)

test_that("ocheck", {
  expect_true(ocheck(oa1))
  expect_true(ocheck(st1))
  expect_true(ocheck(as.data.frame(oa1)))
  expect_output(ocheck(nullcase, verbose = TRUE))
})

test_that("ocheck3", {
  temp <- OSOAs_LiuLiu(DoE.base::L81.3.10, optimize=FALSE)
  expect_false(ocheck3(st1))
  expect_true(ocheck3(temp))
  expect_true(ocheck3(as.data.frame(temp)))

  expect_true(ocheck3(oa1))
  expect_output(ocheck3(oa1, verbose = TRUE))
})

test_that("soacheck2D", {
  expect_false(soacheck2D(nullcase, s=2))
  expect_true(soacheck2D(st1, s=2, el=3, t=4))

  # test when el = 2 and t = 4
  expect_message(expect_true(soacheck2D(lhs::createBose(4, ncol=3), s=2, el=2, t=4)))
  # test when min(OA) == 1
  expect_true(soacheck2D(lhs::createBose(8, ncol=3) + 1, s=2, el=3, t=4))
  # test verbose output
  expect_output(expect_true(soacheck2D(st1, s=2, el=3, t=4, verbose=TRUE)))
})

test_that("soacheck3D", {
  expect_false(soacheck3D(nullcase, s=2))
  expect_true(soacheck3D(st1, s=2, el=3, t=4))

  # the result of this test depends on the random outcome of createBose;
  # apparently, soacheck3D is sometimes TRUE
  set.seed(12345)
  expect_false(soacheck3D(lhs::createBose(8, ncol=3) + 1, s=2, el=3, t=4))
  # test verbose output
  expect_output(expect_true(soacheck3D(st1, s=2, el=3, t=4, verbose=TRUE)))
  # test error when min(OA) == 2 (not permitted)
  expect_error(soacheck2D(lhs::createBose(8, ncol=3) + 2, s=2, el=3, t=4))
  expect_error(soacheck3D(lhs::createBose(8, ncol=3) + 2, s=2, el=3, t=4))
})

test_that("count_npairs", {
  temp <- count_npairs(oa1)
  expect_true(is.list(temp))
  expect_equal(names(temp), c("paircounts", "columnpaircounts"))
})

test_that("count_nallpairs", {
  expect_equal(count_nallpairs(c(2,3,4)), c(2*3, 2*4, 3*4))
})
