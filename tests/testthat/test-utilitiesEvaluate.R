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

test_that("count_npairs", {
  temp <- count_npairs(oa1)
  expect_true(is.list(temp))
  expect_equal(names(temp), c("paircounts", "columnpaircounts"))
})

test_that("count_nallpairs", {
  expect_equal(count_nallpairs(c(2,3,4)), c(2*3, 2*4, 3*4))
})

