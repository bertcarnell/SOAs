oa1 <- lhs::createBusht(3, 4, 3)
nullcase <- matrix(0:7, nrow=8, ncol=4)
## Shi and Tang strength 3+ construction in 7 8-level factors for 32 runs
st1 <- SOAs8level(32, optimize=FALSE)

test_that("ocheck.default", {
  expect_true(ocheck(oa1))
})

test_that("ocheck.SOA", {
  expect_false(ocheck(st1))
})

test_that("soacheck2D", {
  expect_false(soacheck2D(nullcase, s=2))
  expect_true(soacheck2D(st1, s=2, el=3, t=4))
})

test_that("soacheck3D", {
  expect_false(soacheck3D(nullcase, s=2))
  expect_true(soacheck3D(st1, s=2, el=3, t=4))
})
