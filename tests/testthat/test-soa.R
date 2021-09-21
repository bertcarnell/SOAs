test_that("soa", {
  ## t=5, random
  oa <- (desnum(FrF2::FrF2(128, 9, randomize=FALSE))+1)/2
  temp <- soa(oa, t=5)
  expect_equal(dim(temp), c(128, 4))
  expect_equal(range(temp), c(0, 31))
  expect_error(soa(oa, t=5, m=5), regexp="morig <= m is not TRUE", fixed=TRUE)
  temp <- soa(oa, t=5, permlist=list(list(0:1, 0:1, 0:1, 0:1, 0:1),
                                       list(0:1, 0:1, 0:1, 0:1, 0:1),
                                       list(0:1, 0:1, 0:1, 0:1, 0:1),
                                       list(0:1, 0:1, 0:1, 0:1, 0:1)))
  expect_snapshot_output(temp)
})
