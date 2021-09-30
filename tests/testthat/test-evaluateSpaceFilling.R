
test_that("mindist", {
  D1 <- matrix(1:7,7,5)
  dfD1 <- as.data.frame(D1)
  expect_equal(mindist(D1), 5)
  expect_equal(mindist(dfD1), 5)
  expect_equal(mindist(D1, "euclidean"), sqrt(5))
  expect_error(mindist(D1, "jaccard"))
})
