test_that("utilitiesCreate", {
  ## s=3, k=3
temp1 <- createAB(3)
## list of three 27x6 matrices
expect_snapshot_output(temp1)

## s=2 is special case
temp2 <- createAB(2, k=4)
## list of three 16x10 matrices
expect_snapshot_output(temp2)

temp <- BcolsFromBcolllist(list(1:2, 1:3, 4:6, 1))
expect_equal(temp[1], 2)
expect_equal(temp[2], 3)
expect_equal(temp[3], 4)
expect_equal(temp[4], 1)
})
