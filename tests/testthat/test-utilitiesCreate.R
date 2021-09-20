## s=3, k=3
temp <- createAB(3)
expect_snapshot_output(temp)

## s=2 is special case
temp <- createAB(2, k=4)
expect_snapshot_output(temp)

temp <- BcolsFromBcolllist(list(1:2, 1:3, 4:6, 1))
expect_equal(temp[1], 2)
expect_equal(temp[2], 3)
expect_equal(temp[3], 4)
expect_equal(temp[4], 1)
