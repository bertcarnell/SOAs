## s=3, k=3
temp <- createAB(3)
expect_snapshot_output(temp)

## s=2 is special case
temp <- createAB(2, k=4)
expect_snapshot_output(temp)
