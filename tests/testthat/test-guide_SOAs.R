## guidance_testcases

test_that("guide_SOAs", {

    ## errors
    expect_error(guide_SOAs(s=2, el=2, m=20, n=16),
               regexp="m >= n is not permitted",
               fixed=TRUE)
    expect_error(guide_SOAs(s=3, el=3, n=9),
               regexp="more levels than runs",
               fixed=TRUE)
    expect_error(guide_SOAs(s=2, el=3, n=9),
                 regexp="n is not a multiple of s^el",
                 fixed=TRUE)
    expect_error(guide_SOAs(s=2.5, el=3, n=9),
                 regexp="s%%1 == 0 is not TRUE",
                 fixed=TRUE)
    expect_error(guide_SOAs(s=2, el=3.2, n=9),
                 regexp="el%%1 == 0 is not TRUE",
                 fixed=TRUE)
    expect_error(guide_SOAs(s=2, el=2, n=8.2),
                 regexp="n%%1 == 0 || is.null(n) is not TRUE",
                 fixed=TRUE)
    expect_error(guide_SOAs(s=2, el=2, m=8.2),
                 regexp="m%%1 == 0 || is.null(m) is not TRUE",
                 fixed=TRUE)
    expect_error(guide_SOAs(s=2, el=4),
                 regexp="el %in% 2:3 is not TRUE",
                 fixed=TRUE)
    expect_error(guide_SOAs(s=3, el=2, n=9),
                 regexp="n < s^3 is not permitted",
                 fixed=TRUE)


    ## el=2
    ## s=2
    temp1 <- guide_SOAs(s=2, el=2)
    expect_snapshot_output(temp1)
    temp2 <- guide_SOAs(s=2, el=2, n=32, m=21)
    expect_snapshot_output(temp2)

    ## n only given, m=22 obtained for HCT and ZT
    temp3 <- guide_SOAs(s=2, el=2, n=32)
    expect_snapshot_output(temp3)
    ## m only given, n=32 for HCT, n=48 for ZT
    temp4 <- guide_SOAs(s=2, el=2, m=20)
    expect_snapshot_output(temp4)
    ## 64 runs
    temp5 <- guide_SOAs(s=4, el=2, n=80)
    expect_snapshot(temp5)
    ## n=256, m=45 (HCT) or m=21 (ZT)
    temp6 <- guide_SOAs(s=4, el=2, n=256)
    expect_snapshot(temp6)
    temp7 <- guide_SOAs(s=4, el=2, m=15)
    expect_snapshot_output(temp7)

    ## s=3
    ## m=6 (HCT) or m=4 (ZT)
    temp8 <- guide_SOAs(s=3, el=2)
    expect_snapshot_output(temp8)
    temp9 <- guide_SOAs(s=3, el=2, n=27)
    expect_snapshot_output(temp9)
    ## n=81 for HCT, n=243 for ZT
    temp10 <- guide_SOAs(s=3, el=2, m=20)
    expect_snapshot_output(temp10)
    ## n=81, m=25 (HCT) or m=13 (ZT)
    temp11 <- guide_SOAs(s=3, el=2, n=81)
    expect_snapshot_output(temp11)
    ## n=729, m=18 (HCT) or m=10 (ZT)
    temp12 <- guide_SOAs(s=9, el=2, n=729)
    expect_snapshot_output(temp12)
    ## n=512, mmax=16 (HCT) or n=512 mmax=9 (ZT)
    temp13 <- guide_SOAs(s=8, el=2, m=8)
    expect_snapshot_output(temp13)

    ## el=3
    ## s=2
    temp14 <- guide_SOAs(s=2, el=3)
    expect_snapshot_output(temp14)
    temp15 <- guide_SOAs(s=2, el=3, n=32, m=12)
    expect_snapshot_output(temp15)
    ## this is an exception for Shi and Tang's family 1
    ##    m only 9 and not 10
    temp16 <- guide_SOAs(s=2, el=3, n=32)
    expect_snapshot_output(temp16)
    temp17 <- guide_SOAs(s=2, el=3, m=20)
    expect_snapshot_output(temp17)
    ## creates Shi and Tang family 2 (one too many for fam 3)
    temp18 <- guide_SOAs(s=2, el=3, m=32)
    expect_snapshot_output(temp18)

    ## n=256, m=20 LLY
    temp19 <- guide_SOAs(s=4, el=3, n=256)
    expect_snapshot_output(temp19)
    ## n=256, mmax=20 LLY
    temp20 <- guide_SOAs(s=4, el=3, m=15)
    expect_snapshot_output(temp20)
    ## n=4096, mmax=340 LLY
    temp21 <- guide_SOAs(s=4, el=3, m=100)
    expect_snapshot_output(temp21)
    ## n=625, mmax=30 LLY
    temp22 <- guide_SOAs(s=5, el=3, m=10)
    expect_snapshot_output(temp22)
    ## empty output without error
    temp23 <- guide_SOAs(s=2, el=2, n=8)
    expect_snapshot_output(temp23)
})
