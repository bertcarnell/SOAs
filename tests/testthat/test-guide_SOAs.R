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


    ## el=2
    ## s=2
    ## n only given, m=22 obtained for HCT and ZT
    temp <- guide_SOAs(s=2, el=2, n=32)
    expect_snapshot_output(temp)
    ## m only given, n=32 for HCT, n=48 for ZT
    temp <- guide_SOAs(s=2, el=2, m=20)
    expect_snapshot_output(temp)
    ## n=16, m=1
    temp <- guide_SOAs(s=4, el=2, n=32)
    expect_snapshot(temp)
    ## n=256, m=45 (HCT) or m=21 (ZT)
    temp <- guide_SOAs(s=4, el=2, n=256)
    ## n=256, mmax=45 (HCT) or mmax=21 (ZT)
    temp <- guide_SOAs(s=4, el=2, m=15)
    expect_snapshot_output(temp)

    ## s=3
    ## m=6 (HCT) or m=4 (ZT)
    temp <- guide_SOAs(s=3, el=2, n=27)
    expect_snapshot_output(temp)
    ## n=81 for HCT, n=243 for ZT
    temp <- guide_SOAs(s=3, el=2, m=20)
    expect_snapshot_output(temp)
    ## n=81, m=25 (HCT) or m=13 (ZT)
    temp <- guide_SOAs(s=3, el=2, n=81)
    expect_snapshot_output(temp)
    ## n=729, m=18 (HCT) or m=10 (ZT)
    temp <- guide_SOAs(s=9, el=2, n=729)
    expect_snapshot_output(temp)
    ## n=512, mmax=16 (HCT) or n=512 mmax=9 (ZT)
    temp <- guide_SOAs(s=8, el=2, m=8)
    expect_snapshot_output(temp)

    ## el=3
    ## s=2
    ## m=22
    temp <- guide_SOAs(s=2, el=3, n=32)
    expect_snapshot_output(temp)
    ## n=32 for HCT, n=48 for ZT
    temp <- guide_SOAs(s=2, el=3, m=20)
    expect_snapshot_output(temp)

    ## n=256, m=20 LLY
    temp <- guide_SOAs(s=4, el=3, n=256)
    expect_snapshot_output(temp)
    ## n=256, mmax=20 LLY
    temp <- guide_SOAs(s=4, el=3, m=15)
    expect_snapshot_output(temp)
    ## n=4096, mmax=340 LLY
    temp <- guide_SOAs(s=4, el=3, m=100)
    expect_snapshot_output(temp)
    ## n=625, mmax=30 LLY
    temp <- guide_SOAs(s=5, el=3, m=10)
    expect_snapshot_output(temp)
})
