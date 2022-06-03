## guidance_testcases

test_that("guide_SOAs_from_OA", {

    ## el=2
    ## errors
    expect_error(guide_SOAs_from_OA(s=2, tOA=2, el=2, mOA=20, nOA=16),
               regexp="mOA >= nOA is not possible",
               fixed=TRUE)
    expect_error(guide_SOAs_from_OA(s=3, el=3, nOA=9, mOA=8, tOA=2),
               regexp="more levels than runs",
               fixed=TRUE)
    expect_error(guide_SOAs_from_OA(s=3, el=2, nOA=18, mOA=8, tOA=2),
                 regexp="There must be something wrong, an OA like this cannot exist.",
                 fixed=TRUE)
    expect_error(guide_SOAs_from_OA(s=2.2, tOA=2, el=2, mOA=8, nOA=16),
                 regexp="s%%1 == 0 is not TRUE",
                 fixed=TRUE)
    expect_error(guide_SOAs_from_OA(s=2, tOA=2.2, el=2, mOA=8, nOA=16),
                 regexp="tOA%%1 == 0 && tOA >= 2 is not TRUE",
                 fixed=TRUE)
    expect_error(guide_SOAs_from_OA(s=2, tOA=2, el=2.2, mOA=8, nOA=16),
                 regexp="el%%1 == 0 is not TRUE",
                 fixed=TRUE)
    expect_error(guide_SOAs_from_OA(s=2, tOA=3, el=4, mOA=8, nOA=16),
                 regexp="The chosen combination of tOA and el is not possible.",
                 fixed=TRUE)
    expect_error(guide_SOAs_from_OA(s=2, tOA=2, el=2, mOA=8.3, nOA=16),
                 regexp="mOA%%1 == 0 is not TRUE",
                 fixed=TRUE)
    expect_error(guide_SOAs_from_OA(s=2, tOA=2, el=2, mOA=8, nOA=16.5),
                 regexp="nOA%%1 == 0 is not TRUE",
                 fixed=TRUE)
    expect_error(guide_SOAs_from_OA(s=2, tOA=6, el=6, mOA=20, nOA=1024),
                 regexp="el %in% 2:5 is not TRUE",
                 fixed=TRUE)
    expect_error(guide_SOAs_from_OA(s=2, el=4, tOA=2, mOA=4, nOA=8),
                 regexp="OA has more levels than runs - not possible",
                 fixed=TRUE)

  ## three variants: two strength 2 with n=16,
  ## one strength 2+ or 3- with n=32
  ## HT or LL or ZT
  temp1 <- guide_SOAs_from_OA(2, 16, 15, 2, el=2)
  expect_snapshot_output(temp1)
  ## for el=3: LLY in 32 runs with m=14
  temp2 <- guide_SOAs_from_OA(2, 16, 15, 2, el=3)
  expect_snapshot_output(temp2)
  ## for el=3: tOA=3
  temp3 <- guide_SOAs_from_OA(2, 16, 8, 3, el=3)
  expect_snapshot_output(temp3)
  ## for el=4: m=2 81-level columns in 81 runs
  ## both with HT and with LL
  temp4 <- guide_SOAs_from_OA(3, 81, 5, 4, el=4)
  expect_snapshot_output(temp4)
  ## for el=2 in spite of tOA=4:
  ## m=5 9-level columns in 81 runs HT
  ## or m=4 orthogonal 9-level columns LL
  ## or n=243 runs with 5 9-level columns ZT
  temp5 <- guide_SOAs_from_OA(3, 81, 5, 4, el=2)
  expect_snapshot_output(temp5)
  temp6 <- guide_SOAs_from_OA(s=2, tOA=6, el=5, mOA=20, nOA=1024)
  expect_snapshot_output(temp6)
})
