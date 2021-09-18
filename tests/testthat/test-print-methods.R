test_that("print.SOA", {
  temp <- OSOAs_regular(s=3, k=3, optimize=FALSE)
  expect_output(print(temp), regexp = "OSOA")
})

test_that("print.MDLE", {
  suppressMessages(temp <- MDLEs(DoE.base::L16.4.5, 2, noptim.rounds = 1))
  expect_output(print(temp), regexp = "MDLE")
})
