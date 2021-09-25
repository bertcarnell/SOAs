gf <- lhs::create_galois_field(9)

expect_all_equal <- function(x, y, ...)
{
  temp <- mapply(function(x, y) {expect_equal(x, y, ...)}, x, y)
  invisible(temp)
}

test_that("int2poly", {
  expect_all_equal(int2poly(0, gf), c(0,0))
  expect_all_equal(int2poly(0, gf), c(0,0))
  expect_all_equal(int2poly(1, gf), c(1,0))
  expect_all_equal(int2poly(2, gf), c(2,0))
  expect_all_equal(int2poly(3, gf), c(0,1))
  expect_all_equal(int2poly(4, gf), c(1,1))
  expect_all_equal(int2poly(5, gf), c(2,1))
  expect_all_equal(int2poly(6, gf), c(0,2))
  expect_all_equal(int2poly(7, gf), c(1,2))
  expect_all_equal(int2poly(8, gf), c(2,2))

  expect_equal(int2poly(c(1,4,7), gf), matrix(c(1,0,1,1,1,2), nrow = 3, byrow=TRUE))

  expect_error(int2poly(9, gf))
})

test_that("poly2int", {
  expect_equal(c(poly2int(c(0,0), gf)), 0)
  expect_equal(c(poly2int(c(1,0), gf)), 1)
  expect_equal(c(poly2int(c(2,0), gf)), 2)
  expect_equal(c(poly2int(c(0,1), gf)), 3)
  expect_equal(c(poly2int(c(0,2), gf)), 6)

  expect_all_equal(c(poly2int(as.matrix(expand.grid(0:2, 0:2)), gf)), 0:8)
})

test_that("gf_sum", {
  expect_equal(gf_sum(1, 2, gf), 0)
  expect_equal(gf_sum(6, 8, gf), 5)
  expect_all_equal(gf_sum(c(3,5), c(0,1), gf), c(3, 3))
})

test_that("gf_prod", {
  expect_equal(gf_prod(0, 2, gf), 0)
  expect_equal(gf_prod(1, 2, gf), 2)
  expect_equal(gf_prod(3, 2, gf), 6)
  expect_equal(gf_prod(1, 3, gf), 3)
  expect_equal(gf_prod(8, 8, gf), 5)
  expect_all_equal(gf_prod(c(0,1,3), c(2,2,2), gf), c(0, 2, 6))
})

test_that("gf_sum_list", {
  expect_all_equal(gf_sum_list(list(a = c(1,0), b=c(2,5), d=c(0,8)), gf, checks = TRUE), c(0, 1))
  expect_error(gf_sum_list(list(a = rep(0:9, each = 8),
                           b=rep(0:9, times = 8)), gf, checks = TRUE),
                           regexp="invalid numbers occur in ll",
               fixed=TRUE)
  expect_error(gf_sum_list(list(a = rep(0:8, each = 8),
                                b=rep(0:8, times = 9)), gf, checks = TRUE),
               regexp="all elements of ll must have the same length",
               fixed=TRUE)

  temp <- gf_sum_list(list(a = rep(0:8, each = 81),
                           b=rep(0:8, each = 9, times = 9),
                           d=rep(0:8, times = 81)), gf, checks = TRUE)
  counter <- 1
  for (i in 0:8)
  {
    for (j in 0:8)
    {
      for (k in 0:8)
      {
        expect_true(gf_sum(gf_sum(i, j, gf), k, gf) == temp[counter])
        counter <- counter + 1
      }
    }
  }
})

test_that("gf_matmult", {
  M1 <- matrix(0:8, ncol = 1)
  M2 <- matrix(0:8, nrow = 1)
  expect_all_equal(c(gf_matmult(M1, M2, gf, checks = TRUE)), c(gf$times))

  M1 <- matrix(0:3, ncol = 2, nrow = 2)
  M2 <- matrix(0:3, ncol = 2, nrow = 2)
  expect_all_equal(c(gf_matmult(M1, M2, gf, checks = TRUE)), c(2, 3, 6, 6))
})

