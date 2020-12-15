context("Test Basic Logic Functions")

test_that("is.nan returns FALSE for NA", {
  expect_equal(is.nan(NA), FALSE)
})

test_that("is.nan returns TRUE for vector of NaNs", {
  expect_equal(all(is.nan(rep(NaN, 5))), TRUE)
})

test_that("is.prime returns TRUE for 13", {
  expect_equal(is.prime(13), TRUE)
})

test_that("is.prime returns FALSE for 25", {
  expect_equal(is.prime(25), FALSE)
})

test_that("is.prime returns TRUE for 0", {
  expect_equal(is.prime(0), FALSE)
})

test_that("is.within returns TRUE for 5 in c(0,10)", {
  expect_equal(is.within(5, c(0,10)), TRUE)
})

test_that("is.within returns TRUE for 'e' in c('a','z')", {
  expect_equal(is.within('e', c('a','z')), TRUE)
})

test_that("is.prime finds 25 primes between 1 and 100 (inclusive)", {
  expect_equal(sum(is.prime(1:100)), 25)
})

test_that("is.prime finds 168 primes between 1 and 1000 (inclusive)", {
  expect_equal(sum(is.prime(1:1000)), 168)
})

test_that("which.unique finds 10 unique entries in 1:10", {
  expect_equal(length(which.unique(1:10)), 10)
})

test_that("which.unique finds 1 unique entries in repeated list", {
  expect_equal(length(which.unique(rep('a', 5))), 1)
})

test_that("which.within finds 6 entries inside 1:100 between 5 and 10", {
  expect_equal(length(which.within(c(1:100), c(5,10))), 6)
})

