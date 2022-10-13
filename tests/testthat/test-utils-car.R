test_that("nu.d function",{

  p <- 0.3

  expect_equal(nu.d("simple", p_trt=1/2), 1/4)
  expect_equal(nu.d("simple", p_trt=p), p*(1-p))
  expect_equal(nu.d("permuted-block"), 0)
  expect_equal(nu.d("biased-coin"), 0)
  expect_equal(nu.d("urn"), NA)
  expect_equal(nu.d("whatever"), NA)
})
