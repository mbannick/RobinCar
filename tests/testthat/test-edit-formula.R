test_that("edit formula", {
  form <- "y ~ A:(x1+x2)"
  edited <- .edit.formula(form, response_col="y", treat_col="A")
  expect_equal(edited[[1]], "response ~ treat:(x1 + x2)")
  expect_equal(edited[[2]], c("x1", "x2"))
})

test_that("edit formula", {
  form <- "y ~ A*(x1+x2)"
  edited <- .edit.formula(form, response_col="y", treat_col="A")
  expect_equal(edited[[1]], "response ~ treat * (x1 + x2)")
  expect_equal(edited[[2]], c("x1", "x2"))
})

test_that("edit formula with odd covariates 1", {
  form <- "y ~ A + y1"
  edited <- .edit.formula(form, response_col="y", treat_col="A")
  expect_equal(edited[[1]], "response ~ treat + y1")
  expect_equal(edited[[2]], c("y1"))
})

test_that("edit formula with odd covariates 2", {
  form <- "y ~ A + A1"
  edited <- .edit.formula(form, response_col="y", treat_col="A")
  expect_equal(edited[[1]], "response ~ treat + A1")
  expect_equal(edited[[2]], c("A1"))
})

