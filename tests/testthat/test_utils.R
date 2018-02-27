library(TWAS)
context("utils")
test_that("vec2str creates a truncated string of its arguments", {
  expect_equal(vec2str('foo'),'foo')
})
