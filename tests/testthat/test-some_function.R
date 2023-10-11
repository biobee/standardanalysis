# Name can be very useful in identifying what test fails
# tests will always use most recent version of functions so no need for rebuilding
test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("addition works", {
  expect_equal(2 + 6, 8)
  actual_addition <- 7 + 4
  expected_result <- 11
  expect_equal(actual_addition, expected_result)
  expect_true(actual_addition == expected_result)
  expect_false(actual_addition != expected_result)
  # expect_error()
})

# Tip from Barb, start with failing tests
#test_that("dimensions of matrix are ok", {
#  expect_equal(dim(make_groups(classmates)), c(3,2))
#})
