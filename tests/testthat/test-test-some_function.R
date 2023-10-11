test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})




# Write multiple tests for a single function, to test that the function deals well with multiple kinds of input.
# Remember you can put multiple assertions in a single test (an assertion is stating an expectation, e.g. expect_equal, expect_true).
# Check if your function produces specific error messages (using e.g. expect_error).

# Write a test for (part of) your workflow, in which multiple functions are used (remember that this is called integration testing).

# Use covr::report() to identify parts of your package that could benefit from more tests.
