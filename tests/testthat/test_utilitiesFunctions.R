context("Test utilities functions")

#create df for testing
test_df = data.frame(col1=factor(c("a", "a","a", "b","b", 1))
                     ,col2= factor(c(1, 3, 2, 5, 3, 3))
                    )

test_that(".changeFactor", {
  expect_equal(.changeFactor(test_df)$col1, c("a", "a","a", "b", "b", 1))
  expect_equal(.changeFactor(test_df)$col2, c("1", "3", "2", "5", "3", "3"))
  expect_equal(class(test_df$col2), "factor")
  expect_equal(class(.changeFactor(test_df)$col2), "character")
})

test_that(".noNA", {
  expect_equal(.noNA(NA), "")
  expect_equal(.noNA("hello"), "hello")
  expect_equal(.noNA(""), "")
  expect_equal(.noNA("FALSE"), "FALSE")
  expect_equal(.noNA(FALSE), FALSE)
})

test_that("%notin%", {
  expect_equal(test_df$col1 %notin% test_df$col2, c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE))
  expect_equal("" %notin% c("A", "B"), TRUE)
  expect_equal(NA %notin% c("A", "B"), TRUE)
  expect_equal(NULL %notin% c("A", "B"), logical(0))
})
  
test_that("%eq%", {
  expect_equal(as.character(test_df$col1) %eq% as.character(test_df$col2), c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
  expect_equal("A" %eq% c("A", "B"), c(TRUE, FALSE))
  expect_equal(c("A","B") %eq% c("A", "B"), c(TRUE, TRUE))
  expect_equal(c("a","B") %eq% c("A", "B"), c(FALSE, TRUE))
  expect_equal("" %eq% c("A", "B"), c(FALSE, FALSE))
  expect_equal(NA %eq% c("A", "B"), c(FALSE, FALSE))
  expect_warning(c("A","B", "A", 1, 2) %eq% c("A", "B"))
  expect_equal(c("A", "B")  %eq% c("B", NA), c(FALSE, FALSE))
  expect_equal(c("A", "B")  %eq% c("A", NA), c(TRUE, FALSE))
})

test_that("%noteq%", {
  expect_equal(c("A", "B")  %noteq% c("B", NA), c(TRUE, TRUE))
  expect_equal(c("A", "B")  %noteq% c("A", NA), c(FALSE, TRUE))
})
  
