context("miscellaneous")

test_that("fast_read works", {
  return()
  
  dat <- dummy_data()
  write.csv(dat, "data.csv", row.names = FALSE)
  dat2 <- fast_read("data.csv")
  expect_equal(dat, dat2)
  unlink("data.csv")
})
