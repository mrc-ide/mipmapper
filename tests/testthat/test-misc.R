context("miscellaneous")

test_that("`dummy_data()` works", {

  dat <- dummy_data()
  expect_equal(dim(dat), c(11, 13))

  dat2 <- filter_misc(dat)
  expect_equal(dim(dat2), c(5, 13))
  expect_identical(dat2$Alt[1], "G,T")

  dat2 <- filter_misc(dat, SNP_only = FALSE)
  expect_equal(dim(dat2), c(7, 13))

  dat2 <- filter_misc(dat, group_Alt = FALSE)
  expect_equal(dim(dat2), c(6, 13))
  expect_identical(dat2$Alt[1], "G")

  dat2 <- filter_misc(dat, drop_irregular = FALSE)
  expect_equal(dim(dat2), c(8, 13))

  })

test_that("fast_read works", {
  
  dat <- dummy_data()
  write.csv(dat, "data.csv", row.names = FALSE)
  dat2 <- fast_read("data.csv")
  expect_equal(dat, dat2)
  unlink("data.csv")
})
