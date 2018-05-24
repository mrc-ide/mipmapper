context("plot")

test_that("`plot_coverage()` works", {
  dat <- dummy_data()
  dat2 <- filter_misc(dat)
  plot_coverage(dat2)
})
