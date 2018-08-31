context("analysis")

test_that("`pca_mip_data()` works", {
  return()
  
  dat <- dummy_data()
  dat <- filter_misc(dat = dat)
  dat <- filter_coverage(dat = dat, min_coverage = 2)
  dat <- melt_mip_data(dat)
  dat <- impute_mip_data(dat = dat)
  pca <- pca_mip_data(dat)

  expect_equal(pca$var[1], 100)
  expect_equal(pca$loadings[1, 1], 1)
})
