context("plot")

test_that("`plot_coverage()` works", {
  dat <- dummy_data()
  dat2 <- filter_misc(dat)
  expect_identical(plot_coverage(dat2), dat2)
})

test_that("`plot_pca_variance()` works", {
   dat <- dummy_data()
   dat <- filter_misc(dat = dat)
   dat <- filter_coverage(dat = dat, min_coverage = 2)
   dat <- melt_mip_data(dat = dat)
   dat <- impute_mip_data(dat = dat)
   pca <- pca_mip_data(dat = dat)
   out <- plot_pca_variance(pca, num_components = 3)
   expect_true(inherits(out, "plotly"))
})

test_that("`plot_pca()` works", {
   dat <- dummy_data()
   dat <- filter_misc(dat = dat)
   dat <- filter_coverage(dat = dat, min_coverage = 2)
   dat <- melt_mip_data(dat = dat)
   dat <- impute_mip_data(dat = dat)
   pca <- pca_mip_data(dat = dat)
   plot_pca(pca, num_components = 2, meta_var = "Country")
   plot_pca(pca, num_components = 3, meta_var = "Study")
   expect_error(plot_pca(pca, num_components = 3, meta_var = "gibberish"))
   expect_message(plot_pca(pca, num_components = 4, meta_var = "Country"))
})
