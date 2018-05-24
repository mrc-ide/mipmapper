context("Main")

context("`filter_coverage()`")

test_that("`filter_coverage()` works", {
  dat <- dummy_data()
  dat <- filter_misc(dat)
  dat2 <- filter_coverage(dat, min_coverage = 5)
  expect_true(dim(dat2)[1] < dim(dat)[1])
})


context("`filter_misc()`")

test_that("`filter_misc()` works", {
  dat <- dummy_data()
  dat2 <- filter_misc(dat)
  expect_true(dim(dat2)[1] < dim(dat)[1])
})

test_that("`melt_mip_data()` works", {

  dat <- data.frame(
    Sample_ID = c(rep("a", 3), rep("b", 2)),
    Chrom = c("chr1", "chr1", "chr2", "chr1", "chr1"),
    Pos = c(100, 200, 50, 100, 200),
    Coverage = c(47, 95, 100, 52, 100),
    Barcode_Count = c(47, 0, 40, 52, 70)
  )

  expect_equal(
    melt_mip_data(dat = dat),
    data.frame("Sample_ID" = c("a", "b"),
               "chr1_100" = c(1, 1),
               "chr1_200" = c(0.0, 0.7),
               "chr2_50" = c(0.4, NA))
  )

})


test_that("`impute_mip_data()` works", {

  dat <- data.frame(
    Sample_ID = c(rep("a", 3), rep("b", 2)),
    Chrom = c("chr1", "chr1", "chr2", "chr1", "chr1"),
    Pos = c(100, 200, 50, 100, 200),
    Coverage = c(47, 95, 100, 52, 100),
    Barcode_Count = c(47, 0, 40, 52, 70)
  )

  dat <- melt_mip_data(dat = dat)
  dat <- impute_mip_data(dat = dat)
  expect_equal(
    dat,
    structure(c(1, 2, 1, 1, 0, 0.7, 0.4, 0.4),
              .Dim = c(2L, 4L),
              .Dimnames = list(NULL,
                               c("", "chr1_100", "chr1_200", "chr2_50")
              )
    )
  )

})
