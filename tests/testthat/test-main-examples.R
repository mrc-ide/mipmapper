context("Main")

test_that("`filter_coverage()` works", {
  
})

test_that("`filter_misc()` works", {
  
})

test_that("`melt_mip_data()` works", {
  
  data <- data.frame(
    Sample_ID = c(rep("a", 3), rep("b", 2)), 
    Chrom = c("chr1","chr1", "chr2", "chr1", "chr1"), 
    Pos = c(100, 200, 50, 100, 200), 
    Coverage = c(47, 95, 100, 52, 100), 
    Barcode_Count = c(47, 0, 40, 52, 70)
  )
  
  expect_equal(
    melt_mip_data(data),
    data.frame("Sample_ID" = c("a","b"),
               "chr1_100" = c(1, 1),
               "chr1_200" = c(0.0, 0.7),
               "chr2_50" = c(0.4, NA))
  )
  
})
