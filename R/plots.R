
#------------------------------------------------
# generate any number of divergent colours
#' @noRd
divergent_colours <- function(K) {
  
  # generate palette and colours
  #raw_cols <- c("#D73027", "#FC8D59", "#FEE090", "#E0F3F8", "#91BFDB", "#4575B4")
  raw_cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
  my_palette <- colorRampPalette(raw_cols)
  
  # simple case if small K
  if (K<=2) {
    return(my_palette(K))
  }
  
  # some logic to choose a palette size and sequence of colours that is
  # consistent across different values of K
  ncol <- 3
  while(ncol<K) {
    ncol <- ncol+(ncol-1)
  }
  dist_mat <- matrix(1:ncol, ncol, ncol)
  dist_mat <- abs(t(dist_mat)-dist_mat)
  x <- rep(FALSE, ncol)
  
  col_index <- 1
  for (i in 2:K) {
    x[col_index] <- TRUE
    s <- apply(dist_mat[which(x),,drop=FALSE], 2, min)
    next_index <- which.max(s)
    col_index <- c(col_index, next_index)
  }
  col_index
  ret <- my_palette(ncol)[col_index]
  
  return(ret)
}

#------------------------------------------------
#' @title Plot raster image of coverage matrix
#'
#' @description Plot raster image of coverage matrix.
#'
#' @param project a \code{mipmapper_project} with data already loaded
#'
#' @export

plot_coverage_matrix <- function(project) {
  
  # check inputs
  assert_custom_class(project, "mipmapper_project")
  if (is.null(project$data_processed$data_coverage)) {
    stop("no data loaded")
  }
  
  # plot raster image
  mat_raster <- raster(project$data_processed$data_coverage)
  raster::plot(mat_raster, col = plasma(100), axes = FALSE, box = FALSE,
               xlab = "locus", ylab = "sample", legend.args = list(text = "Coverage"))
}

#------------------------------------------------
#' @title Plot raster image of within-sample allele frequencies
#'
#' @description Plot raster image of within-sample allele frequencies.
#'
#' @param project a \code{mipmapper_project} with data already loaded
#'
#' @export

plot_WSAF_matrix <- function(project) {
  
  # check inputs
  assert_custom_class(project, "mipmapper_project")
  if (is.null(project$data_processed$data_coverage)) {
    stop("no data loaded")
  }
  
  # plot raster image
  mat_raster <- raster(project$data_processed$data_WSAF)
  raster::plot(mat_raster, col = plasma(100), axes = FALSE, box = FALSE,
               xlab = "locus", ylab = "sample", legend.args = list(text = "WSAF"))
  
}

#------------------------------------------------
#' @title Plot PCA
#'
#' @description Plots either the first 2 or 3 principal components, with the
#'   data points labelled accordng to chosen meta data.
#'
#' @param project a \code{mipmapper_project} containing pca analysis output
#' @param num_components numeric for number of components used
#' @param grouping_var character for the desired meta variable to be used for 
#'   labelling the scatterplot
#'
#' @export

plot_pca <- function(project, num_components = 2, grouping_var = NULL) {
  
  # check inputs
  assert_custom_class(project, "mipmapper_project")
  if (length(project$analyses$pca) == 0) {
    stop("project must have a pca analysis already completed")
  }
  assert_in(num_components, c(2,3))
  if (!is.null(grouping_var)) {
    assert_single_string(grouping_var)
    assert_in(grouping_var, names(project$data_processed$data_sample))
    assert_string(project$data_processed$data_sample[,grouping_var])
  }
  
  # check number of components
  pca <- project$analyses$pca
  n <- ncol(pca$x)
  
  # create colour vector
  group <- rep("",n)
  if (!is.null(grouping_var)) {
    group <- project$data_processed$data_sample[,grouping_var]
  }
  col_vec <- rev(divergent_colours(length(unique(group))))
  col_vec <- setNames(col_vec, unique(group))
  
  # scatterplot of first 2 principal components
  if (num_components == 2) {
    ret <- plot_ly(as.data.frame(pca$x), x = ~PC1, y = ~PC2,
                   color = group, type = "scatter", colors = col_vec,
                   mode = "markers", marker = list(size = 5))
  }
  
  # scatterplot of first 3 principal components
  if (num_components == 3) {
    ret <- plot_ly(as.data.frame(pca$x), x = ~PC1, y = ~PC2, z = ~PC3,
                   color = group, type = "scatter3d", colors = col_vec,
                   mode = "markers", marker = list(size = 3))
  }
  
  # return
  return(ret)
}

#------------------------------------------------
#' @title Plot variance explained by PCA components
#'
#' @description Plot the variance explained by each component. The number of 
#'   components shown is controlled by \code{num_components}, with up to the
#'   first 10 components shown by default.
#'
#' @param project a \code{mipmapper_project} containing pca analysis output
#' @param num_components number of components to be shown
#'
#' @export

plot_pca_variance <- function(project, num_components = 10) {
  
  # check inputs
  assert_custom_class(project, "mipmapper_project")
  if (length(project$analyses$pca) == 0) {
    stop("project must have a pca analysis already completed")
  }
  assert_single_pos_int(num_components, zero_allowed = FALSE)
  
  # check number of components
  pca <- project$analyses$pca
  nc <- min(length(colnames(pca$x)), num_components)
  
  # x has to be factored in order to preserve PC order
  x <- factor(colnames(pca$x)[seq_len(nc)],
              levels = colnames(pca$x)[seq_len(nc)])
  
  ret <- plot_ly(x = x, y = pca$var[seq_len(nc)], type = "bar") %>%
    plotly::layout(title = "Percentage of total variance explaned by PCA",
                   xaxis = list(title = "Principal Component"),
                   yaxis = list(title = "Percentage of variance explained"))
  
  # return
  return(ret)
}

#------------------------------------------------
#' @title Plot loading values of PCA
#'
#' @description Plot loading values of PCA for a given component.
#'
#' @param project a \code{mipmapper_project} containing pca analysis output
#' @param component which principal component to plot
#' @param colour optional name of variable to colour by
#' @param plot_resistance whether to overlay known drug resistance loci
#'
#' @export

plot_pca_loadings <- function(project, component = 1, colour = NULL, plot_resistance = TRUE) {
  
  # check inputs
  assert_custom_class(project, "mipmapper_project")
  if (length(project$analyses$pca) == 0) {
    stop("project must have a pca analysis already completed")
  }
  assert_single_pos_int(component, zero_allowed = FALSE)
  assert_single_string(colour)
  assert_single_logical(plot_resistance)
  
  # create dataframe of locus info and loadings
  df <- project$data_processed$data_locus
  loading <- project$analyses$pca$loadings[,component]
  df$Loading <- loading/max(loading)
  
  # add colour variable
  if (!is.null(colour)) {
    assert_in(colour, names(df))
    df$Colour <- df[,colour]
  }
  
  # add chromosome lengths to df
  chrom_length <- c(643292, 947102, 1060087, 1204112, 1343552, 1418244, 1501717, 1419563, 1541723, 1687655, 2038337, 2271478, 2895605, 3291871)
  df$Chrom_length <- chrom_length[df$Chrom]
  
  # create basic empty plot
  plot_load <- ggplot()
  plot_load <- plot_load + theme(axis.text.x = element_blank(),
                                 panel.background = element_blank(),
                                 panel.grid.major = element_blank(), 
                                 panel.grid.minor = element_blank())
  plot_load <- plot_load + scale_y_continuous(breaks = 1:14+0.5, labels=14:1)
  plot_load <- plot_load + xlab("position") + ylab("chromosome")
  
  # add grey rect for each chromosome
  plot_load <- plot_load  + geom_rect(aes_(xmin = 0, xmax = ~Chrom_length, ymin = ~(15-Chrom), ymax = ~(15-Chrom + 0.95)), fill = grey(0.97), colour = grey(0.5), size = 0.2, data = df)
  
  # add loci
  if (is.null(colour)) {
    plot_load <- plot_load  + geom_segment(aes_(x = ~Pos, xend = ~Pos, y = ~(15-Chrom), yend = ~(15-Chrom + 0.95*Loading)), data = df)
  } else {
    plot_load <- plot_load  + geom_segment(aes_(x = ~Pos, xend = ~Pos, y = ~(15-Chrom), yend = ~(15-Chrom + 0.95*Loading), colour = ~Colour), data = df)
    plot_load <- plot_load + guides(colour = guide_legend(title = colour))
  }
  
  # optionally overlay known resistance loci
  if (plot_resistance) {
    
    # create separate data frame
    df_resistance <- data.frame(Name = c("kelch K13", "pfcrt", "pfmdr1", "pfmdr2", "pfmrp1", "pfdhps", "pfdhfr", "pfaat1"), Chrom = c(13, 7, 5, 14, 1, 8, 4, 6), Pos = c(1724817, 403222, 957890, 1954601, 464726, 548200, 748088, 1214976))
    
    # add to plot
    plot_load <- plot_load  + geom_point(aes_(x = ~Pos, y = ~(15-Chrom + 0.5), text = ~Name), size = 0.5, data = df_resistance)
  }
  
  # plot using plotly
  plotly::ggplotly(plot_load, tooltip = "text")
}

#------------------------------------------------
#' @title Plot points on map of DRC
#'
#' @description Plot points on map of DRC
#'
#' @param longitude x-coordinates of points to plot
#' @param latitude y-coordinates of points to plot
#' @param v plotting value
#' @param longitude_limits plotting x-limits
#' @param latitude_limits plotting y-limits
#' @param vmin minimum value used in colour scale
#' @param vmax maximum value used in colour scale
#' @param plot_water whether to overlay major water bodies
#' @param raw_cols the cols that make up the colour palette
#'
#' @export

plot_DRC <- function(longitude = NULL, latitude = NULL, v = NULL, longitude_limits = c(10,40), latitude_limits = c(-15,7), vmin = min(v, na.rm = TRUE), vmax = max(v, na.rm = TRUE), plot_water = TRUE, raw_cols = viridis(100)) {
  
  # check inputs
  if (!is.null(longitude)) {
    assert_numeric(longitude)
  }
  if (!is.null(latitude)) {
    assert_numeric(latitude)
  }
  if (!is.null(v)) {
    assert_numeric(v)
  }
  assert_same_length_multiple(longitude, latitude, v)
  assert_numeric(longitude_limits)
  assert_length(longitude_limits, 2)
  assert_numeric(latitude_limits)
  assert_length(latitude_limits, 2)
  assert_single_numeric(vmin)
  assert_single_numeric(vmax)
  assert_gr(vmax, vmin)
  assert_single_logical(plot_water)
  
  # load country border shapefiles
  country_borders <- mipmapper_file("country_borders.rds")
  
  # empty plot
  plot(0, type = "n", xlim = longitude_limits, ylim = latitude_limits, xaxs = "i", yaxs = "i",
       xlab = "longitude", ylab = "latitude")
  
  # add country shapes
  country_codes <- c("ANG", "BUR", "CAM", "CAR", "CHA", "CNG", "EQG", "ETH", "GAB", "KEN", "MAA", "MOZ", "NIR", "RWA", "SUD", "TAN", "UGA", "ZAI", "ZAM")
  for (i in 1:length(country_codes)) {
    plot(country_borders[[country_codes[i]]], col = grey(0.2), border = "white", add = TRUE)
  }
  
  # add major water bodies
  if (plot_water) {
    
    # load shapefile
    DRC_water <- mipmapper_file("DRC_water.rds")
    
    # add to plot
    plot(DRC_water, col = grey(0.4), border = NA, add = TRUE)
  }
  
  # add points
  if (!is.null(longitude)) {
    
    # convert points to unit scale
    v_scale <- (v-vmin)/(vmax-vmin)
    
    # add to plot
    my_pal <- colorRampPalette(raw_cols)
    col_vec <- my_pal(101)[floor(v_scale*100)+1]
    points(longitude, latitude, col = col_vec, pch = 20, cex = 0.5)
  }
  
}
