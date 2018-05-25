
# run this after pkgdown::build_site to correct image links

# change lines in README.md to make link point to tools
lines <- readLines("README.md")
lines[grep("(.*png)",lines)] <- gsub("!\\[\\]\\(R","!\\[\\](tools/R",lines[grep("(.*png)",lines)])
writeLines(lines, "README.md")

# get all png files in root directory
png_files <- grep("png",list.files(),value = TRUE)

# copy into tools/ and into docs/. Delete original
file.copy(png_files ,paste0("tools/", png_files), overwrite = TRUE)
file.copy(png_files, paste0("docs/", png_files), overwrite = TRUE)
unlink(png_files)

# copy everything from tools to docs/tools
dir.create("docs/tools", showWarnings = FALSE)
file.copy("tools","docs",recursive = TRUE, overwrite = TRUE)

# change lines in infex.html to point to tools
l <- readLines("docs/index.html")
fun <- grep(".png",l,value=TRUE, fixed=TRUE)
files <- strsplit(fun,"/|\"") %>% lapply(function(x) grep("png",x,value=TRUE)) %>% unlist
for(i in 1:length(fun)) {
  files[i] <- gsub("\"(.*png)\"",paste0("\"","tools/",files[i],"\""),fun[i])
}
l[grepl(".png",l, fixed=TRUE)] <- files
writeLines(l,"docs/index.html")

