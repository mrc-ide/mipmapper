# run this after pkgdown::build_site to correct image links
# fix link

fix_pkgdown <- function(){

# first knit the Rmd
knitr::knit("README.Rmd")

# then build the site
pkgdown::build_site()

# correct pkgdowns
lines <- readLines("README.md")
lines[grep("(.*png)",lines)] <- gsub("!\\[\\]\\(R","!\\[\\](tools/R",lines[grep("(.*png)",lines)])
writeLines(lines, "README.md")
file.copy(grep("png",list.files(),value = TRUE),paste0("tools/",grep("png",list.files(),value = TRUE)), overwrite = TRUE)
file.copy(grep("png",list.files(),value = TRUE),paste0("docs/",grep("png",list.files(),value = TRUE)), overwrite = TRUE)
unlink(grep("png",list.files(),value = TRUE))

dir.create("docs/tools")
file.copy("tools","docs",recursive = TRUE, overwrite = TRUE)

l <- readLines("docs/index.html")
fun <- grep(".png",l,value=TRUE, fixed=TRUE)
files <- strsplit(fun,"/|\"") %>% lapply(function(x) grep("png",x,value=TRUE)) %>% unlist
for(i in 1:length(fun)) {
 files[i] <- gsub("\"(.*png)\"",paste0("\"","tools/",files[i],"\""),fun[i])
}
l[grepl(".png",l, fixed=TRUE)] <- files
writeLines(l,"docs/index.html")

}
