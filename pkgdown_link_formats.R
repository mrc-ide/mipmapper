# run this after pkgdown::build_site to correct image links
# fix link
lines <- readLines("README.md")
lines[grep("(.*png)",lines)] <- gsub("!\\[\\]\\(","!\\[\\](tools/",lines[grep("(.*png)",lines)])
writeLines(lines, "README.md")
file.copy(grep("png",list.files(),value = TRUE),paste0("tools/",grep("png",list.files(),value = TRUE)))
file.copy(grep("png",list.files(),value = TRUE),paste0("docs/",grep("png",list.files(),value = TRUE)))
unlink(grep("png",list.files(),value = TRUE))