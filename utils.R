save.xlsx <- function (file, ...)
{
  require(xlsx, quietly = TRUE)
  objects <- list(...)
  fargs <- as.list(match.call(expand.dots = TRUE))
  objnames <- as.character(fargs)[-c(1, 2)]
  nobjects <- length(objects)
  for (i in 1:nobjects) {
    if (i == 1)
      write.xlsx(objects[[i]], file, sheetName = objnames[i], row.names = FALSE)
    else write.xlsx(objects[[i]], file, sheetName = objnames[i],
                    append = TRUE, row.names = FALSE)
  }
  print(paste("Workbook", file, "has", nobjects, "worksheets."))
}

install_and_load_libraries <- function(){
  
  # CRAN
  list.of.packages <- c("ggplot2", "reshape2","doParallel", "pheatmap", "igraph", "seqinr", "foreach", "plotly", "rtracklayer")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  # bioconductor
  # source("https://bioconductor.org/biocLite.R")
  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install()
  
  list.of.packages <- c("Biostrings", "VariantAnnotation","BSgenome")
  #new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  #if(length(new.packages)) biocLite(new.packages)
  
  BiocManager::install(list.of.packages)
  
  require(seqinr)
  require(Biostrings)
  require(VariantAnnotation)
  require(foreach)
  require(doParallel)
  require(ggplot2)
  require(BSgenome)
  require(VariantAnnotation)
  require(plotly)
  require(rtracklayer)
}
