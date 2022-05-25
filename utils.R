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


bQTL_scatterplot <- function(postTotal=postTotal){
  
  #  postTotal <- l.selection[[1]]
  
  ### scatter plot of the postfrequencies
  v.offset <- numeric(n.chromosomes)
  v.start <- numeric(n.chromosomes)
  
  for(i in 1:n.chromosomes){
    postTotal.i <- subset(postTotal, postTotal$contig == i)
    v.offset[i] <-  max(postTotal.i$position)
    v.start[i] <- sum(v.offset[1:i])
  }
  
  v.start <- c(0, v.start)
  
  # v.offset <- v.offset[]
  
  # directionality plot 
  df.scatterplot <- subset(postTotal, postTotal$contig == 1)
  
  
  for(i in 1:(n.chromosomes - 1)){
    postTotal.i <- subset(postTotal, postTotal$contig == (i + 1))
    
    postTotal.i$position <- postTotal.i$position + sum(v.offset[1:i])
    
    df.scatterplot <- rbind(df.scatterplot, postTotal.i)
  }
  
  
  p6 <- ggplot(df.scatterplot, aes(x = position, y = POSTfreq, fill = POSTfreq)) +
    geom_point(shape = 21, size = 1,  alpha = 0.5, stroke = 0.0) + theme_bw() 
  
  p6 <- p6 + scale_fill_continuous(low = "blue", high = "green2")
  p6 <- p6 + geom_vline(xintercept = v.start, col='black', lwd=0.5, linetype="dashed")
  
  # plot as pdf
  #(p6 <- ggplotly(p6))
  
  p6
}

# only display selection - figure paper 
bQTL_scatterplot_chr <- function(postTotal=postTotal, chr = 3){
  
  #  postTotal <- l.selection[[1]]
  
  ### scatter plot of the postfrequencies
  v.offset <- numeric(n.chromosomes)
  v.start <- numeric(n.chromosomes)
  
  for(i in 1:n.chromosomes){
    postTotal.i <- subset(postTotal, postTotal$contig == i)
    v.offset[i] <-  max(postTotal.i$position)
    v.start[i] <- sum(v.offset[1:i])
  }
  
  v.start <- c(0, v.start)
  
  # v.offset <- v.offset[]
  
  # directionality plot 
  df.scatterplot <- subset(postTotal, postTotal$contig == chr)
  
  
  
  # 
  # for(i in 1:(n.chromosomes - 1)){
  #   postTotal.i <- subset(postTotal, postTotal$contig == (i + 1))
  #   # remove artifacts
  #   if(FALSE){
  #     message("artifact removal")
  #     if(i == 2){
  #       df.postTotal_both.selection.i <- subset(df.postTotal_both.selection, df.postTotal_both.selection$contig == 3)
  #       postTotal.i <- subset(postTotal.i, postTotal.i$position < min(df.postTotal_both.selection.chromosome.3$position) |
  #                                          postTotal.i$position > max(df.postTotal_both.selection.chromosome.3$position))
  # 
  #     }
  #   }
  # 
  #   postTotal.i$position <- postTotal.i$position + sum(v.offset[1:i])
  # 
  #   df.scatterplot <- rbind(df.scatterplot, postTotal.i)
  # }
  
  
  p6 <- ggplot(df.scatterplot, aes(x = position, y = POSTfreq, fill = POSTfreq)) +
    geom_point(shape = 21, size = 1,  alpha = 0.5, stroke = 0.0) + theme_bw() 
  
  p6 <- p6 + scale_fill_continuous(low = "blue", high = "green2")
  # p6 <- p6 + geom_vline(xintercept = v.start, col='black', lwd=0.5, linetype="dashed")
  
  # plot as pdf
  #(p6 <- ggplotly(p6))
  
  p6
}

