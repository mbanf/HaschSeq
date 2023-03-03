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


write_bed <- function(df.QTLs){
  
  df.QTLs.bed <- df.QTLs[,c("B73-chr", "B73-pos", "B73-pos")]
  colnames(df.QTLs.bed) <- c('chrom', 'chromStart', 'chromEnd')
  
  for(i in 1:nrow(df.QTLs)){
    if(df.QTLs$POSTfreq[i] > 0.5){
      df.QTLs.bed$chrom[i] = df.QTLs$`B73-chr`[i]
      df.QTLs.bed$chromStart[i] = df.QTLs$`B73-pos`[i] - 1
      df.QTLs.bed$chromEnd[i] = df.QTLs$`B73-pos`[i]
      
    }else if(df.QTLs$POSTfreq[i] < 0.5){
      df.QTLs.bed$chrom[i] = df.QTLs$`Mo17-chr`[i]
      df.QTLs.bed$chromStart[i] = df.QTLs$`Mo17-pos`[i] - 1
      df.QTLs.bed$chromEnd[i] = df.QTLs$`Mo17-pos`[i]
    }
  }
  
  df.QTLs.bed
}


extract_sequences <- function(df.bQTLs,
                              genome_sequences,
                              s.env_bps = 100,
                              n.chromosomes = 10
){
  
  df.bQTLs["id"] <- seq(1, nrow(df.bQTLs))
  si <- seqinfo(genome_sequences)
  message("extract sequences surrounding ASBs from diploid genome")
  
  v.sequences <- character(length = nrow(df.bQTLs))
  
  for(i in 1:n.chromosomes){ 
    print(i)
    
    chr_name <- paste("B73-chr", i, sep = "")
    chr_seq <- getSeq(genome_sequences, GRanges(chr_name, IRanges(1, seqlengths(si)[chr_name])))
    
    df.bQTLs.B73.i <- subset(df.bQTLs, df.bQTLs$`B73-chr` == chr_name)
    df.bQTLs.B73.i <- subset(df.bQTLs.B73.i, df.bQTLs.B73.i$POSTfreq > 0.5)
    
    for(j in 1:nrow(df.bQTLs.B73.i)){
      id <- df.bQTLs.B73.i$id[j]
      v.sequences[id] <- subseq(chr_seq, start=df.bQTLs.B73.i$`B73-pos`[j] - s.env_bps, end=df.bQTLs.B73.i$`B73-pos`[j] + s.env_bps)  
    }
    
    chr_name <- paste("Mo17-chr", i, sep = "")
    chr_seq <- getSeq(genome_sequences, GRanges(chr_name, IRanges(1, seqlengths(si)[chr_name])))
    
    df.bQTLs.Mo17.i <- subset(df.bQTLs, df.bQTLs$`Mo17-chr` == chr_name)
    df.bQTLs.Mo17.i <- subset(df.bQTLs.Mo17.i, df.bQTLs.Mo17.i$POSTfreq < 0.5)
    
    for(j in 1:nrow(df.bQTLs.Mo17.i)){
      id <- df.bQTLs.Mo17.i$id[j]
      v.sequences[id] <- subseq(chr_seq, start=df.bQTLs.Mo17.i$`Mo17-pos`[j] - s.env_bps, end=df.bQTLs.Mo17.i$`Mo17-pos`[j] + s.env_bps)  
    }
  }
  
  descs <-  paste("asb", df.bQTLs$id, sep = "")
  asb_environment <- DNAStringSet(v.sequences)
  names(asb_environment) <- descs
  
  asb_environment
}


load_gene_annotation <- function(path.gene_annotation){
  
  df.gff <- read.table(path.gene_annotation, header = FALSE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE, quote = "")
  v.genePartitions <- c("gene", "five_prime_UTR", "CDS", "three_prime_UTR", "exon")
  df.gff <- subset(df.gff, df.gff$V3 %in% v.genePartitions)
  
  ids.list <- sapply(df.gff$V9, function(m) {strsplit(m, ";")})
  
  df.gff["id"] <- sapply(df.gff$V9, function(m) {gsub("ID=", "", unlist(strsplit(m, ";"))[1])})
  df.gff$id <- gsub("gene:","", df.gff$id)
  df.gff$id <- gsub("Parent=","", df.gff$id)
  df.gff$id <- gsub("transcript:","", df.gff$id) 
  df.gff$id <- gsub("CDS:","", df.gff$id)
  df.gff$id <- gsub("\\_.*","", df.gff$id)
  
  names(df.gff) <- c("chr", "source", "partition", "pos.start", "pos.stop", "aux1", "strand", "aux2", "gene_meta", "gene.ID")
  
  df.gff
}


bQTL_scatterplot <- function(df.bQTLs, n.chromosomes = 10){
  
  ### scatter plot of the postfrequencies - should be updated to > , < 0.5
  v.offset <- numeric(n.chromosomes)
  v.start <- numeric(n.chromosomes)
  
  for(i in 1:n.chromosomes){
    chr <- paste("B73-chr", i, sep = "")
    df.bQTLs.i <- subset(df.bQTLs, df.bQTLs$`B73-chr` == chr)
    v.offset[i] <-  max(df.bQTLs.i$`B73-pos`)
    v.start[i] <- sum(v.offset[1:i])
  }
  
  v.start <- c(0, v.start)
  
  # directionality plot 
  df.scatterplot <- subset(df.bQTLs, df.bQTLs$`B73-chr` == "B73-chr1")
  
  for(i in 1:(n.chromosomes - 1)){
    chr <- paste("B73-chr", (i+1), sep = "")
    df.bQTLs.i <- subset(df.bQTLs, df.bQTLs$`B73-chr` == chr)
    df.bQTLs.i$`B73-pos` <- df.bQTLs.i$`B73-pos` + sum(v.offset[1:i])
    df.scatterplot <- rbind(df.scatterplot, df.bQTLs.i)
  }
  
  p6 <- ggplot(df.scatterplot, aes(x = `B73-pos`, y = POSTfreq, fill = POSTfreq)) +
    geom_point(shape = 21, size = 1,  alpha = 0.5, stroke = 0.0) + theme_bw() 
  
  p6 <- p6 + scale_fill_continuous(low = "blue", high = "green2")
  p6 <- p6 + geom_vline(xintercept = v.start, col='black', lwd=0.5, linetype="dashed")
  
  p6
}



genomic_distribution <- function(df.ASBs, 
                                 df.bgSNPs, 
                                 v.partitions){
  
  df.distribution <- data.frame(ASBs = apply(df.ASBs[,v.partitions], 2, table)["yes",],
                                bgSNPs = apply(df.bgSNPs[,v.partitions], 2, table)["yes",])
  df.distribution["%ASBs"] <-  round(df.distribution$ASBs / sum(df.distribution$ASBs) * 100, 1)
  df.distribution["%bgSNPs"] <-  round(df.distribution$bgSNPs / sum(df.distribution$bgSNPs) * 100, 1)
  df.distribution["multiplier"] <- round(df.distribution$bgSNPs / df.distribution$ASBs, 1)
  
  print(df.distribution)
  
}


define_path <- function(folder, file){
  paste(folder, file, sep = "/")
}

# 
# bQTL_scatterplot <- function(postTotal=postTotal){
#   
#   #  postTotal <- l.selection[[1]]
#   
#   ### scatter plot of the postfrequencies
#   v.offset <- numeric(n.chromosomes)
#   v.start <- numeric(n.chromosomes)
#   
#   for(i in 1:n.chromosomes){
#     postTotal.i <- subset(postTotal, postTotal$contig == i)
#     v.offset[i] <-  max(postTotal.i$position)
#     v.start[i] <- sum(v.offset[1:i])
#   }
#   
#   v.start <- c(0, v.start)
#   
#   # v.offset <- v.offset[]
#   
#   # directionality plot 
#   df.scatterplot <- subset(postTotal, postTotal$contig == 1)
#   
#   
#   for(i in 1:(n.chromosomes - 1)){
#     postTotal.i <- subset(postTotal, postTotal$contig == (i + 1))
#     
#     postTotal.i$position <- postTotal.i$position + sum(v.offset[1:i])
#     
#     df.scatterplot <- rbind(df.scatterplot, postTotal.i)
#   }
#   
#   
#   p6 <- ggplot(df.scatterplot, aes(x = position, y = POSTfreq, fill = POSTfreq)) +
#     geom_point(shape = 21, size = 1,  alpha = 0.5, stroke = 0.0) + theme_bw() 
#   
#   p6 <- p6 + scale_fill_continuous(low = "blue", high = "green2")
#   p6 <- p6 + geom_vline(xintercept = v.start, col='black', lwd=0.5, linetype="dashed")
#   
#   # plot as pdf
#   #(p6 <- ggplotly(p6))
#   
#   p6
# }
# 
# # only display selection - figure paper 
# bQTL_scatterplot_chr <- function(postTotal=postTotal, chr = 3){
#   
#   #  postTotal <- l.selection[[1]]
#   
#   ### scatter plot of the postfrequencies
#   v.offset <- numeric(n.chromosomes)
#   v.start <- numeric(n.chromosomes)
#   
#   for(i in 1:n.chromosomes){
#     postTotal.i <- subset(postTotal, postTotal$contig == i)
#     v.offset[i] <-  max(postTotal.i$position)
#     v.start[i] <- sum(v.offset[1:i])
#   }
#   
#   v.start <- c(0, v.start)
#   
#   # v.offset <- v.offset[]
#   
#   # directionality plot 
#   df.scatterplot <- subset(postTotal, postTotal$contig == chr)
#   
#   
#   
#   # 
#   # for(i in 1:(n.chromosomes - 1)){
#   #   postTotal.i <- subset(postTotal, postTotal$contig == (i + 1))
#   #   # remove artifacts
#   #   if(FALSE){
#   #     message("artifact removal")
#   #     if(i == 2){
#   #       df.postTotal_both.selection.i <- subset(df.postTotal_both.selection, df.postTotal_both.selection$contig == 3)
#   #       postTotal.i <- subset(postTotal.i, postTotal.i$position < min(df.postTotal_both.selection.chromosome.3$position) |
#   #                                          postTotal.i$position > max(df.postTotal_both.selection.chromosome.3$position))
#   # 
#   #     }
#   #   }
#   # 
#   #   postTotal.i$position <- postTotal.i$position + sum(v.offset[1:i])
#   # 
#   #   df.scatterplot <- rbind(df.scatterplot, postTotal.i)
#   # }
#   
#   
#   p6 <- ggplot(df.scatterplot, aes(x = position, y = POSTfreq, fill = POSTfreq)) +
#     geom_point(shape = 21, size = 1,  alpha = 0.5, stroke = 0.0) + theme_bw() 
#   
#   p6 <- p6 + scale_fill_continuous(low = "blue", high = "green2")
#   # p6 <- p6 + geom_vline(xintercept = v.start, col='black', lwd=0.5, linetype="dashed")
#   
#   # plot as pdf
#   #(p6 <- ggplotly(p6))
#   
#   p6
# }
# 
