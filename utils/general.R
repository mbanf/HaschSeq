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
  BiocManager::install(version = "3.11")
  
  #list.of.packages <- c("Biostrings", "VariantAnnotation","BSgenome")
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



get_GOSlim_annotations <- function(filename=filename){
  library(GO.db)
  df.GO.annot <- read.table(filename, header = FALSE, sep = "\t", na.strings=NA, stringsAsFactors = FALSE)
  names(df.GO.annot) <- c("acc.", "go")
  df.GO.annot["Term"] <- NA
  df.GO.annot["Ontology"] <- NA
  for(j in 1:nrow(df.GO.annot)){
    cat("Processing... ", round(j/nrow(df.GO.annot) * 100, digits = 2) , "%", "\r") 
    flush.console()
    if(!is.null(GOTERM[[df.GO.annot$go[j]]])){
      df.GO.annot$Term[j] <- Term(GOTERM[[df.GO.annot$go[j]]])
      df.GO.annot$Ontology[j] <- Ontology(GOTERM[[df.GO.annot$go[j]]])
    }
  }
  df.GO.annot$acc. <- gsub("\\|.*", "", df.GO.annot$acc.)
  return(df.GO.annot)
}

create_BSGenome_AGP4 <- function(){
  
  homeFolder <- "/Users/michaelbanf/Documents/postdoctoral_work/Projects/HashSeq/AGPv4/"
  setwd(homeFolder)
  
  # setwd("/Users/michaelbanf/Documents/postdoctoral_work/Projects/")
  # setwd("/shared/Labs/Rhee Lab/Everyone/Michael/Thomas/hybrid_chip/MaizeGenome4/")
  #source("https://bioconductor.org/biocLite.R")
  #biocLite("BSgenome")
  
  library(BSgenome)
  forgeBSgenomeDataPkg("zmays_seed")
  
  # install.packages("BSgenome.Zmays.NCBI.AGPv4", repos = NULL, type="source")
  # assert ? 
  
  homeFolder <- "/Users/michaelbanf/Documents/postdoctoral_work/Projects/HashSeq/"
  setwd(homeFolder)
  
  install.packages("BSgenome.Zmays.NCBI.AGPv4", repos = NULL, type="source")
  library("BSgenome.Zmays.NCBI.AGPv4")
  
  
  
  # snps einlesen als vcf 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}



# df.snp_positions - needs to have identical syntax
addNearestGene <- function(df.snp_positions, pos_peak, chr_Nr,  df.gene_annotation = df.gene_annotation, 
                           df.gene_annotation.AGPv3 = NULL, 
                           df.gene_function = NULL){
  
  postTotal.significant <- df.ChipSeq.filtered
  postTotal.significant <- l.selection[[1]]
  postTotal.significant <- df.snp_positions
  postTotal.significant <- df.peaks
  
  
  df.gene_conversion.AGPv3_to_AGPv4 <- df.geneID_conversion
  df.gene_function <- df.gene_function
  
  # ASBs
  pos_peak = "position" 
  chr_Nr = "contig"
  
  # Peaks
  pos_peak = "posPeak" 
  chr_Nr = "seqnames"
  
  # df.ChipSeq.filtered
  pos_peak = "pos" 
  chr_Nr = "chr"
  
  
  idx.pos <- which( colnames(postTotal.significant)== pos_peak)
  idx.chr <- which( colnames(postTotal.significant) == chr_Nr)
  
  
  message("identifying gene partitions - positions of ASBs with ")
  
  df.bQTL_gene_partitioning <- c()
  
  strt<-Sys.time() 
  cl<-makeCluster(min(n.chromosomes, n.cpus))
  registerDoParallel(cl)
  
  l.bQTL_gene_partitioning_per_chromosome <-  foreach(i = 1:n.chromosomes, .packages=c("seqinr", "VariantAnnotation", "Biostrings")) %dopar% { 
    #for(i in 1:n.chromosomes){
    
    message(paste("processing chromosome", i))
    
    postTotal.significant.i <- subset(postTotal.significant, postTotal.significant[,idx.chr] == i)
    
    postTotal.significant.i["promoter_5kb"] <- "no"
    postTotal.significant.i["promoter_1kb"] <- "no"
    postTotal.significant.i["gene"] <- "no"
    postTotal.significant.i["five_prime_UTR"] <- "no"
    postTotal.significant.i["exon"] <- "no"
    postTotal.significant.i["intron"] <- "no"
    postTotal.significant.i["three_prime_UTR"] <- "no"
    postTotal.significant.i["post_gene_1kb"] <- "no"
    postTotal.significant.i["post_gene_5kb"] <- "no"
    
    postTotal.significant.i["gene.ID"] <- NA
    
    df.gff.i <- subset(df.gene_annotation, df.gene_annotation$chr == i)
    
    l.df.gff.i <- vector(mode = "list", length = length(v.genePartitions))
    l.df.gff.i[[1]] <- subset(df.gff.i, df.gff.i$partition %in% v.genePartitions[1])
    l.df.gff.i[[2]] <- subset(df.gff.i, df.gff.i$partition %in% v.genePartitions[2])
    l.df.gff.i[[3]] <- subset(df.gff.i, df.gff.i$partition %in% v.genePartitions[3])
    l.df.gff.i[[4]] <- subset(df.gff.i, df.gff.i$partition %in% v.genePartitions[4])
    l.df.gff.i[[5]] <- subset(df.gff.i, df.gff.i$partition %in% v.genePartitions[5])
    
    # preselect all in gene 
    pb <- txtProgressBar(min = 0, max = nrow(postTotal.significant.i), style = 3)
    
    for(j in 1:nrow(postTotal.significant.i)){
      
      setTxtProgressBar(pb, j)
      
      # in gene 
      i.set <- which(l.df.gff.i[[1]]$pos.start < postTotal.significant.i[j,idx.pos] & postTotal.significant.i[j,idx.pos] < l.df.gff.i[[1]]$pos.stop)
      
      if(length(i.set) > 0){
        
        postTotal.significant.i$gene[j] <- "yes"
        postTotal.significant.i$gene.ID[j] <- l.df.gff.i[[1]]$gene.ID[i.set[1]]
        
        # in 5 prime utr
        i.set <- which(l.df.gff.i[[2]]$pos.start < postTotal.significant.i[j,idx.pos] & postTotal.significant.i[j,idx.pos] < l.df.gff.i[[2]]$pos.stop)
        if(length(i.set) > 0){
          postTotal.significant.i$five_prime_UTR[j] <- "yes"
        }else{
          # in 3 prime utr
          i.set <- which(l.df.gff.i[[4]]$pos.start < postTotal.significant.i[j,idx.pos] & postTotal.significant.i[j,idx.pos] < l.df.gff.i[[4]]$pos.stop)
          if(length(i.set) > 0){
            postTotal.significant.i$three_prime_UTR[j] <- "yes"
          }else{
            # in CDS (exon)
            i.set <- which(l.df.gff.i[[5]]$pos.start < postTotal.significant.i[j,idx.pos] & postTotal.significant.i[j,idx.pos] < l.df.gff.i[[5]]$pos.stop)
            if(length(i.set) > 0){
              postTotal.significant.i$exon[j] <- "yes"
            }else{
              # i.set <- which(l.df.gff.i[[3]]$pos.start < postTotal.significant.i$position[j] & postTotal.significant.i$position[j] < l.df.gff.i[[3]]$pos.stop)
              #if(length(i.set) > 0){  
              postTotal.significant.i$intron[j] <- "yes"
              #}
            }
          } 
        }
      }
    }
    close(pb)
    
    # identify bQTL in 5kb and 1kb distance 
    l.df.gff.i[[1]]["up_5kb"] <- 0
    l.df.gff.i[[1]]["up_1kb"] <- 0
    l.df.gff.i[[1]]["start_gene"] <- 0
    l.df.gff.i[[1]]["end_gene"] <- 0
    l.df.gff.i[[1]]["down_1kb"] <- 0
    l.df.gff.i[[1]]["down_5kb"] <- 0
    
    for(k in 1:nrow(l.df.gff.i[[1]])){
      if(l.df.gff.i[[1]]$strand[k] == "+"){
        l.df.gff.i[[1]]$up_5kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k] - 5000)
        l.df.gff.i[[1]]$up_1kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k] - 1000)
        l.df.gff.i[[1]]$start_gene[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k])
        l.df.gff.i[[1]]$end_gene[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k])
        l.df.gff.i[[1]]$down_1kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k] + 1000)
        l.df.gff.i[[1]]$down_5kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k] + 5000)
      }else if(l.df.gff.i[[1]]$strand[k] == "-"){
        l.df.gff.i[[1]]$up_1kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k] + 1000)
        l.df.gff.i[[1]]$up_5kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k] + 5000)
        l.df.gff.i[[1]]$start_gene[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k])
        l.df.gff.i[[1]]$end_gene[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k])
        l.df.gff.i[[1]]$down_1kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k] - 1000)
        l.df.gff.i[[1]]$down_5kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k] - 5000)
      }
    }
    
    
    
    i.set_upDown <- which(postTotal.significant.i$gene == "no")
    
    pb <- txtProgressBar(min = 0, max = length(i.set_upDown), style = 3)
    
    for(j in 1:length(i.set_upDown)){
      
      setTxtProgressBar(pb, j)
      
      # in 1 KB upstream 
      i.set <- which(l.df.gff.i[[1]]$up_1kb < postTotal.significant.i[i.set_upDown[j],idx.pos] & postTotal.significant.i[i.set_upDown[j],idx.pos] < l.df.gff.i[[1]]$start_gene 
                     | l.df.gff.i[[1]]$start_gene < postTotal.significant.i[i.set_upDown[j],idx.pos] & postTotal.significant.i[i.set_upDown[j],idx.pos] < l.df.gff.i[[1]]$up_1kb)
      
      if(length(i.set)  > 0){
        postTotal.significant.i$promoter_1kb[i.set_upDown[j]] <- "yes"
        postTotal.significant.i$gene.ID[i.set_upDown[j]] <- l.df.gff.i[[1]]$gene.ID[i.set[1]]
      }else{
        # in 5 KB upstream 
        i.set <- which(l.df.gff.i[[1]]$up_5kb < postTotal.significant.i[i.set_upDown[j],idx.pos] & postTotal.significant.i[i.set_upDown[j],idx.pos] < l.df.gff.i[[1]]$start_gene
                       | l.df.gff.i[[1]]$start_gene < postTotal.significant.i[i.set_upDown[j],idx.pos] & postTotal.significant.i[i.set_upDown[j],idx.pos] < l.df.gff.i[[1]]$up_5kb)
        if(length(i.set)  > 0){
          postTotal.significant.i$promoter_5kb[i.set_upDown[j]] <- "yes"
          postTotal.significant.i$gene.ID[i.set_upDown[j]] <- l.df.gff.i[[1]]$gene.ID[i.set[1]]
        }
      }
      
      
      # now also duoble partition assignments 
      i.set <- which(l.df.gff.i[[1]]$end_gene < postTotal.significant.i[i.set_upDown[j],idx.pos] & postTotal.significant.i[i.set_upDown[j],idx.pos] < l.df.gff.i[[1]]$down_1kb
                     | l.df.gff.i[[1]]$down_1kb < postTotal.significant.i[i.set_upDown[j],idx.pos] & postTotal.significant.i[i.set_upDown[j],idx.pos] < l.df.gff.i[[1]]$end_gene)
      
      if(length(i.set)  > 0){
        
        postTotal.significant.i$post_gene_1kb[i.set_upDown[j]] <- "yes"
        postTotal.significant.i$gene.ID[i.set_upDown[j]] <- l.df.gff.i[[1]]$gene.ID[i.set[1]]
        
      }else{
        i.set <- which(l.df.gff.i[[1]]$end_gene < postTotal.significant.i[i.set_upDown[j],idx.pos] & postTotal.significant.i[i.set_upDown[j],idx.pos] < l.df.gff.i[[1]]$down_5kb
                       | l.df.gff.i[[1]]$down_5kb < postTotal.significant.i[i.set_upDown[j],idx.pos] & postTotal.significant.i[i.set_upDown[j],idx.pos] < l.df.gff.i[[1]]$end_gene)
        if(length(i.set) > 0){
          postTotal.significant.i$post_gene_5kb[i.set_upDown[j]] <- "yes"
          postTotal.significant.i$gene.ID[i.set_upDown[j]] <- l.df.gff.i[[1]]$gene.ID[i.set[1]]
        }
      }
    }
    close(pb)
    postTotal.significant.i
  }
  stopCluster(cl)
  print(Sys.time()-strt)
  
  
  df.bQTL_gene_partitioning <- c()
  for(i in 1:n.chromosomes){
    df.bQTL_gene_partitioning <- rbind(df.bQTL_gene_partitioning, l.bQTL_gene_partitioning_per_chromosome[[i]])
  }
  
  
  df.bQTL_gene_partitioning["non_genic"] <- "no"
  i.set <- apply(df.bQTL_gene_partitioning, 1, function(m) {all(m[v.partitions] == "no")})
  df.bQTL_gene_partitioning$non_genic[i.set] <- "yes"
  for(j in 1:length(v.partitions)){
    print(v.partitions[j])
    tmp <- subset(df.bQTL_gene_partitioning, df.bQTL_gene_partitioning[,v.partitions[j]] == "yes")
    print(mean(abs(tmp$POSTfreq - 0.5)))
  }
  
  # numbers - summary
  df.partitions[1,s+1] <- table(df.bQTL_gene_partitioning$promoter_5kb)[2]
  df.partitions[2,s+1] <- table(df.bQTL_gene_partitioning$promoter_1kb)[2]
  df.partitions[3,s+1] <- table(df.bQTL_gene_partitioning$gene)[2]
  df.partitions[4,s+1] <- table(df.bQTL_gene_partitioning$five_prime_UTR)[2]
  df.partitions[5,s+1] <- table(df.bQTL_gene_partitioning$exon)[2]
  df.partitions[6,s+1] <- table(df.bQTL_gene_partitioning$intron)[2]
  df.partitions[7,s+1] <- table(df.bQTL_gene_partitioning$three_prime_UTR)[2]
  df.partitions[8,s+1] <- table(df.bQTL_gene_partitioning$post_gene_1kb)[2]
  df.partitions[9,s+1] <- table(df.bQTL_gene_partitioning$post_gene_5kb)[2]
  
  df.partitions[10,s+1] <- nrow(postTotal.significant) - sum(df.partitions[1:9,s+1])
  df.partitions[11,s+1] <- nrow(postTotal.significant)
  
  # saveRDS(df.partitions, "df.partitions.rds")
  # test <- apply(df.bQTL_gene_partitioning, 1, function(m) {length(which(m[c(17:18,20:24)] == "yes"))})
  
  # i.set <- apply( l.bQTL_gene_partitioning[[2]],1, function(m) {all(m[4:12] == "no")})
  
  print(df.partitions)
  
  length(unique(df.bQTL_gene_partitioning$gene.ID))
  
  # add the AGPv3 gene annotation
  if(!is.null(df.gene_conversion.AGPv3_to_AGPv4)){
    df.bQTL_gene_partitioning["gene.ID.AGPv3"] <- NA
    #idx.geneID.AGPv3 <- which(df.bQTL_gene_partitioning$gene.ID %in% df.gene_conversion.AGPv3_to_AGPv4$gene.ID.AGPv4)
    for(i in 1:nrow(df.bQTL_gene_partitioning)){
      idx.AGPv3 <- which(df.gene_conversion.AGPv3_to_AGPv4$gene.ID.AGPv4 == df.bQTL_gene_partitioning$gene.ID[i])
      if(length(idx.AGPv3) > 0){
        df.bQTL_gene_partitioning$gene.ID.AGPv3[i] <-  df.gene_conversion.AGPv3_to_AGPv4$gene.ID.AGPv3[idx.AGPv3[1]]
      }
    }
  }
  
  
  message("gene function annotation based on AGPv4")
  if(!is.null(df.gene_function)){
    df.bQTL_gene_partitioning["gene.function"] <- NA
    #idx.geneID_to_function <- which(df.bQTL_gene_partitioning$gene.ID %in% df.gene_function$gene.ID)
    for(i in 1:nrow(df.bQTL_gene_partitioning)){
      idx.geneID_to_function <- which(df.gene_function$gene.ID == df.bQTL_gene_partitioning$gene.ID[i])
      if(length(idx.geneID_to_function) > 0){
        df.bQTL_gene_partitioning$gene.function[i] <-  df.gene_function$gene.function[idx.geneID_to_function[1]]
      }
    }
  }
  
  message("add arabidopsis ortholog")
  if(!is.null(df.gene_orthologs)){
    df.bQTL_gene_partitioning["Arabidopsis_ortholog"] <- NA
    for(i in 1:nrow(df.bQTL_gene_partitioning)){
      idx.AGPv3 <- which(df.gene_orthologs$gene.ID == df.bQTL_gene_partitioning$gene.ID[i])
      if(length(idx.AGPv3) > 0){
        df.bQTL_gene_partitioning$Arabidopsis_ortholog[i] <-  df.gene_orthologs$Ara_gene.ID[idx.AGPv3[1]]
      }
    }
    
  }
  
  # TODO: create Table 4
  
  
  write.csv(df.bQTL_gene_partitioning, "output/df.bQTL_gene_partitioning_peaks.csv", row.names = FALSE)
  
  
  # gene function enrichment tests 
  # BARPLOTS mit STERNCHEN (if p < 0.05) - identify biological processes
  
  
  df.snp_genes_with_functions <- unique(df.bQTL_gene_partitioning[,c("gene.ID", "gene.function")])
  df.snp_genes_with_functions <- subset(df.snp_genes_with_functions, !is.na(df.snp_genes_with_functions$gene.ID))
  
  tb.snp_genes_with_functions  <- table(df.genes_with_functions$gene.function)
  v.snp_genes_with_functions  <- unique(df.genes_with_functions$gene.function)
  v.snp_genes_with_functions <- v.snp_genes_with_functions[!is.na(v.snp_genes_with_functions)]
  
  tb.gene_function <- table(df.gene_function$gene.function)
  
  
  # df.GO.annot.set <- subset(df.GO.annot, df.GO.annot$acc. %in% df.bQTL_gene_partitioning_with_gwas$gene.ID)
  # 
  # df.tb.annot <- as.data.frame(table(df.GO.annot$Term[df.GO.annot$Term != ""]), stringsAsFactors = FALSE)
  # df.tb.annot.set <- as.data.frame(table(df.GO.annot.set$Term[df.GO.annot.set$Term != ""]), stringsAsFactors = FALSE)
  # df.tb.annot.set <- merge(df.tb.annot.set, df.tb.annot, by = "Var1")
  # 
  ## filter go - evidence codes, biological process
  df.gene_function_enrichment <- data.frame(gene_function = v.snp_genes_with_functions, 
                                            p.val = rep(1, length(v.snp_genes_with_functions)),
                                            n.genes = numeric(length(v.snp_genes_with_functions)), 
                                            foldchange = numeric(length(v.snp_genes_with_functions)), 
                                            stringsAsFactors = FALSE)
  
  for(j in 1:length(v.snp_genes_with_functions)){  
    
    ### global enrichment test - gene basis
    hitInSample <- tb.snp_genes_with_functions[v.snp_genes_with_functions[j]]
    sampleSize <- sum(tb.snp_genes_with_functions)
    
    hitInPop <- tb.gene_function[v.snp_genes_with_functions[j]] #sum(tb.rate_limiting_domains$Freq)
    popSize <- sum(tb.gene_function)
    
    failInPop <- popSize - hitInPop #(nrow(df.global.domains) - hitInPop)
    p.val <- print(phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE))
    
    #fisher.test(matrix(c(hitInSample-1, hitInPop, failInPop, sampleSize), 2, 2), alternative='less');
    foldChange <- (hitInSample / sampleSize) / (hitInPop / popSize)
    #   mat.count <- matrix(c(n.inset,n.genes.sset - n.inset,  n.genomewide ,n.genes - n.genomewide), ncol = 2, byrow = FALSE)
    #   p.val <- fisher.test(mat.count)$p.value
    df.gene_function_enrichment$p.val[j] <- p.val
    df.gene_function_enrichment$n.genes[j] = hitInSample
    df.gene_function_enrichment$foldchange[j] = foldChange
    
  }
  
  # integrate multiple hypothesis testing
  df.gene_function_enrichment <- df.gene_function_enrichment[order(df.gene_function_enrichment$p.val),]
  df.gene_function_enrichment <- subset(df.BP_enrichment, df.BP_enrichment$p.val <= 0.05)
  
  
  head(df.bQTL_gene_partitioning)

  # TODO: create Table 4
  write.csv(df.bQTL_gene_partitioning, "output/df.bQTL_gene_partitioning_peaks.csv", row.names = FALSE)
  
  df.bQTL_gene_partitioning.subset <- subset(df.bQTL_gene_partitioning, df.bQTL_gene_partitioning$nonGenic == "no")
  length(unique(df.bQTL_gene_partitioning.subset$gene.ID))
  
  length(unique(df.bQTL_gene_partitioning.subset$Arabidopsis_ortholog))
  
  
  df.ChipSeq.gene_partitioning.subset <- subset(df.ChipSeq.gene_partitioning, !is.na(df.ChipSeq.gene_partitioning$Arabidopsis_ortholog))
  df.ChipSeq.gene_partitioning.subset$Arabidopsis_ortholog <- gsub("\\..*","", df.ChipSeq.gene_partitioning.subset$Arabidopsis_ortholog)
  
  
  
}
#write.csv(df.bQTL_gene_partitioning, "output/df.bQTL_gene_partitioning_ChipSeq.csv", row.names = FALSE)

# B73 Chipseq und ASBs with genes

# TODO: remove misleading names .. 

df.ChipSeq.gene_partitioning = read.table("D:/junkDNA.ai/Projects/HASCHSEQ/HaschSeqBasePrototype/output/gene_partitioning/S2.txt", header = T)

#df.ChipSeq.gene_partitioning <- df.bQTL_gene_partitioning # TODO: remove 
comparative_intersection_analysis <- function(df.ChipSeq.gene_partitioning=df.ChipSeq.gene_partitioning){
  
  message("paper only - overlap analysis")
  
  # load all arabidopsis datasets (old)
  
  
  library(VennDiagram)
  library(xlsx)
  
  
  
  df.GO.annot <- readRDS("data/ArabidopsisScriptsAndDatasets/df.GO.annot.Ath.rds")
  
  #df.GO.annot <- readRDS("df.GO_Maize_Trimmed.rds")
  df.GO.annot$acc. <- gsub("\\..*","", df.GO.annot$acc.)
  df.GO.annot <- unique(df.GO.annot)
  
  
  
  # ath chip seq - dark
  df.At_BZR1_ChipSeqTargets <- read.table("data/ArabidopsisScriptsAndDatasets/At_BZR1_ChipSeqTargets.txt", header = TRUE, sep ="\t", quote = "", stringsAsFactors = FALSE)
  names(df.At_BZR1_ChipSeqTargets) <- "locus"
  
  # ath chIP chip - light
  df.At_BZR1_ChIPchipTargets <- read.table("data/ArabidopsisScriptsAndDatasets/At_BZR1_ChIPchip.txt", header = TRUE, sep ="\t", quote = "", stringsAsFactors = FALSE)
  
  df.ChipSeq.gene_partitioning.subset <- subset(df.ChipSeq.gene_partitioning, !is.na(df.ChipSeq.gene_partitioning$Arabidopsis_ortholog))
  df.ChipSeq.gene_partitioning.subset$Arabidopsis_ortholog <- gsub("\\..*","", df.ChipSeq.gene_partitioning.subset$Arabidopsis_ortholog)
  
  ######## 
  
  
  #df.homologies.mapping <- subset(df.homologies.subset, df.homologies.subset$locus %in% df.At_BZR1_ChIPchipTargets$locus)
  #df.homologies.mapping <- subset(df.homologies.subset, df.homologies.subset$locus %in% df.At_BZR1_ChipSeqTargets$locus)
  
  
  #length(unique(intersect(intersect(df.At_BZR1_ChipSeqTargets$locus, df.At_BZR1_ChIPchipTargets$locus), df.ChipSeq.gene_partitioning.subset$Arabidopsis_ortholog)))
  
  a1 <- length(unique(df.At_BZR1_ChIPchipTargets$locus))
  a2 <- length(unique(df.At_BZR1_ChipSeqTargets$locus))
  a3 <- length(unique(df.ChipSeq.gene_partitioning.subset$Arabidopsis_ortholog)) 
  
  a12 <- length(intersect(unique(df.At_BZR1_ChIPchipTargets$locus), unique(df.At_BZR1_ChipSeqTargets$locus)))
  a23 <- length(intersect(unique(df.At_BZR1_ChipSeqTargets$locus), unique(df.ChipSeq.gene_partitioning.subset$Arabidopsis_ortholog)))
  a13 <- length(intersect(unique(df.At_BZR1_ChIPchipTargets$locus), unique(df.ChipSeq.gene_partitioning.subset$Arabidopsis_ortholog)))
  a123 <- length(intersect(intersect(unique(df.At_BZR1_ChIPchipTargets$locus), unique(df.At_BZR1_ChipSeqTargets$locus)), unique(df.ChipSeq.gene_partitioning.subset$Arabidopsis_ortholog)))
  
  
  lst.sets <- vector(mode = "list", length = 7)
  lst.sets[[1]] <- intersect(intersect(unique(df.At_BZR1_ChIPchipTargets$locus), unique(df.At_BZR1_ChipSeqTargets$locus)), unique(df.ChipSeq.gene_partitioning.subset$Arabidopsis_ortholog))
  
  lst.sets[[2]] <- intersect(unique(df.At_BZR1_ChIPchipTargets$locus), unique(df.At_BZR1_ChipSeqTargets$locus)) 
  lst.sets[[3]] <- intersect(unique(df.At_BZR1_ChipSeqTargets$locus), unique(df.ChipSeq.gene_partitioning.subset$Arabidopsis_ortholog))
  lst.sets[[4]] <- intersect(unique(df.At_BZR1_ChIPchipTargets$locus), unique(df.ChipSeq.gene_partitioning.subset$Arabidopsis_ortholog))
  lst.sets[[2]] <- lst.sets[[2]][!lst.sets[[2]] %in% lst.sets[[1]]]
  lst.sets[[3]] <- lst.sets[[3]][!lst.sets[[3]] %in% lst.sets[[1]]]
  lst.sets[[4]] <- lst.sets[[4]][!lst.sets[[4]] %in% lst.sets[[1]]]
  
  lst.sets[[5]] <- unique(df.At_BZR1_ChIPchipTargets$locus)
  lst.sets[[6]] <- unique(df.At_BZR1_ChipSeqTargets$locus)
  lst.sets[[7]] <- unique(df.ChipSeq.gene_partitioning.subset$Arabidopsis_ortholog)
  
  lst.sets[[5]] <- lst.sets[[5]][!lst.sets[[5]] %in% c(unique(df.At_BZR1_ChipSeqTargets$locus), unique(df.ChipSeq.gene_partitioning.subset$Arabidopsis_ortholog))]
  lst.sets[[6]] <- lst.sets[[6]][!lst.sets[[6]] %in% c(unique(df.At_BZR1_ChIPchipTargets$locus), unique(df.ChipSeq.gene_partitioning.subset$Arabidopsis_ortholog))]
  lst.sets[[7]] <- lst.sets[[7]][!lst.sets[[7]] %in% c(unique(df.At_BZR1_ChIPchipTargets$locus), unique(df.At_BZR1_ChipSeqTargets$locus))]
  
  names(lst.sets) <- c("a123_466", "a12_1038","a23_525","a13_420", "a1_1487","a2_2272","a3_2784")
  
  
  # write genes 
  a123_466 <- lst.sets[[1]]
  a12_1038 <- lst.sets[[2]]
  a23_525 <- lst.sets[[3]]
  a13_420 <- lst.sets[[4]]
  a1_1487 <- lst.sets[[5]]
  a2_2272 <- lst.sets[[6]]
  a3_2784 <- lst.sets[[7]]
  save.xlsx("GeneSets.xlsx", a123_466, a12_1038, a23_525, a13_420, a1_1487, a2_2272, a3_2784)
  
  
  ## remark: do venn diagramm
  grid.newpage()
  draw.triple.venn(area1 = a1, area2 = a2, area3 = a3, n12 = a12, n23 = a23, n13 = a13, cex = 1.5,  col = 1, cat.cex = 1.5, SetNames=c( "A", "B","A", "B","A", "B","A"),
                   n123 = a123, category = c("AtBZR1 ChipChip", "AtBZR1 ChipSeq", "ZmBZR1 Seq"), lty = 1, 
                   fill = c("blue", "red", "green"))
  
  
}




