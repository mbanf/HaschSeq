add_nearest_gene <- function(postTotal.significant, 
                           pos_peak, chr_Nr,  
                           df.gene_annotation = NULL,
                           n.cpus = 1){
  
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
    
    for(k in 1:nrow(l.df.gff.i[[1]])){
      if(l.df.gff.i[[1]]$strand[k] == "+"){
        l.df.gff.i[[1]]$up_5kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k] - 5000)
        l.df.gff.i[[1]]$up_1kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k] - 1000)
        l.df.gff.i[[1]]$start_gene[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k])
        l.df.gff.i[[1]]$end_gene[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k])
        l.df.gff.i[[1]]$down_1kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k] + 1000)
      }else if(l.df.gff.i[[1]]$strand[k] == "-"){
        l.df.gff.i[[1]]$up_1kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k] + 1000)
        l.df.gff.i[[1]]$up_5kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k] + 5000)
        l.df.gff.i[[1]]$start_gene[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k])
        l.df.gff.i[[1]]$end_gene[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k])
        l.df.gff.i[[1]]$down_1kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k] - 1000)
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
  
  return(df.bQTL_gene_partitioning)
  
}

