

add_genomic_location <- function(postTotal.significant, 
                                 pos_peak,
                                 chr_Nr,
                                 df.gene_annotation,
                                 v.genePartitions = c("gene", "five_prime_UTR",  "CDS", "three_prime_UTR", "exon"),
                                 n.cpus = 1){
  
  idx.pos <- which( colnames(postTotal.significant)== pos_peak)
  idx.chr <- which( colnames(postTotal.significant) == chr_Nr)
  
  df.bQTL_gene_partitioning <- c()
  
  strt<-Sys.time() 
  cl<-makeCluster(min(n.chromosomes, n.cpus))
  registerDoParallel(cl)
  
  l.bQTL_gene_partitioning_per_chromosome <-  foreach(i = 1:n.chromosomes, .packages=c("seqinr", "VariantAnnotation", "Biostrings")) %dopar% { 
    
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
    
    postTotal.significant.i["non_genic"] <- "no"
    
    postTotal.significant.i["gene.ID"] <- NA
    postTotal.significant.i["transcript.ID"] <- NA
    postTotal.significant.i["strand"] <- NA
    
    postTotal.significant.i["distance_to_gene"] <- NA
    
    
    df.gff.i <- subset(df.gene_annotation, df.gene_annotation$chr == i)
    
    l.df.gff.i <- vector(mode = "list", length = length(v.genePartitions))
    l.df.gff.i[[1]] <- subset(df.gff.i, df.gff.i$partition %in% v.genePartitions[1])
    l.df.gff.i[[2]] <- subset(df.gff.i, df.gff.i$partition %in% v.genePartitions[2])
    l.df.gff.i[[3]] <- subset(df.gff.i, df.gff.i$partition %in% v.genePartitions[3])
    l.df.gff.i[[4]] <- subset(df.gff.i, df.gff.i$partition %in% v.genePartitions[4])
    l.df.gff.i[[5]] <- subset(df.gff.i, df.gff.i$partition %in% v.genePartitions[5])
    
    # nearest gene distance 
    
    # preselect all in gene 
    pb <- txtProgressBar(min = 0, max = nrow(postTotal.significant.i), style = 3)
    for(j in 1:nrow(postTotal.significant.i)){
      
      setTxtProgressBar(pb, j)
      
      # non strand sensitive part #
      
      # in gene
      i.set <- which(l.df.gff.i[[1]]$pos.start < postTotal.significant.i[j,idx.pos] & postTotal.significant.i[j,idx.pos] < l.df.gff.i[[1]]$pos.stop)
      
      if(length(i.set) > 0){
        
        postTotal.significant.i$gene[j] <- "yes"
        postTotal.significant.i$gene.ID[j] <- l.df.gff.i[[1]]$gene.ID[i.set[1]]
        postTotal.significant.i$strand[j] <- l.df.gff.i[[1]]$strand[i.set[1]]
        
        
        postTotal.significant.i$distance_to_gene[j] <- 0
        
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
    
    # strand sensitive part #
    
    # identify bQTL in 5kb and 1kb distance 
    l.df.gff.i[[1]]["up_5kb"] <- 0
    l.df.gff.i[[1]]["up_1kb"] <- 0
    l.df.gff.i[[1]]["start_gene"] <- 0
    l.df.gff.i[[1]]["end_gene"] <- 0
    l.df.gff.i[[1]]["down_1kb"] <- 0
    l.df.gff.i[[1]]["down_5kb"] <- 0
    
    for(k in 1:nrow(l.df.gff.i[[1]])){
      
      # strand sensitive 
      if(l.df.gff.i[[1]]$strand[k] == "+"){# promoter before (-) the gene
        
        l.df.gff.i[[1]]$up_5kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k] - 5000)
        l.df.gff.i[[1]]$up_1kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k] - 1000)
        l.df.gff.i[[1]]$start_gene[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k])
        
        l.df.gff.i[[1]]$end_gene[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k])
        l.df.gff.i[[1]]$down_1kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k] + 1000)
        l.df.gff.i[[1]]$down_5kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k] + 5000)
        
      }else if(l.df.gff.i[[1]]$strand[k] == "-"){# promoter after (+) the gene - perspective from positive strand 
        
        l.df.gff.i[[1]]$up_1kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k] + 1000)
        l.df.gff.i[[1]]$up_5kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k] + 5000)
        
        l.df.gff.i[[1]]$start_gene[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k])
        
        l.df.gff.i[[1]]$end_gene[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k])
        l.df.gff.i[[1]]$down_1kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k] - 1000)
        l.df.gff.i[[1]]$down_5kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k] - 5000)
        
        
      }
    }
    
    
    # search only for genes?
    i.set_upDown <- which(postTotal.significant.i$gene == "no")
    
    pb <- txtProgressBar(min = 0, max = length(i.set_upDown), style = 3)
    
    for(j in 1:length(i.set_upDown)){
      
      setTxtProgressBar(pb, j)
      
      # in 1 KB upstream - both strands 
      i.set <- which(l.df.gff.i[[1]]$up_1kb < postTotal.significant.i[i.set_upDown[j], idx.pos] & postTotal.significant.i[i.set_upDown[j], idx.pos] < l.df.gff.i[[1]]$start_gene 
                     | l.df.gff.i[[1]]$start_gene < postTotal.significant.i[i.set_upDown[j], idx.pos] & postTotal.significant.i[i.set_upDown[j], idx.pos] < l.df.gff.i[[1]]$up_1kb)
      
      if(length(i.set) > 0){
        
        if(postTotal.significant.i$gene[i.set_upDown[j]] == "no"){
          postTotal.significant.i$promoter_1kb[i.set_upDown[j]] <- "yes"
        }
        
        postTotal.significant.i$gene.ID[i.set_upDown[j]] <- l.df.gff.i[[1]]$gene.ID[i.set[1]]
        postTotal.significant.i$strand[i.set_upDown[j]] <- l.df.gff.i[[1]]$strand[i.set[1]]
        
        if(postTotal.significant.i$strand[i.set_upDown[j]] == "+"){
          postTotal.significant.i$distance_to_gene[i.set_upDown[j]] <- l.df.gff.i[[1]]$pos.start[i.set[1]] - postTotal.significant.i[i.set_upDown[j], idx.pos]
        }else{
          postTotal.significant.i$distance_to_gene[i.set_upDown[j]] <- postTotal.significant.i[i.set_upDown[j], idx.pos] - l.df.gff.i[[1]]$pos.stop[i.set[1]]
        }
        
        
        # take the closest
        # postTotal.significant.i$distance_to_gene[i.set_upDown[j]] <- 
        
      }else{
        
        # in 5 KB upstream - consider both strands
        i.set <- which(l.df.gff.i[[1]]$up_5kb < postTotal.significant.i[i.set_upDown[j], idx.pos] & postTotal.significant.i[i.set_upDown[j],idx.pos] < l.df.gff.i[[1]]$start_gene
                       | l.df.gff.i[[1]]$start_gene < postTotal.significant.i[i.set_upDown[j], idx.pos] & postTotal.significant.i[i.set_upDown[j], idx.pos] < l.df.gff.i[[1]]$up_5kb)
        
        if(length(i.set)  > 0){
          
          if(postTotal.significant.i$gene[i.set_upDown[j]] == "no"){
            postTotal.significant.i$promoter_5kb[i.set_upDown[j]] <- "yes"  
          }
          
          postTotal.significant.i$gene.ID[i.set_upDown[j]] <- l.df.gff.i[[1]]$gene.ID[i.set[1]]
          postTotal.significant.i$strand[i.set_upDown[j]] <- l.df.gff.i[[1]]$strand[i.set[1]]
          
          if(postTotal.significant.i$strand[i.set_upDown[j]] == "+"){
            postTotal.significant.i$distance_to_gene[i.set_upDown[j]] <- l.df.gff.i[[1]]$pos.start[i.set[1]] - postTotal.significant.i[i.set_upDown[j], idx.pos]
          }else{
            postTotal.significant.i$distance_to_gene[i.set_upDown[j]] <- postTotal.significant.i[i.set_upDown[j], idx.pos] - l.df.gff.i[[1]]$pos.stop[i.set[1]]
          }
          
        }
        
      }
      
      
      # now also duoble partition assignments 
      i.set <- which(l.df.gff.i[[1]]$end_gene < postTotal.significant.i[i.set_upDown[j], idx.pos] & postTotal.significant.i[i.set_upDown[j], idx.pos] < l.df.gff.i[[1]]$down_1kb
                     | l.df.gff.i[[1]]$down_1kb < postTotal.significant.i[i.set_upDown[j], idx.pos] & postTotal.significant.i[i.set_upDown[j], idx.pos] < l.df.gff.i[[1]]$end_gene)
      
      if(length(i.set)  > 0){
        
        if(postTotal.significant.i$gene[i.set_upDown[j]] == "no"){
          postTotal.significant.i$post_gene_1kb[i.set_upDown[j]] <- "yes"
        }
        
        postTotal.significant.i$gene.ID[i.set_upDown[j]] <- l.df.gff.i[[1]]$gene.ID[i.set[1]]
        postTotal.significant.i$strand[i.set_upDown[j]] <- l.df.gff.i[[1]]$strand[i.set[1]]
        
        if(postTotal.significant.i$strand[i.set_upDown[j]] == "-"){ # inverse because post gene
          postTotal.significant.i$distance_to_gene[i.set_upDown[j]] <- l.df.gff.i[[1]]$pos.start[i.set[1]] - postTotal.significant.i[i.set_upDown[j], idx.pos]
        }else{
          postTotal.significant.i$distance_to_gene[i.set_upDown[j]] <- postTotal.significant.i[i.set_upDown[j], idx.pos] - l.df.gff.i[[1]]$pos.stop[i.set[1]]
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
  
  
  # if not found anywehere near or in gene
  i.set <- apply(df.bQTL_gene_partitioning, 1, function(m) {all(m[v.partitions] == "no")})
  i.set <- which(i.set == TRUE)
  df.bQTL_gene_partitioning$non_genic[i.set] <- "yes"
  
  return(df.bQTL_gene_partitioning)
  
}



get_partitions <- function(df.bQTL_gene_partitioning,
                           v.genePartitions = c("gene", "five_prime_UTR",  "CDS", "three_prime_UTR", "exon")){
  
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
  # df.partitions[9,s+1] <- table(df.bQTL_gene_partitioning$post_gene_5kb)[2]
  
  # adapted - no post 5 kb
  df.partitions[9,s+1] <- nrow(df.bQTL_gene_partitioning) - sum(df.partitions[1:8,s+1])
  df.partitions[10,s+1] <- nrow(df.bQTL_gene_partitioning)
  
  df.partitions["percentage_significant"] <- df.partitions[,2] / df.partitions[10,2]
  df.partitions["percentage_non_significant"] <- df.partitions[,3] / df.partitions[10,3]
  
  return(df.partitions)
}


