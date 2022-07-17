add_genomic_location <- function(df.bQTLs, 
                                 chr = "B73-chr",
                                 pos = "B73-pos",
                                 df.gene_annotation,
                                 v.gene_partitions = c("gene", "five_prime_UTR",  "CDS", "three_prime_UTR", "exon"),
                                 n.chromosomes = 10,
                                 n.cpus = 1){
  
  strt<-Sys.time() 
  cl<-makeCluster(min(n.chromosomes, n.cpus))
  registerDoParallel(cl)
  
  l.bQTL_genomic_location_per_chromosome <-  foreach(i = 1:n.chromosomes, .packages=c("seqinr", "VariantAnnotation", "Biostrings")) %dopar% { 
    
    chr.i <- paste(chr, i, sep = "")
    
    df.bQTLs.i <- subset(df.bQTLs, df.bQTLs[,chr] == chr.i)
    
    df.bQTLs.i["promoter_5kb"] <- "no"
    df.bQTLs.i["promoter_1kb"] <- "no"
    df.bQTLs.i["gene"] <- "no"
    df.bQTLs.i["five_prime_UTR"] <- "no"
    df.bQTLs.i["exon"] <- "no"
    df.bQTLs.i["intron"] <- "no"
    df.bQTLs.i["three_prime_UTR"] <- "no"
    df.bQTLs.i["post_gene_1kb"] <- "no"
    
    df.bQTLs.i["non_genic"] <- "no"
    
    df.bQTLs.i["gene.ID"] <- NA
    df.bQTLs.i["transcript.ID"] <- NA
    df.bQTLs.i["strand"] <- NA
    
    df.bQTLs.i["distance_to_gene"] <- NA
    
    
    df.gff.i <- subset(df.gene_annotation, df.gene_annotation$chr == chr.i)
    
    l.df.gff.i <- vector(mode = "list", length = length(v.gene_partitions))
    l.df.gff.i[[1]] <- subset(df.gff.i, df.gff.i$partition %in% v.gene_partitions[1])
    l.df.gff.i[[2]] <- subset(df.gff.i, df.gff.i$partition %in% v.gene_partitions[2])
    l.df.gff.i[[3]] <- subset(df.gff.i, df.gff.i$partition %in% v.gene_partitions[3])
    l.df.gff.i[[4]] <- subset(df.gff.i, df.gff.i$partition %in% v.gene_partitions[4])
    l.df.gff.i[[5]] <- subset(df.gff.i, df.gff.i$partition %in% v.gene_partitions[5])
    
    # non strand sensitive part - in gene
    pb <- txtProgressBar(min = 0, max = nrow(df.bQTLs.i), style = 3)
    for(j in 1:nrow(df.bQTLs.i)){
      setTxtProgressBar(pb, j)

      # in gene
      i.set <- which(l.df.gff.i[[1]]$pos.start < df.bQTLs.i[j,pos] & df.bQTLs.i[j,pos] < l.df.gff.i[[1]]$pos.stop)
      
      if(length(i.set) > 0){
        
        df.bQTLs.i$gene[j] <- "yes"
        df.bQTLs.i$gene.ID[j] <- l.df.gff.i[[1]]$gene.ID[i.set[1]]
        df.bQTLs.i$strand[j] <- l.df.gff.i[[1]]$strand[i.set[1]]
        
        df.bQTLs.i$distance_to_gene[j] <- 0
        
        # in 5 prime utr
        i.set <- which(l.df.gff.i[[2]]$pos.start < df.bQTLs.i[j,pos] & df.bQTLs.i[j,pos] < l.df.gff.i[[2]]$pos.stop)
        if(length(i.set) > 0){
          df.bQTLs.i$five_prime_UTR[j] <- "yes"
        }else{
          # in 3 prime utr
          i.set <- which(l.df.gff.i[[4]]$pos.start < df.bQTLs.i[j,pos] & df.bQTLs.i[j,pos] < l.df.gff.i[[4]]$pos.stop)
          if(length(i.set) > 0){
            df.bQTLs.i$three_prime_UTR[j] <- "yes"
          }else{
            # in CDS (exon)
            i.set <- which(l.df.gff.i[[5]]$pos.start < df.bQTLs.i[j,pos] & df.bQTLs.i[j,pos] < l.df.gff.i[[5]]$pos.stop)
            if(length(i.set) > 0){
              df.bQTLs.i$exon[j] <- "yes"
            }else{
              df.bQTLs.i$intron[j] <- "yes"
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
      if(l.df.gff.i[[1]]$strand[k] == "+"){  # promoter before (-) the gene
        
        l.df.gff.i[[1]]$up_5kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k] - 5000)
        l.df.gff.i[[1]]$up_1kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k] - 1000)
        l.df.gff.i[[1]]$start_gene[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k])
        
        l.df.gff.i[[1]]$end_gene[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k])
        l.df.gff.i[[1]]$down_1kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k] + 1000)
        l.df.gff.i[[1]]$down_5kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k] + 5000)
        
      }else if(l.df.gff.i[[1]]$strand[k] == "-"){  # promoter after (+) the gene - perspective from positive strand 
        
        l.df.gff.i[[1]]$up_1kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k] + 1000)
        l.df.gff.i[[1]]$up_5kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k] + 5000)
        
        l.df.gff.i[[1]]$start_gene[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k])
        
        l.df.gff.i[[1]]$end_gene[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k])
        l.df.gff.i[[1]]$down_1kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k] - 1000)
        l.df.gff.i[[1]]$down_5kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k] - 5000)
        
      }
    }
    
    
    # search only for genes?
    i.set_upDown <- which(df.bQTLs.i$gene == "no")
    
    pb <- txtProgressBar(min = 0, max = length(i.set_upDown), style = 3)
    
    for(j in 1:length(i.set_upDown)){
      
      setTxtProgressBar(pb, j)
      
      # in 1 KB upstream - both strands 
      i.set <- which(l.df.gff.i[[1]]$up_1kb < df.bQTLs.i[i.set_upDown[j], pos] & df.bQTLs.i[i.set_upDown[j], pos] < l.df.gff.i[[1]]$start_gene 
                     | l.df.gff.i[[1]]$start_gene < df.bQTLs.i[i.set_upDown[j], pos] & df.bQTLs.i[i.set_upDown[j], pos] < l.df.gff.i[[1]]$up_1kb)
      
      if(length(i.set) > 0){
        
        if(df.bQTLs.i$gene[i.set_upDown[j]] == "no"){
          df.bQTLs.i$promoter_1kb[i.set_upDown[j]] <- "yes"
        }
        
        df.bQTLs.i$gene.ID[i.set_upDown[j]] <- l.df.gff.i[[1]]$gene.ID[i.set[1]]
        df.bQTLs.i$strand[i.set_upDown[j]] <- l.df.gff.i[[1]]$strand[i.set[1]]
        
        if(df.bQTLs.i$strand[i.set_upDown[j]] == "+"){
          df.bQTLs.i$distance_to_gene[i.set_upDown[j]] <- l.df.gff.i[[1]]$pos.start[i.set[1]] - df.bQTLs.i[i.set_upDown[j], pos]
        }else{
          df.bQTLs.i$distance_to_gene[i.set_upDown[j]] <- df.bQTLs.i[i.set_upDown[j], pos] - l.df.gff.i[[1]]$pos.stop[i.set[1]]
        }
        
      }else{
        
        # in 5 KB upstream - consider both strands
        i.set <- which(l.df.gff.i[[1]]$up_5kb < df.bQTLs.i[i.set_upDown[j], pos] & df.bQTLs.i[i.set_upDown[j],pos] < l.df.gff.i[[1]]$start_gene
                       | l.df.gff.i[[1]]$start_gene < df.bQTLs.i[i.set_upDown[j], pos] & df.bQTLs.i[i.set_upDown[j], pos] < l.df.gff.i[[1]]$up_5kb)
        
        if(length(i.set)  > 0){
          
          if(df.bQTLs.i$gene[i.set_upDown[j]] == "no"){
            df.bQTLs.i$promoter_5kb[i.set_upDown[j]] <- "yes"  
          }
          
          df.bQTLs.i$gene.ID[i.set_upDown[j]] <- l.df.gff.i[[1]]$gene.ID[i.set[1]]
          df.bQTLs.i$strand[i.set_upDown[j]] <- l.df.gff.i[[1]]$strand[i.set[1]]
          
          if(df.bQTLs.i$strand[i.set_upDown[j]] == "+"){
            df.bQTLs.i$distance_to_gene[i.set_upDown[j]] <- l.df.gff.i[[1]]$pos.start[i.set[1]] - df.bQTLs.i[i.set_upDown[j], pos]
          }else{
            df.bQTLs.i$distance_to_gene[i.set_upDown[j]] <- df.bQTLs.i[i.set_upDown[j], pos] - l.df.gff.i[[1]]$pos.stop[i.set[1]]
          }
          
        }
        
      }
      
      
      # now also duoble partition assignments 
      i.set <- which(l.df.gff.i[[1]]$end_gene < df.bQTLs.i[i.set_upDown[j], pos] & df.bQTLs.i[i.set_upDown[j], pos] < l.df.gff.i[[1]]$down_1kb
                     | l.df.gff.i[[1]]$down_1kb < df.bQTLs.i[i.set_upDown[j], pos] & df.bQTLs.i[i.set_upDown[j], pos] < l.df.gff.i[[1]]$end_gene)
      
      if(length(i.set)  > 0){
        
        if(df.bQTLs.i$gene[i.set_upDown[j]] == "no"){
          df.bQTLs.i$post_gene_1kb[i.set_upDown[j]] <- "yes"
        }
        
        df.bQTLs.i$gene.ID[i.set_upDown[j]] <- l.df.gff.i[[1]]$gene.ID[i.set[1]]
        df.bQTLs.i$strand[i.set_upDown[j]] <- l.df.gff.i[[1]]$strand[i.set[1]]
        
        if(df.bQTLs.i$strand[i.set_upDown[j]] == "-"){ # inverse because post gene
          df.bQTLs.i$distance_to_gene[i.set_upDown[j]] <- l.df.gff.i[[1]]$pos.start[i.set[1]] - df.bQTLs.i[i.set_upDown[j], pos]
        }else{
          df.bQTLs.i$distance_to_gene[i.set_upDown[j]] <- df.bQTLs.i[i.set_upDown[j], pos] - l.df.gff.i[[1]]$pos.stop[i.set[1]]
        }
      }
      
    }
    close(pb)
    
    df.bQTLs.i
    
  }
  
  stopCluster(cl)
  print(Sys.time()-strt)
  
  
  df.bQTL_genomic_location <- c()
  for(i in 1:n.chromosomes){
    df.bQTL_genomic_location <- rbind(df.bQTL_genomic_location, l.bQTL_genomic_location_per_chromosome[[i]])
  }
  
  # if not found anywhere near or in gene
  i.set <- apply(df.bQTL_genomic_location, 1, function(m) {all(m[v.partitions] == "no")})
  i.set <- which(i.set == TRUE)
  df.bQTL_genomic_location$non_genic[i.set] <- "yes"
  
  return(df.bQTL_genomic_location)
}



add_genomic_location_bQTLs <- function(df.bQTLs, 
                                       df.gene_annotation,
                                       v.gene_partitions = c("gene", "five_prime_UTR",  "CDS", "three_prime_UTR", "exon"),
                                       n.chromosomes = 10,
                                       n.cpus = 1){
  # added >= 0.5 postfrequency for genomic location to include 0.5 allelic bias SNPs (will not affect ASBs anyways)
  df.bQTL_genomic_location <- add_genomic_location(subset(df.bQTLs, df.bQTLs$POSTfreq >= 0.5), 
                                                   chr = "B73-chr",
                                                   pos = "B73-pos",
                                                   df.gene_annotation,
                                                   v.gene_partitions = v.gene_partitions,
                                                   n.chromosomes = n.chromosomes,
                                                   n.cpus = n.cpus)
  
  
  df.bQTL_genomic_location <- rbind(df.bQTL_genomic_location, 
                                    add_genomic_location(subset(df.bQTLs, df.bQTLs$POSTfreq < 0.5), 
                                                         chr = "Mo17-chr",
                                                         pos = "Mo17-pos",
                                                         df.gene_annotation,
                                                         v.gene_partitions = v.gene_partitions,
                                                         n.chromosomes = n.chromosomes,
                                                         n.cpus = n.cpus))
        
  # TODO: plot ratio of annotation !! - table
  
  return(df.bQTL_genomic_location)
}



get_partitions <- function(df.bQTL_genomic_location,
                           v.genePartitions = c("gene", "five_prime_UTR",  "CDS", "three_prime_UTR", "exon")){
  
  for(j in 1:length(v.partitions)){
    print(v.partitions[j])
    tmp <- subset(df.bQTL_genomic_location, df.bQTL_genomic_location[,v.partitions[j]] == "yes")
    print(mean(abs(tmp$POSTfreq - 0.5)))
  }
  
  # numbers - summary
  df.partitions[1,s+1] <- table(df.bQTL_genomic_location$promoter_5kb)[2]
  df.partitions[2,s+1] <- table(df.bQTL_genomic_location$promoter_1kb)[2]
  df.partitions[3,s+1] <- table(df.bQTL_genomic_location$gene)[2]
  df.partitions[4,s+1] <- table(df.bQTL_genomic_location$five_prime_UTR)[2]
  df.partitions[5,s+1] <- table(df.bQTL_genomic_location$exon)[2]
  df.partitions[6,s+1] <- table(df.bQTL_genomic_location$intron)[2]
  df.partitions[7,s+1] <- table(df.bQTL_genomic_location$three_prime_UTR)[2]
  df.partitions[8,s+1] <- table(df.bQTL_genomic_location$post_gene_1kb)[2]
  # df.partitions[9,s+1] <- table(df.bQTL_gene_partitioning$post_gene_5kb)[2]
  
  # adapted - no post 5 kb
  df.partitions[9,s+1] <- nrow(df.bQTL_genomic_location) - sum(df.partitions[1:8,s+1])
  df.partitions[10,s+1] <- nrow(df.bQTL_genomic_location)
  
  df.partitions["percentage_significant"] <- df.partitions[,2] / df.partitions[10,2]
  df.partitions["percentage_non_significant"] <- df.partitions[,3] / df.partitions[10,3]
  
  return(df.partitions)
}


