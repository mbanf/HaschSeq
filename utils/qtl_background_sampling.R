create_background_QTLs <- function(df.ASBs_genomic_location,
                                  df.snps_genomic_location, 
                                  multiplyer = 10, 
                                  v.partitions = v.partitions, 
                                  n.chromosomes = 10,
                                  seed = 1234){

  ###### create_background_distribution_with_similar_gene_partitioning ####
  
  # b.duplicateRemoval = FALSE, 
  # duplicate removal in regions
  # v.partitions.duplicates <- v.partitions[c(1,2,8)]
  # 
  # if(b.duplicateRemoval){
  #   # identify duplicate partitions and remove => enter into general pipeline
  #   for(i in 1:2){
  #     i.set <- which(apply(l.bQTL_gene_partitioning[[i]][,v.partitions.duplicates], 1, function(m) length(which(m == "yes"))) == 2)
  #     i.set <- (!1:nrow(l.bQTL_gene_partitioning[[i]]) %in% i.set)
  #     l.bQTL_gene_partitioning[[i]] <- l.bQTL_gene_partitioning[[i]][i.set, ]
  #   }
  # }
  
  # equivalent distribution according to multiplyer 
  set.seed(seed)
  df.bgSNPs <- c()
  
  df.ASBs.B73 <- subset(df.ASBs_genomic_location, df.ASBs_genomic_location$POSTfreq > 0.5)
  df.snps.B73 <- subset(df.snps_genomic_location, df.snps_genomic_location$POSTfreq > 0.5)
  v.chromosomes <- unique(df.ASBs.B73$`B73-chr`)
  
  for(chr in v.chromosomes){
    
    df.ASBs.B73.chr <- subset(df.ASBs.B73, df.ASBs.B73$`B73-chr` == chr)
    df.snps.B73.chr <- subset(df.snps.B73, df.snps.B73$`B73-chr` == chr)
    df.snps.B73.chr <- subset(df.snps.B73.chr, !df.snps.B73.chr$`B73-pos` %in% df.ASBs.B73.chr$`B73-pos`)
    
    for(j in 1:length(v.partitions)){
      if(v.partitions[j] != "gene" & any(df.ASBs.B73.chr[,v.partitions[j]] == "yes")){
        n.samples <- length(which(as.character(df.ASBs.B73.chr[,v.partitions[j]]) == "yes")) * multiplyer
        df.snps.B73.j <- df.snps.B73.chr[which(as.character(df.snps.B73.chr[,v.partitions[j]]) == "yes"),]
        i.set <- sample(nrow(df.snps.B73.j), n.samples, replace = FALSE)
        df.bgSNPs <- rbind(df.bgSNPs, df.snps.B73.j[i.set,])
      }
    }
  
  }
  
  
  df.ASBs.Mo17 <- subset(df.ASBs_genomic_location, df.ASBs_genomic_location$POSTfreq < 0.5)
  df.snps.Mo17 <- subset(df.snps_genomic_location, df.snps_genomic_location$POSTfreq < 0.5)
  v.chromosomes <- unique(df.snps.Mo17$`Mo17-chr`)
  
  for(chr in v.chromosomes){
    
    df.ASBs.Mo17.chr <- subset(df.ASBs.Mo17, df.ASBs.Mo17$`Mo17-chr` == chr)
    df.snps.Mo17.chr <- subset(df.snps.Mo17, df.snps.Mo17$`Mo17-chr` == chr)
    df.snps.Mo17.chr <- subset(df.snps.Mo17.chr, !df.snps.Mo17.chr$`Mo17-pos` %in% df.ASBs.Mo17.chr$`Mo17-pos`)
    
    for(j in 1:length(v.partitions)){
      if(v.partitions[j] != "gene" & any(df.ASBs.Mo17.chr[,v.partitions[j]] == "yes")){
        n.samples <- length(which(as.character(df.ASBs.Mo17.chr[,v.partitions[j]]) == "yes")) * multiplyer
        df.snps.Mo17.j <- df.snps.Mo17.chr[which(as.character(df.snps.Mo17.chr[,v.partitions[j]]) == "yes"),]
        i.set <- sample(nrow(df.snps.Mo17.j), n.samples)
        df.bgSNPs <- rbind(df.bgSNPs, df.snps.Mo17.j[i.set,])
      }
    }
    
  }
  
  
  df.bgSNPs <- unique(df.bgSNPs)
  
  # plot ratio of distribution
  
  v.distribution <- (apply(df.ASBs_genomic_location[,v.partitions], 2, table)["yes",])
  v.distribution.bg_equivalent <- (apply(df.bgSNPs[,v.partitions], 2, table)["yes",])
  df.distribution <- data.frame(ASBs = v.distribution, bgSNPs = v.distribution.bg_equivalent, multiplyer = (v.distribution.bg_equivalent / v.distribution))
  print(df.distribution)
 
  return(list(df.bgSNPs = df.bgSNPs, df.distribution = df.distribution))
  
}

