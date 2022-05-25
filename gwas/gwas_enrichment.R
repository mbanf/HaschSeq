
gwas_enrichment <- function(l.bQTL_gene_partitioning, v.traits, output_file, s.dist_ASB_to_GWAS=2000){
  
  message("Phenotypic GWAS enrichment")
    
  l.bQTL_gene_partitioning_with_gwas <- vector(mode = "list", length = 2)
  l.number_per_trait <- vector(mode = "list", length = 2)
  l.number_per_trait[[1]] <- l.number_per_trait[[2]] <- numeric(length(v.traits))
  names(l.number_per_trait[[1]]) <- names(l.number_per_trait[[2]]) <- v.traits
  
  for(s in 1:2){
    
    if(s == 1){
      message("ASBs")
    }else{
      message("population")
    }
    
    df.bQTL_gene_partitioning <- l.bQTL_gene_partitioning[[s]]
    
    # add trait annotation to ASBs
    if(s == 1)
      df.bQTL_gene_partitioning[v.traits] <- "no"
    
    df.bQTL_gene_partitioning["dist_ASB_to_gwas"] <- NA
    
    
    # dynamic extension for traits
    df.bQTL_gene_partitioning_with_gwas <- c() 
    
    message("number of phenotype gwas: ", nrow(df.phenotype_gwas))
    
    pb <- txtProgressBar(min = 0, max = n.chromosomes, style = 3)
    
    for(i in 1:n.chromosomes){
      
      setTxtProgressBar(pb, i)
      
      df.phenotype_gwas.i <- subset(df.phenotype_gwas, df.phenotype_gwas$chr == i)  
      df.bQTL_gene_partitioning.i <- subset(df.bQTL_gene_partitioning, df.bQTL_gene_partitioning$contig == i)
      
      for(k in 1:length(v.traits)){ # trait
        
        df.phenotype_gwas.i.k <- subset(df.phenotype_gwas.i, df.phenotype_gwas.i$trait == v.traits[k])
        
        if(nrow(df.phenotype_gwas.i.k) > 0){
          
          for(j in 1:nrow(df.bQTL_gene_partitioning.i)){
            
            v.dist <- abs(df.bQTL_gene_partitioning.i$position[j] - df.phenotype_gwas.i.k$pos_start)
            j.set <- which(v.dist <= s.dist_ASB_to_GWAS)
            
            if(length(j.set) > 0){
              
              l.number_per_trait[[s]][unique(df.phenotype_gwas.i.k$trait[j.set])] <- l.number_per_trait[[s]][unique(df.phenotype_gwas.i.k$trait[j.set])] + length(j.set)
              if(s == 1){
                df.bQTL_gene_partitioning.i[j, unique(df.phenotype_gwas.i.k$trait[j.set])] <- "yes"
              }
              df.bQTL_gene_partitioning.i$dist_ASB_to_gwas[j] <- mean(v.dist[j.set])   
            }
          }
        }
      }
      df.bQTL_gene_partitioning_with_gwas <- rbind(df.bQTL_gene_partitioning_with_gwas, df.bQTL_gene_partitioning.i)
    }
    l.bQTL_gene_partitioning_with_gwas[[s]] <- df.bQTL_gene_partitioning_with_gwas
  }
  
  saveRDS(l.bQTL_gene_partitioning_with_gwas, "../tmp/l.bQTL_gene_partitioning_with_gwas.rds")
  saveRDS(l.number_per_trait, "../tmp/l.number_per_trait.rds")
  
  v.distribution <- apply(l.bQTL_gene_partitioning[[1]][,v.partitions], 2, table)["yes",]
  n.sampleSize <- sum(v.distribution) - v.distribution["gene"]
  
  df.enrichment <- data.frame(id = character(length(v.traits)), p.value = numeric(length(v.traits)), foldchange = numeric(length(v.traits)), stringsAsFactors = FALSE)
  
  for(k in 1:length(v.traits)){
    
    hitInSample <- l.number_per_trait[[1]][v.traits[k]]
    sampleSize <- n.sampleSize
    hitInPop <- l.number_per_trait[[2]][v.traits[k]]
    popSize <- nrow(l.bQTL_gene_partitioning[[2]])
    
    failInPop <- popSize - hitInPop
    
    foldchange <- (hitInSample/sampleSize)/(hitInPop/popSize)
    
    df.enrichment$id[k] <- v.traits[k]
    df.enrichment$p.value[k] <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
    df.enrichment$foldchange[k] <-  as.numeric(foldchange)
    
  }
  
  write.table(df.enrichment, output_file,  row.names = FALSE, quote = FALSE, sep ="\t")
  # write.table(df.enrichment, paste(folder_output, "/gwas/S_1000.txt", sep = ""),  row.names = FALSE, quote = FALSE, sep ="\t")
}