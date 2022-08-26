
gwas_enrichment <- function(l.bQTL_gene_partitioning, df.phenotype_gwas, output_file, s.dist_ASB_to_GWAS=2000){
  
  message("Phenotypic GWAS enrichment")
    
  v.traits <- unique(df.phenotype_gwas$trait)
  
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





message("GWAS analysis - comparing distances of ASBs and bgSNPs to nearest GWAS")

snp_distance_to_nearest_gwas <- function(l.bQTL_gene_partitioning, df.phenotype_gwas){
  
  l.distances <- vector(mode = "list", length = 2)
  
  for(s in 1:length(l.bQTL_gene_partitioning)){
    
    df.bQTL_gene_partitioning <- l.bQTL_gene_partitioning[[s]]
    
    # analyze non gene snps
    df.bQTL_gene_partitioning <- subset(df.bQTL_gene_partitioning, df.bQTL_gene_partitioning$gene == "no") 
    
    #strt<-Sys.time() 
    #cl<-makeCluster(min(n.chromosomes, n.cpus))
    #registerDoParallel(cl)
    
    v.dist_bQTL_to_nearest_GWAS <- c()
    
    #l.snp_to_gene_distances_per_chromosome <-  foreach(i = 1:n.chromosomes, .packages=c("seqinr", "VariantAnnotation", "Biostrings")) %dopar% { 
    for(i in 1:n.chromosomes){
      
      df.bQTL_gene_partitioning.i <- subset(df.bQTL_gene_partitioning, df.bQTL_gene_partitioning$contig == i)  
      df.phenotype_gwas.i <- subset(df.phenotype_gwas, df.phenotype_gwas$chr == i)
      
      v.dist_bQTL_to_nearest_GWAS.i <- numeric(nrow(df.bQTL_gene_partitioning.i))
      
      pb <- txtProgressBar(min = 0, max = nrow(df.bQTL_gene_partitioning.i), style = 3)
      
      # find closest gene to snp
      for(j in 1:nrow(df.bQTL_gene_partitioning.i)){
        
        setTxtProgressBar(pb, j)
        
        v.dist_bQTL_to_nearest_GWAS.i[j] <- min(abs(df.phenotype_gwas.i$pos_start - df.bQTL_gene_partitioning.i$position[j]))
        
      }
      
      close(pb)
      
      v.dist_bQTL_to_nearest_GWAS <- c(v.dist_bQTL_to_nearest_GWAS, v.dist_bQTL_to_nearest_GWAS.i)
      
    }
    
    l.distances[[s]] <- v.dist_bQTL_to_nearest_GWAS
    
  }
  
  
  message("plot log2 distances")
  
  norm <- as.data.frame(cbind(l.distances[[2]], c(l.distances[[1]], rep(NA, length(l.distances[[2]]) - length(l.distances[[1]])))))
  names(norm) <- c("bgSNPs", "ASB")
  
  v.bgSNPs <- norm$bgSNPs
  v.ASB <- norm$ASB
  
  norm$bgSNPs = norm$bgSNPs / sum(norm$bgSNPs, na.rm=T)
  norm$ASB = norm$ASB / sum(norm$ASB, na.rm=T)
  
  
  #Divide into breaks
  breaks=2^(0:20)	#Divide into bins by powers of 2 (upper limit determined empirically from data)
  v.bgSNPs.cuts= cut(v.bgSNPs, breaks=breaks, labels=F, include.lowest=T)
  v.ASB.cuts= cut(v.ASB, breaks=breaks, labels=F, include.lowest=T)
  
  
  final=data.frame(bgSNPs=rep(0, length(breaks)), ASB=rep(0, length(breaks)))
  for(group in 1:length(breaks)){
    final$bgSNPs[group]=sum(norm$bgSNPs[v.bgSNPs.cuts==group], na.rm=T)
    final$ASB[group]=sum(norm$ASB[v.ASB.cuts==group], na.rm=T)
  }
  
  ylim=range(final)
  xlim=range(log2(breaks))
  lwd=2
  plot(x=1,y=1, type="n", xlim=xlim, ylim=ylim, main="Non-genic SNP distances to nearest GWAS", xlab="Log2 distance (bp)", ylab="Proportion of non-genic SNPs")
  lines(x=log2(breaks), y=final$bgSNPs, col="black", lwd=lwd)
  lines(x=log2(breaks), y=final$ASB, col="red", lwd=lwd)
  # axis(1, xaxp=c(10, 200, 19), las=2)
  legend(x="topleft", legend=c("background SNPs","GWAS hits"), col=c("black","red"), lwd=2)
  #dev.off()
  
  message("comparative distribution test")
  wilcox.test(v.ASB, v.bgSNPs[sample(length(v.bgSNPs), 10000)])
  
  
  message("ASB vs background SNPs")
  
  norm <- as.data.frame(cbind(l.distances[[2]], c(l.distances[[1]], rep(NA, length(l.distances[[2]]) - length(l.distances[[1]])))))
  names(norm) <- c("bgSNPs", "ASB")
  
  v.bgSNPs <- norm$bgSNPs
  v.ASB <- norm$ASB
  
  
  v.ASB <- l.distances[[1]]
  names(v.ASB) <- l.positions[[1]]
  
  v.ASB <- v.ASB[v.ASB > 0]
  
  v.bgSNPs <- l.distances[[2]]
  names(v.bgSNPs) <- l.positions[[2]]
  v.bgSNPs <- v.bgSNPs[v.bgSNPs > 0]
  
  
  v.ASB <- v.ASB[names(v.ASB) == "-"]
  v.bgSNPs <- v.bgSNPs[names(v.bgSNPs) == "-"]
  
  hist(v.ASB, breaks = 5000, xlim = c(0,15000))
  hist(v.bgSNPs, breaks = 10000, xlim = c(0,15000))
  
}


message("Phenotypic GWAS enrichment ")

# v.traits <- unique(df.geneSet.gwas.set$trait)
gwas_enrichment <- function(){
  
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
              # df.bQTL_gene_partitioning.i$gwas_trait[j] <- paste(df.bQTL_gene_partitioning.i$gwas_trait[j], paste(df.phenotype_gwas.i.k$trait[j.set], collapse = ", "), collapse = "; ")
              # df.bQTL_gene_partitioning.i$dist_bQTL_to_GWAS_snp[j] <- paste(df.bQTL_gene_partitioning.i$dist_bQTL_to_GWAS_snp[j], paste(abs(df.bQTL_gene_partitioning.i$position[j] - df.phenotype_gwas.i.k$pos[j.set]), collapse = ", "), collapse = "; ")  
              
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
    #df.bQTL_gene_partitioning_with_gwas <- subset(df.bQTL_gene_partitioning_with_gwas, df.bQTL_gene_partitioning_with_gwas$gwas_trait != "")
    l.bQTL_gene_partitioning_with_gwas[[s]] <- df.bQTL_gene_partitioning_with_gwas
  }
  
  saveRDS(l.bQTL_gene_partitioning_with_gwas, "tmp/l.bQTL_gene_partitioning_with_gwas.rds")
  saveRDS(l.number_per_trait, "tmp/l.number_per_trait.rds")
  
  
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
    
    #tab <- matrix(c(hitInSample, hitInPop, sampleSize - hitInSample, popSize - hitInPop), nrow = 2)
    #print(paste("equal size comparison: FC:", foldchange," , p-value", fisher.test(tab)$p.value))
    
    df.enrichment$id[k] <- v.traits[k]
    #df.domain_enrichment$p.value[k] <- fisher.test(tab)$p.value
    df.enrichment$p.value[k] <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
    df.enrichment$foldchange[k] <-  as.numeric(foldchange)
    
  }
  
  # subset(df.enrichment, df.enrichment$p.value < 0.05)
  
  #df.enrichment <- subset(df.enrichment, df.enrichment$p.value < 0.1)
  #df.enrichment
  
  # write.csv(df.enrichment, paste("C:/Users/Michael/Desktop/HaschSeq_data/F7b.csv", sep = ""))  
  write.table(df.enrichment, paste(folder_output, "/gwas/S0_1.txt", sep = ""),  row.names = FALSE, quote = FALSE, sep ="\t")
  
}

# log2


message("average sizes - per transcript and gene separately") 

## todo: comparative plot divided over average sizes 


if(!b.loadGenePartitioning){
  
  v.geneIDs <- unique(df.gene_annotation$gene.ID)
  # v.geneIDs <- v.geneIDs[sample(length(v.geneIDs), 10000)] # sample representative number of genes etc.
  
  l.partitionSets <- vector(mode = "list", length = 5)
  for(i in 1:5)
    l.partitionSets[[i]] <- numeric(length(v.geneIDs))
  
  pb <- txtProgressBar(min = 0, max = length(v.geneIDs), style = 3)
  
  for(i in 1:length(v.geneIDs)){
    
    setTxtProgressBar(pb, i)
    
    df.gene_annotation.i <- subset(df.gene_annotation, df.gene_annotation$gene.ID == v.geneIDs[i])
    
    if(df.gene_annotation.i$partition[1] == "gene"){
      
      s.gene_length <- abs(df.gene_annotation.i$pos.start[1] - df.gene_annotation.i$pos.stop[1])
      
      l.partitionSets[[1]][i] <- s.gene_length
      
      # df.average_partition_sizes[1,2] <- df.average_partition_sizes[1,2] + s.gene_length
      # df.average_partition_sizes[1,3] <- df.average_partition_sizes[1,3] + 1
      
    }else{ # per transcript 
      
      i.mRNA <- which(df.gene_annotation.i$partition == "mRNA")
      i.5utr <- which(df.gene_annotation.i$partition == "five_prime_UTR")
      i.3utr <- which(df.gene_annotation.i$partition == "three_prime_UTR")
      i.cds <- which(df.gene_annotation.i$partition == "CDS")
      
      if(length(i.5utr) > 0 & length(i.3utr) > 0 & length(i.cds) > 0 & length(i.mRNA) > 0){
        
        s.mRNA_length = abs(df.gene_annotation.i$pos.start[i.mRNA] - df.gene_annotation.i$pos.stop[i.mRNA])
        
        s.5utr_length = abs(df.gene_annotation.i$pos.start[i.5utr] - df.gene_annotation.i$pos.stop[i.5utr])
        
        s.cds_length <- 0
        
        for(j in 1:length(i.cds)){
          s.cds_length = s.cds_length + abs(df.gene_annotation.i$pos.start[i.cds[j]] - df.gene_annotation.i$pos.stop[i.cds[j]])
        }
        
        s.3utr_length = abs(df.gene_annotation.i$pos.start[i.3utr] - df.gene_annotation.i$pos.stop[i.3utr])
        
        s.intron_length = s.mRNA_length - s.5utr_length - s.cds_length - s.3utr_length
        
        if(FALSE){
          
          s.cds_length <- s.cds_length / length(i.cds)
          
          div.factor <- length(i.cds)
          
          if(abs(df.gene_annotation.i$pos.stop[i.5utr] - df.gene_annotation.i$pos.start[i.cds[1]]) <= 1)
            div.factor <- div.factor + 1
          
          if(abs(df.gene_annotation.i$pos.start[i.3utr] - df.gene_annotation.i$pos.stop[i.cds[length(i.cds)]]) <= 1)
            div.factor <- div.factor + 1
          
          s.intron_length <- s.intron_length / div.factor
          
        }
        #s.intron_length = abs(df.gene_annotation.i$pos.start[i.5utr] - df.gene_annotation.i$pos.start[i.cds]) + abs(df.gene_annotation.i$pos.stop[i.3utr] - df.gene_annotation.i$pos.stop[i.cds])
        
        l.partitionSets[[2]][i] <- s.5utr_length
        l.partitionSets[[3]][i] <- s.cds_length
        l.partitionSets[[4]][i] <- s.intron_length
        l.partitionSets[[5]][i] <- s.3utr_length
        
      }
      
    }
    
  }
  
  close(pb)
  
  saveRDS(l.partitionSets, paste("tmp/l.partitionSets_", timeStamp, ".rds", sep = ""))
}else{
  l.partitionSets <- readRDS(paste("tmp/l.partitionSets_", timeStamp, ".rds", sep = ""))
}


message("Average sizes")

df.average_partition_sizes <- data.frame(partition = c("gene", "five_prime_UTR", "exon", "intron", "three_prime_UTR", "non_genic"), 
                                         mean = numeric(6), 
                                         standard_deviation = numeric(6))

for(i in 1:5){
  df.average_partition_sizes[i, 2] <- (mean(l.partitionSets[[i]][l.partitionSets[[i]] > 0]))
  df.average_partition_sizes[i, 3] <- (sd(l.partitionSets[[i]][l.partitionSets[[i]] > 0]))
}

df.gene_annotation.genes_only <- subset(df.gene_annotation, df.gene_annotation$partition == "gene")
v.nongenic_length <- numeric()

for(i in 1:n.chromosomes){
  df.gene_annotation.genes_only.i <- subset(df.gene_annotation.genes_only, df.gene_annotation.genes_only$chr == i)
  df.gene_annotation.genes_only.i <- df.gene_annotation.genes_only.i[order(-df.gene_annotation.genes_only.i$pos.start),]
  
  v.nongenic_length.i <- abs(df.gene_annotation.genes_only.i$pos.stop[2:nrow(df.gene_annotation.genes_only.i)] - df.gene_annotation.genes_only.i$pos.start[1:(nrow(df.gene_annotation.genes_only.i) - 1)])
  v.nongenic_length.i <- v.nongenic_length.i - 10000 # 5kb upstream and downstream
  
  v.nongenic_length <- c(v.nongenic_length, v.nongenic_length.i)
}

df.average_partition_sizes[6, 2] <- mean(v.nongenic_length)
df.average_partition_sizes[6, 3] <- sd(v.nongenic_length)


# write.csv(df.average_partition_sizes, paste("tmp/df.average_partition_sizes_", timeStamp, ".rds", sep = ""))
# df.average_partition_sizes <- read.csv(paste("tmp/df.average_partition_sizes_", timeStamp, ".rds", sep = ""))

###



#### 








df.partitions


df.average_partition_sizes[,2]/ df.average_partition_sizes[,3]

v.partitions

#   #write.csv(df.postTotal_both.selection, paste("/Shared/Everyone/Michael_Thomas/Parental_effects/df.bQTL_Peak_parental_genes.",as.character(Sys.time()), ".csv", sep = ""))
#   saveRDS(l.bQTL_gene_partitioning, "l.bQTL_gene_partitioning.rds")
#   
#   
#   write.csv(l.bQTL_gene_partitioning[[1]], paste("/Shared/Everyone/Michael_Thomas/df.bQTL_gene_partitioning.significant.",as.character(Sys.time()), ".csv", sep = ""))
#   write.csv(l.bQTL_gene_partitioning[[2]], paste("/Shared/Everyone/Michael_Thomas/df.bQTL_gene_partitioning.nonsignificant.",as.character(Sys.time()), ".csv", sep = ""))
#   
#   
#   df.rnaseq.test <- read.table("/Shared/Everyone/Michael_Thomas/B73_Mo17_rnaseq/ABS_significant_expression_qTELLER.txt", header = TRUE, sep = "")
#   
#   v.diff <- df.rnaseq.test$B73_Shoot / df.rnaseq.test$Mo17_Shoot
#   names(v.diff) <- df.rnaseq.test$gene_name
#   v.diff <- v.diff[v.diff >= 1.5]
#   v.diff <- v.diff[!is.na(v.diff)]
#   v.diff <- v.diff[v.diff != Inf]
#   
#   write.csv(names(v.diff), paste("/Shared/Everyone/Michael_Thomas/B73_Mo17_rnaseq/genes_ratio_B73_Mo17_ratio.",as.character(Sys.time()), ".csv", sep = ""))
#   
#   test <- subset(l.bQTL_gene_partitioning[[1]], l.bQTL_gene_partitioning[[1]]$promoter_5kb == "yes" | l.bQTL_gene_partitioning[[1]]$five_prime_UTR == "yes"  | l.bQTL_gene_partitioning[[1]]$post_gene_1kb == "yes")
#   #write.csv(unique(test$gene.ID), paste("/Shared/Everyone/Michael_Thomas/df.bQTL_genes.significant.",as.character(Sys.time()), ".csv", sep = ""))
#   #df.rnaseq.test <- subset(df.rnaseq.test, df.rnaseq.test$gene_name %in% unique(test$gene.ID))
#   
#   
#   
#   
#   
#   
#   test <- subset(l.bQTL_gene_partitioning[[2]], l.bQTL_gene_partitioning[[2]]$promoter_5kb == "yes" | l.bQTL_gene_partitioning[[2]]$five_prime_UTR == "yes" | l.bQTL_gene_partitioning[[2]]$post_gene_1kb == "yes")
#   #write.csv(unique(test$gene.ID), paste("/Shared/Everyone/Michael_Thomas/df.bQTL_genes.nonsignificant.",as.character(Sys.time()), ".csv", sep = ""))
#   





#   df.bQTL_gene_partitioning.set <- df.bQTL_gene_partitioning[,c(1:24,26)]
#   df.bQTL_gene_partitioning.set <- subset(df.bQTL_gene_partitioning.set, length(df.bQTL_gene_partitioning.set[,17:24] == "yes"))
#   
#   # write.csv(df.bQTL_gene_partitioning, "/Shared/Everyone/Michael_Thomas/df.bQTL_gene_partitioning_0111.csv")
#   
#   
#   test <- apply(df.bQTL_gene_partitioning.set, 1, function(m) {length(which(m[c(17:18,20:24)] == "yes"))})
#   v.gns.bQTL <- unique(df.bQTL_gene_partitioning.set$gene.ID)
#saveRDS(l.bQTL_gene_partitioning, "l.bQTL_gene_partitioning.rds")
# message("bQTL to gene distribution ")

# df.rnaseq.thomas.k <- subset(df.rnaseq.thomas, df.rnaseq.thomas$id %in% v.gns.DE.k)


### begin no partition - bqtl in eqtl - targets - rnaseq 


##### GWAS trait specific analysis (01.19.17)

if(b.GWAS){
  
  if(!b.RNAseq){ # 
    df.bQTL_gene_partitioning <- l.bQTL_gene_partitioning[[1]]
    # df.bQTL_gene_partitioning <- subset(df.bQTL_gene_partitioning, df.bQTL_gene_partitioning$post_gene_5kb == "no")
    # only annotate - exclusive 
    df.gene_function_partitioning <- merge(df.bQTL_gene_partitioning, df.gene_function, by = "gene.ID", all.x = TRUE)
  }
  
  ####
  
  # v.traits <- unique(df.geneSet.gwas.set$trait)
  
  l.bQTL_gene_partitioning_with_gwas <- vector(mode = "list", length = 2)
  l.number_per_trait <- vector(mode = "list", length = 2)
  l.number_per_trait[[1]] <- l.number_per_trait[[2]] <- numeric(length(v.traits))
  names(l.number_per_trait[[1]]) <- names(l.number_per_trait[[2]]) <- v.traits
  
  
  for(s in 1:2){
    
    df.bQTL_gene_partitioning <- l.bQTL_gene_partitioning[[s]]
    
    # add trait annotation to ASBs
    if(s == 1)
      df.bQTL_gene_partitioning[v.traits] <- "no"
    
    
    df.bQTL_gene_partitioning["dist_ASB_to_gwas"] <- NA
    
    
    # dynamic extension for traits
    df.bQTL_gene_partitioning_with_gwas <- c() 
    
    message("number of phenotype gwas: ", nrow(df.phenotype_gwas))
    
    # pb <- txtProgressBar(min = 0, max = length(v.traits), style = 3)
    
    #for(k in 1:length(v.traits)){ # trait
    
    # setTxtProgressBar(pb, k)
    
    for(i in 1:n.chromosomes){
      
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
              # df.bQTL_gene_partitioning.i$gwas_trait[j] <- paste(df.bQTL_gene_partitioning.i$gwas_trait[j], paste(df.phenotype_gwas.i.k$trait[j.set], collapse = ", "), collapse = "; ")
              # df.bQTL_gene_partitioning.i$dist_bQTL_to_GWAS_snp[j] <- paste(df.bQTL_gene_partitioning.i$dist_bQTL_to_GWAS_snp[j], paste(abs(df.bQTL_gene_partitioning.i$position[j] - df.phenotype_gwas.i.k$pos[j.set]), collapse = ", "), collapse = "; ")  
              
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
    #df.bQTL_gene_partitioning_with_gwas <- subset(df.bQTL_gene_partitioning_with_gwas, df.bQTL_gene_partitioning_with_gwas$gwas_trait != "")
    l.bQTL_gene_partitioning_with_gwas[[s]] <- df.bQTL_gene_partitioning_with_gwas
  }
  
  
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
    
    #tab <- matrix(c(hitInSample, hitInPop, sampleSize - hitInSample, popSize - hitInPop), nrow = 2)
    #print(paste("equal size comparison: FC:", foldchange," , p-value", fisher.test(tab)$p.value))
    
    df.enrichment$id[k] <- v.traits[k]
    #df.domain_enrichment$p.value[k] <- fisher.test(tab)$p.value
    df.enrichment$p.value[k] <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
    df.enrichment$foldchange[k] <-  as.numeric(foldchange)
    
  }
  
  subset(df.enrichment, df.enrichment$p.value < 0.05)
  
  df.enrichment <- subset(df.enrichment, df.enrichment$p.value < 0.1)
  df.enrichment
  
  write.csv(df.enrichment, paste("C:/Users/Michael/Desktop/HaschSeq_data/F7.csv", sep = ""))  
  
  df.enrichment <- subset(df.enrichment, df.enrichment$foldchange != Inf)
  
  df.bQTL_gene_partitioning_with_gwas <- df.bQTL_gene_partitioning
  
  if(b.RNAseq){
    write.csv(l.bQTL_gene_partitioning_with_gwas[[1]], paste("output/df.bQTL_gene_partitioning_with_RNAseq_with_gwas.",timeStamp, ".csv", sep = ""))
  }else{
    write.csv(l.bQTL_gene_partitioning_with_gwas[[1]], paste("output/df.bQTL_gene_partitioning_with_gwas.",timeStamp, ".csv", sep = ""))  
  }
  
  # write.csv(df.bQTL_gene_partitioning_with_gwas, paste("/Shared/Everyone/Michael_Thomas/results/df.bQTL_gene_partitioning_with_RNAseq_with_gwas.",as.character(Sys.time()), ".csv", sep = ""))
  
  # add gene ontology
  if(b.GO){
    
    df.bQTL_gene_partitioning_with_gwas["GO"] <- NA
    
    # the big gene ontology dataset
    df.GO.annot <- readRDS("datasets_paper/GeneOntology/df.GO_FULL_Maize_Trimmed.rds")
    #df.GO.annot <- readRDS("df.GO_Maize_Trimmed.rds")
    df.GO.annot <- unique(df.GO.annot)
    df.GO.annot <- subset(df.GO.annot, df.GO.annot$Ontology == "BP")
    
    pb <- txtProgressBar(min = 0, max = nrow(df.bQTL_gene_partitioning_with_gwas), style = 3)
    
    for(j in 1:nrow(df.bQTL_gene_partitioning_with_gwas)){
      
      setTxtProgressBar(pb, j)
      
      df.GO.annot.set <- subset(df.GO.annot, df.GO.annot$acc. == df.bQTL_gene_partitioning_with_gwas$gene.ID[j])
      df.GO.annot.set <- unique(df.GO.annot.set)
      
      df.bQTL_gene_partitioning_with_gwas$GO[j] <- paste(unique(df.GO.annot.set$Term), collapse = ", ")
      
    }
    
    close(pb)
    
    write.csv(df.bQTL_gene_partitioning_with_gwas, paste("output/df.bQTL_gene_partitioning_plus_rnaseq_gwas_go.",timeStamp, ".csv", sep = ""))
    #write.csv(df.bQTL_gene_partitioning_with_gwas, paste("/Shared/Everyone/Michael_Thomas/results/df.bQTL_gene_partitioning_plus_rnaseq_gwas_go.",as.character(Sys.time()), ".csv", sep = ""))
    
    ### most enriched go processes 
    
    
    # BARPLOTS mit STERNCHEN (if p < 0.05) - identify biological processes
    df.GO.annot.set <- subset(df.GO.annot, df.GO.annot$acc. %in% df.bQTL_gene_partitioning_with_gwas$gene.ID)
    
    df.tb.annot <- as.data.frame(table(df.GO.annot$Term[df.GO.annot$Term != ""]), stringsAsFactors = FALSE)
    df.tb.annot.set <- as.data.frame(table(df.GO.annot.set$Term[df.GO.annot.set$Term != ""]), stringsAsFactors = FALSE)
    df.tb.annot.set <- merge(df.tb.annot.set, df.tb.annot, by = "Var1")
    
    ## filter go - evidence codes, biological process
    df.BP_enrichment <- data.frame(BP = character(nrow(df.tb.annot.set)), p.val = numeric(nrow(df.tb.annot.set)),
                                   percent = numeric(nrow(df.tb.annot.set)), genes = numeric(nrow(df.tb.annot.set)), 
                                   foldchange = numeric(nrow(df.tb.annot.set)), stringsAsFactors = FALSE)
    n.genes.sset <- length(unique(df.GO.annot.set$acc.))
    n.genes <- length(unique(df.GO.annot$acc.))
    
    for(j in 1:nrow(df.tb.annot.set)){  
      
      n.inset <- df.tb.annot.set$Freq.x[j]
      n.genomewide <- df.tb.annot.set$Freq.y[j]
      
      #n.genomewide <- df.tb.annot$Freq[df.tb.annot$Var1 == df.tb.annot$Var1[j]]
      
      ### global enrichment test - gene basis
      hitInSample <- n.inset
      sampleSize <- n.genes.sset
      hitInPop <- n.genomewide #sum(tb.rate_limiting_domains$Freq)
      failInPop <- n.genes - hitInPop #(nrow(df.global.domains) - hitInPop)
      p.val <- print(phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE))
      
      #fisher.test(matrix(c(hitInSample-1, hitInPop, failInPop, sampleSize), 2, 2), alternative='less');
      foldChange <- (hitInSample / sampleSize) / (hitInPop / failInPop)
      #   mat.count <- matrix(c(n.inset,n.genes.sset - n.inset,  n.genomewide ,n.genes - n.genomewide), ncol = 2, byrow = FALSE)
      #   p.val <- fisher.test(mat.count)$p.value
      df.BP_enrichment$BP[j] <- df.tb.annot.set$Var1[j]
      df.BP_enrichment$p.val[j] <- p.val
      df.BP_enrichment$percent[j] <- (n.inset / n.genes.sset)
      df.BP_enrichment$genes[j] = n.inset
      
      df.BP_enrichment$foldchange[j] = foldChange
      
    }
    
    df.BP_enrichment <- df.BP_enrichment[order(df.BP_enrichment$p.val),]
    df.BP_enrichment <- df.BP_enrichment[order(df.BP_enrichment$foldchange),]
    df.BP_enrichment <- subset(df.BP_enrichment, df.BP_enrichment$p.val <= 0.05)
    
    
    write.csv(df.BP_enrichment, paste("output/df.EnrichedBiologicalProcesses_bQTL_gene_partitioning_with_RNAseq_with_gwas.",timeStamp, ".csv", sep = ""))
    # write.csv(df.BP_enrichment, paste("/Shared/Everyone/Michael_Thomas/results/df.EnrichedBiologicalProcesses_bQTL_gene_partitioning_with_RNAseq_with_gwas.",as.character(Sys.time()), ".csv", sep = ""))
    
  }
  
  
  
  
  
  
  
  
  # if GWAS - use the target or non target genes - for RNAseq
  
  v.traits <- unique(df.phenotype_gwas$trait)
  postTotal.significant <- l.selection[[1]]
  df.geneSet.gwas <- c()
  for(i in 1:n.chromosomes){
    df.phenotype_gwas.i <- subset(df.phenotype_gwas, df.phenotype_gwas$chr == i)  
    df.phenotype_gwas.i["dist_to_nearest_bqtl"] <- NA
    postTotal.significant.i <- subset(postTotal.significant, postTotal.significant$contig == i)
    start <- df.gene_annotation$pos.start
    end <- df.gene_annotation$pos.stop
    for(j in 1:nrow(df.phenotype_gwas.i)){
      # gwas snp not in genes
      if(all(start < df.phenotype_gwas.i$pos[j] & df.phenotype_gwas.i$pos[j] < end) == FALSE){
        df.phenotype_gwas.i$dist_to_nearest_bqtl[j]  <- min(abs(df.phenotype_gwas.i$pos[j] - postTotal.significant.i$position))
      }
    }
    df.geneSet.gwas <- rbind(df.geneSet.gwas, df.phenotype_gwas.i)
  }
  
  df.geneSet.gwas.set <- df.geneSet.gwas
  df.geneSet.gwas.set <- subset(df.geneSet.gwas, df.geneSet.gwas$dist_to_nearest_bqtl <= 2000)
  
  
  
  ####
  
  
  # trait specific 
  df.enrichment <- data.frame(id = character(length(v.traits)), p.value = numeric(length(v.traits)),  p.value.under = numeric(length(v.traits)), foldchange = numeric(length(v.traits)), stringsAsFactors = FALSE)
  
  pb <- txtProgressBar(min = 0, max = length(v.traits), style = 3)
  
  for(k in 1:length(v.traits)){ # trait
    
    v.sets <- numeric(2)
    
    setTxtProgressBar(pb, k)
    
    for(s in 1:2){
      
      if(s == 1)
        postTotal.significant <- l.selection[[s]]
      else
        postTotal.significant <- l.postTotal[[s]]
      
      df.geneSet.snp <- c()
      df.geneSet2.snp <- c()
      df.geneSet.gwas <- c()
      df.geneSet2.gwas <- c()
      df.snp_trans <- c()
      
      i.set.center <- 0
      v.hs.template <- numeric(1000)
      
      v.sizes.all <- c()
      
      # gwas nicht in genen => sig bqtl in 2000 kp
      # enrichment of traits 
      
      for(i in 1:n.chromosomes){
        
        df.phenotype_gwas.i <- subset(df.phenotype_gwas, df.phenotype_gwas$chr == i)  
        df.phenotype_gwas.i <- subset(df.phenotype_gwas.i, df.phenotype_gwas.i$trait == v.traits[k])
        
        if(nrow(df.phenotype_gwas.i) > 0){
          
          df.phenotype_gwas.i["dist_to_nearest_bqtl"] <- NA
          
          df.gff <-  df.gene_annotation
          
          df.gff.i <- subset(df.gff, df.gff$chr == i)
          df.gff.i["POSTfreq"] <- NA
          df.gff.i["snp"] <- NA
          df.gff.i["gwas"] <- NA
          df.gff.i["trait"] <- NA
          
          if(FALSE){ # subset to differential expressed target genes
            v.gns.tmp <- unique(subset(df.rnaseq.s, df.rnaseq.s$diffExp == TRUE)$gene.ID)
            df.gff.i <- subset(df.gff.i, df.gff.i$gene.ID %in% v.gns.tmp) 
          }
          
          postTotal.significant.i <- subset(postTotal.significant, postTotal.significant$contig == i)
          
          if(TRUE){
            start <- df.gff.i$pos.start
            end <- df.gff.i$pos.stop
            for(j in 1:nrow(df.phenotype_gwas.i)){
              # gwas snp not in genes
              if(any(start < df.phenotype_gwas.i$pos[j] & df.phenotype_gwas.i$pos[j] < end) == FALSE){
                df.phenotype_gwas.i$dist_to_nearest_bqtl[j]  <- min(abs(df.phenotype_gwas.i$pos[j] - postTotal.significant.i$position))
              }
            }
            df.geneSet.gwas <- rbind(df.geneSet.gwas, df.phenotype_gwas.i)
          }
        }
      }
      df.geneSet.gwas <- subset(df.geneSet.gwas, !is.na(df.geneSet.gwas$dist_to_nearest_bqtl))
      df.geneSet.gwas.set <- subset(df.geneSet.gwas, df.geneSet.gwas$dist_to_nearest_bqtl <= 2000) 
      v.sets[[s]] <- nrow(df.geneSet.gwas.set)
    } 
    
    hitInSample <- v.sets[[1]]
    sampleSize <- nrow(l.postTotal[[1]])
    hitInPop <- v.sets[[2]] #457 #874 
    popSize <-  nrow(l.postTotal[[2]])
    
    failInPop <- popSize - hitInPop
    
    foldchange <- (hitInSample/sampleSize)/(hitInPop/popSize)
    #phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE
    
    df.enrichment$id[k] <- v.traits[k]
    df.enrichment$p.value[k] <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
    df.enrichment$p.value.under[k] <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= TRUE)
    df.enrichment$foldchange[k] <-  as.numeric(foldchange)
    
  }
  
  
  df.enrichment.over <- subset(df.enrichment, df.enrichment$p.value <= 0.05)
  df.enrichment.under <- subset(df.enrichment, df.enrichment$p.value.under < 0.05)
  
  
  
  v.traits <- unique(df.geneSet.gwas.set$trait)
  df.enrichment <- data.frame(id = character(length(v.traits)), p.value = numeric(length(v.traits)), foldchange = numeric(length(v.traits)), stringsAsFactors = FALSE)
  
  for(k in 1:length(v.traits)){
    
    hitInSample <- nrow(subset(df.geneSet.gwas.set, df.geneSet.gwas.set$trait == v.traits[k]))
    sampleSize <- nrow(df.geneSet.gwas.set)
    hitInPop <- nrow(subset(df.phenotype_gwas, df.phenotype_gwas$trait == v.traits[k]))
    popSize <- nrow(df.phenotype_gwas)
    
    failInPop <- popSize - hitInPop
    
    foldchange <- (hitInSample/sampleSize)/(hitInPop/popSize)
    
    #tab <- matrix(c(hitInSample, hitInPop, sampleSize - hitInSample, popSize - hitInPop), nrow = 2)
    #print(paste("equal size comparison: FC:", foldchange," , p-value", fisher.test(tab)$p.value))
    
    df.enrichment$id[k] <- v.traits[k]
    #df.domain_enrichment$p.value[k] <- fisher.test(tab)$p.value
    df.enrichment$p.value[k] <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
    df.enrichment$foldchange[k] <-  as.numeric(foldchange)
    
  }
  
  df.enrichment <- subset(df.enrichment, df.enrichment$p.value < 0.5)
  df.enrichment <- subset(df.enrichment, df.enrichment$foldchange != Inf)
  
}


# 
# l.dist_bQTL_in_eQTL <- vector(mode = "list", length = 2)
# l.dist_bQTL_in_eQTL[[1]] <- numeric()
# l.dist_bQTL_in_eQTL[[2]] <- numeric()
# 
# for(s in 1:2){
#   
#   postTotal.significant <- l.selection[[s]]
#   
#   #postTotal.significant <- l.postTotal[[s]]
#   
# 
#   v.dist_bQTL_in_eQTL <- numeric()
#   
#   v.bQTLL_gwas_in_eqtl <- 0
#   
#   print(nrow(postTotal.significant))
#   
#   n.bQTL_in_eQTL <- 0
#   
#   for(i in 1:10){
#     
#     df.eqtls.liu.agp3.i <- subset(df.eqtls.liu.agp3, df.eqtls.liu.agp3$mapped_id == i)
#     df.eqtls.liu.agp3.i <- unique(df.eqtls.liu.agp3.i[,2:4])
#     postTotal.significant.i <- subset(postTotal.significant, postTotal.significant$contig == i)
#     df.phenotype_gwas.i <- subset(df.phenotype_gwas, df.phenotype_gwas$V5 == i) 
#     
#     for(j in 1:nrow(df.eqtls.liu.agp3.i)){
#       i.set <- which(postTotal.significant.i$position > df.eqtls.liu.agp3.i$mapped_start[j] & postTotal.significant.i$position < df.eqtls.liu.agp3.i$mapped_stop[j])
#       n.bQTL_in_eQTL <- n.bQTL_in_eQTL + length(i.set)
#     }
#   
#   }
#   
#   print(n.bQTL_in_eQTL)
#   #l.dist_bQTL_in_eQTL[[s]] <- v.dist_bQTL_in_eQTL
#   
# }
# 
# 
# ###
# 
# l.dist_bQTL_in_eQTL <- vector(mode = "list", length = 2)
# l.dist_bQTL_in_eQTL[[1]] <- numeric()
# l.dist_bQTL_in_eQTL[[2]] <- numeric()
# 
# for(s in 1:2){
#   
#   postTotal.significant <- l.postTotal[[s]]
#   
#   #postTotal.significant <- postTotal.significant[1:2000,]
#   
#   #df.bQTL_in_eQTL <- c()
#   
#   # gwas nicht in genen => sig bqtl in 2000 kp
#   # enrichment of traits 
#   
#   v.dist_bQTL_in_eQTL <- numeric()
#   
#   v.bQTLL_gwas_in_eqtl <- 0
#   
#   for(i in 1:10){
#   
#     df.geneSet.snp_gwas.i <- subset(df.geneSet.snp_gwas, df.geneSet.snp_gwas$V1 == i)
#     snps <- as.numeric(unlist(sapply(df.geneSet.snp_gwas.i$snp, function(m) strsplit(m," ")))  )
#     
#     df.eqtl.agpv3.i <- subset(df.eqtl.agpv3, df.eqtl.agpv3$Lead_SNP_chr_AGPv3 == i)
#     postTotal.significant.i <- subset(postTotal.significant, postTotal.significant$contig == i)
#     
#     start <- df.eqtl.agpv3.i$eQTL_region_start_AGPv3
#     end <- df.eqtl.agpv3.i$eQTL_region_end_AGPv3
#     pos_leadQTL <- df.eqtl.agpv3.i$Lead_SNP_pos_AGPv3
#     
#     for(k in 1:nrow(postTotal.significant.i)){
#       
#       if(any(start < postTotal.significant.i$position[k] & postTotal.significant.i$position[k] < end) == TRUE){
#         
#         v.dist_bQTL_in_eQTL  <- c(v.dist_bQTL_in_eQTL, min(abs(pos_leadQTL - postTotal.significant.i$position[k])))
#         
#         if(postTotal.significant.i$position[k] %in% as.numeric(snps))
#           v.bQTLL_gwas_in_eqtl <- v.bQTLL_gwas_in_eqtl + 1
#         
#       }
#     }
#   }
#   
#   l.dist_bQTL_in_eQTL[[s]] <- v.dist_bQTL_in_eQTL
#   
# }
# 
# 
# 
# 
# lapply(l.dist_bQTL_in_eQTL, length)
# 
# (length(l.dist_bQTL_in_eQTL[[1]][l.dist_bQTL_in_eQTL[[1]] < 2000]) / nrow(l.postTotal[[1]])) / (length(l.dist_bQTL_in_eQTL[[2]][l.dist_bQTL_in_eQTL[[2]] < 2000]) / nrow(l.postTotal[[2]]))
# 
# (length(l.dist_bQTL_in_eQTL[[1]]) / nrow(l.postTotal[[1]])) / (length(l.dist_bQTL_in_eQTL[[2]]) / nrow(l.postTotal[[2]]))
# 
# 
# 
# 
#     
#     start <- df.gff.i$V4
#     end <- df.gff.i$V5
#     
#     
#     #       
#     #       start <- c(1,5)
#     #       end <- c(5, 11)
#     #       any(start < 5 & 5 < end)
#     #       
#     for(j in 1:nrow(df.phenotype_gwas.i)){
#       
#       # gwas snp not in genes
#       if(any(start < df.phenotype_gwas.i$V13[j] & df.phenotype_gwas.i$V13[j] < end) == FALSE){
#         
#         df.phenotype_gwas.i$dist_to_nearest_bqtl[j]  <- min(abs(df.phenotype_gwas.i$V13[j] - postTotal.significant.i$position))
#         
#       }
#       
#     }
#     
#     
#     
#     # eQTL 
#     
#     df.phenotype_gwas.i <- subset(df.phenotype_gwas, df.phenotype_gwas$V5 == i)  
#     df.phenotype_gwas.i <- subset(df.phenotype_gwas.i, df.phenotype_gwas.i$trait == v.traits[k])
#     
#     
#     df.phenotype_gwas.i["dist_to_nearest_bqtl"] <- NA
#     
#     
#     df.eqtl.agpv3
#     
# 
#     postTotal.significant.i
# 
#   }
#   
#   
# }
# 







## trait speciifc 
df.enrichment <- data.frame(id = character(length(v.traits)), p.value = numeric(length(v.traits)),  p.value.under = numeric(length(v.traits)), foldchange = numeric(length(v.traits)), stringsAsFactors = FALSE)

pb <- txtProgressBar(min = 0, max = length(v.traits), style = 3)

for(k in 1:length(v.traits)){
  
  v.sets <- numeric(2)
  
  setTxtProgressBar(pb, k)
  
  for(s in 1:2){
    
    if(s == 1)
      postTotal.significant <- l.selection[[s]]
    else
      postTotal.significant <- l.postTotal[[s]]
    
    
    # postTotal.significant <- df.postTotal_both.selection
    #postTotal.significant <- postTotal.significant[1:2000,]
    
    df.geneSet.snp <- c()
    df.geneSet2.snp <- c()
    df.geneSet.gwas <- c()
    df.geneSet2.gwas <- c()
    df.snp_trans <- c()
    
    i.set.center <- 0
    v.hs.template <- numeric(1000)
    
    v.sizes.all <- c()
    
    # gwas nicht in genen => sig bqtl in 2000 kp
    # enrichment of traits 
    
    for(i in 1:10){
      
      # eQTL 
      
      df.phenotype_gwas.i <- subset(df.phenotype_gwas, df.phenotype_gwas$V5 == i)  
      df.phenotype_gwas.i <- subset(df.phenotype_gwas.i, df.phenotype_gwas.i$trait == v.traits[k])
      
      
      # print(i)
      
      #     df.hypersensitivity.i <- subset(df.hypersensitivity, df.hypersensitivity$chr == i)
      #     
      #     v.sizes <- df.hypersensitivity.i$end - df.hypersensitivity.i$start
      #     v.sizes.all <- c(v.sizes.all, v.sizes)
      #     
      #     df.hypersensitivity.i["start_margin"] <- df.hypersensitivity.i$start + v.sizes / 2 - 500
      #     df.hypersensitivity.i["end_margin"] <- df.hypersen00sitivity.i$end - v.sizes / 2 + 500
      
      
      #     
      #df.report_gwas.i <- subset(df.report_gwas, df.r0eport_gwas$V5 == i)
      df.phenotype_gwas.i <- subset(df.phenotype_gwas, df.phenotype_gwas$V5 == i)  
      df.phenotype_gwas.i <- subset(df.phenotype_gwas.i, df.phenotype_gwas.i$trait == v.traits[k])
      
      
      if(nrow(df.phenotype_gwas.i) > 0){
        
        #     df.phenotype_gwas.i <- test
        #     names(df.phenotype_gwas.i) <- c("V5", "V13")
        #     df.phenotype_gwas.i <- subset(df.phenotype_gwas.i, df.phenotype_gwas.i$V5 == i)  
        #     
        #     df.phenotype_gwas.i <- l.postTotal[[s]][2000:4000,]
        #     df.phenotype_gwas.i <- subset(df.phenotype_gwas.i, df.phenotype_gwas.i$contig == i) 
        #     names(df.phenotype_gwas.i) <- c("V5", "V13")
        
        df.phenotype_gwas.i["dist_to_nearest_bqtl"] <- NA
        
        df.gff <- l.gffSets[[1]]
        
        df.gff.i <- subset(df.gff, df.gff$V1 == i)
        df.gff.i["POSTfreq"] <- NA
        df.gff.i["snp"] <- NA
        df.gff.i["gwas"] <- NA
        df.gff.i["trait"] <- NA
        
        #     df.gff.i2 <- subset(df.gff, df.gff$V1 == i)
        #     df.gff.i2["POSTfreq"] <- NA
        #     df.gff.i2["snp"] <- NA
        #     df.gff.i2["gwas"] <- NA
        #     df.gff.i2["trait"] <- NA
        
        postTotal.significant.i <- subset(postTotal.significant, postTotal.significant$contig == i)
        #postTotal.significant.i <- postTotal.significant.i[order(postTotal.significant.i$position),]
        
        if(TRUE){
          
          start <- df.gff.i$V4
          end <- df.gff.i$V5
          
          for(j in 1:nrow(df.phenotype_gwas.i)){
            # gwas snp not in genes
            if(any(start < df.phenotype_gwas.i$V13[j] & df.phenotype_gwas.i$V13[j] < end) == FALSE){
              df.phenotype_gwas.i$dist_to_nearest_bqtl[j]  <- min(abs(df.phenotype_gwas.i$V13[j] - postTotal.significant.i$position))
            }
          }
          
          #       for(j in 1:nrow(df.gff.i)){
          #         
          #         if(df.gff.i$V7[j] == "+"){
          #         
          #           start <- df.gff.i$V4[j]
          #           end <- df.gff.i$V5[j]
          #           
          #         }else if(df.gff.i$V7[j] == "-"){
          #           
          #           start <- df.gff.i$V4[j]
          #           end <- df.gff.i$V5[j] 
          #           
          #         }
          #         
          #         # in 2 kb of binding qtls
          #         
          #         # gwas snp not in 
          #         if(all(start < df.phenotype_gwas.i$pos & df.phenotype_gwas.i$pos < end) == FALSE){
          #           
          #         }
          #         
          #         
          #         df.phenotype_gwas.i
          #         
          #         
          #         
          #         #i.set <- which(start < df.phenotype_gwas.i$pos & df.phenotype_gwas.i$pos < end)
          #         
          #         # i.set <- which(start < df.report_gwas.i$V13 & df.report_gwas.i$V13 < end)
          #         #i.set <- which(start < postTotal.significant.i$position & postTotal.significant.i$position < end)
          #       
          #         #         if(length(i.set) > 0){
          #         #           df.gff.i$gwas[j] <-  paste(df.phenotype_gwas.i$pos[i.set], collapse = " ")
          #         #           #df.gff.i$trait[j] <- paste(subset(df.phenotype_gwas.i, df.phenotype_gwas.i$pos %in% df.report_gwas.i$V8[i.set])$trait, collapse = " ")
          #         #         }
          #         #         
          #         #         i.set <- which(start < postTotal.significant.i$position & postTotal.significant.i$position < end)
          #         #         
          #         #         if(length(i.set) > 0){
          #         #           df.gff.i$snp[j] <-  paste(postTotal.significant.i$position[i.set], collapse = " ")
          #         #           df.gff.i$POSTfreq[j] <-  paste(round(postTotal.significant.i$POSTfreq[i.set], 2), collapse = " ")
          #         #         }
          #           
          #       }
          
          #df.phenotype_gwas.i <- df.phenotype_gwas.i[i.set,]
          
          #       df.report_gwas.i["trait"] <- df.report_gwas.i$V13 
          #       
          #       i.set <- which(df.phenotype_gwas.i$pos %in% df.report_gwas.i$V8)
          #       
          #       
          #       i.set
          #       
          #       df.phenotype_gwas.i$trait
          
          df.geneSet.gwas <- rbind(df.geneSet.gwas, df.phenotype_gwas.i)
          
          
          
          #df.geneSet2.gwas <- rbind(df.geneSet2.gwas, df.gff.i2)
          #df.geneSet.gwas <- rbind(df.geneSet.gwas, df.gff.i)
          
        }
      }
      
      ## highlight distance of significant and non-significant qtls ## 
      
      
      
      #     df.snp.pos <- readRDS(paste("df.snp.pos_",i, ".rds"))
      #     
      #     # remove heterozygote snps
      #     idx.snps <- which(lapply(df.snp.pos@elementMetadata$ALT, length) == 1)
      #     
      #     vec.snp.pos <- as.numeric(df.snp.pos@ranges@start[idx.snps])
      #     
      #     vec.snp.pos <- vec.snp.pos[sample(length(vec.snp.pos), 10000)]
      #     
      #     v.positions <- vec.snp.pos
      
      if(FALSE){
        
        v.positions <- postTotal.significant.i$position
        
        for(j in 1:length(v.positions)){
          
          i.set <- which(df.hypersensitivity.i$start < v.positions[j] & v.positions[j] < df.hypersensitivity.i$end)  
          
          if(length(i.set) > 0)
            i.set.center <- i.set.center + 1
          
          i.set <- which(df.hypersensitivity.i$start_margin < v.positions[j] & v.positions[j] < df.hypersensitivity.i$end_margin)  
          
          if(length(i.set) > 0){
            
            v.center <- (df.hypersensitivity.i$end_margin[i.set] - df.hypersensitivity.i$start_margin[i.set])/2 + df.hypersensitivity.i$start_margin[i.set]
            
            v.dist.center <- (v.positions[j] - v.center)
            
            i.set <- which(abs(v.dist.center) == min(abs(v.dist.center)))
            
            i.set <- 500 + v.dist.center[i.set]
            
            v.hs.template[i.set] <- v.hs.template[i.set] + 1
            
          }
          
        }
        
        
        
      }
      
      
      #nrow(postTotal.significant.i)
      
      #(355 / 1875) / (3086/ 15988)
      
      #v.hs.template.sig <- v.hs.template
      
      #plot(seq(1:1000), v.hs.template, type = "l")
      
      #plot(seq(1:1000), v.hs.template.sig, type = "l")
      #lines(seq(1:1000),v.hs.template, lwd=2, col="red")
    }
    
    
    #df.geneSet.gwas.set <- df.geneSet.snp_gwas
    df.geneSet.gwas <- subset(df.geneSet.gwas, !is.na(df.geneSet.gwas$dist_to_nearest_bqtl))
    df.geneSet.gwas.set <- subset(df.geneSet.gwas, df.geneSet.gwas$dist_to_nearest_bqtl <= 2000) 
    
    
    v.sets[[s]] <- nrow(df.geneSet.gwas.set)
  } 
  
  # df.geneSet.gwas
  # df.motif_analysis.postprocessed.peak <- l.postTotal[[1]][,c(1,2,2)]
  # df.motif_analysis.postprocessed.peak$Chr <- paste("chr",df.motif_analysis.postprocessed.peak$Chr, sep = "")
  
  hitInSample <- v.sets[[1]]
  sampleSize <- nrow(l.postTotal[[1]])
  hitInPop <- v.sets[[2]] #457 #874 
  popSize <-  nrow(l.postTotal[[2]])
  
  failInPop <- popSize - hitInPop
  
  foldchange <- (hitInSample/sampleSize)/(hitInPop/popSize)
  #phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
  
  
  df.enrichment$id[k] <- v.traits[k]
  df.enrichment$p.value[k] <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
  df.enrichment$p.value.under[k] <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= TRUE)
  df.enrichment$foldchange[k] <-  as.numeric(foldchange)
  
}

close(pb)

df.enrichment.over <- subset(df.enrichment, df.enrichment$p.value <= 0.05)
df.enrichment.under <- subset(df.enrichment, df.enrichment$p.value.under < 0.05)

## 
v.growth_set <- c("Ratio of ear height to total height",  "Leaf length",  "Ear height", "Average internode length (above ear)", "Plant height", "Nodes per plant", "Cob diameter", "Nodes above ear")

hitInSample <- nrow(subset(df.geneSet.gwas.set, df.geneSet.gwas.set$trait %in% v.growth_set))
sampleSize <- nrow(df.geneSet.gwas.set)
hitInPop <- nrow(subset(df.phenotype_gwas, df.phenotype_gwas$trait %in% v.growth_set))
popSize <- nrow(df.phenotype_gwas)

failInPop <- popSize - hitInPop

foldchange <- (hitInSample/sampleSize)/(hitInPop/popSize)
phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)

#nrow(subset(df.geneSet.gwas.set, df.geneSet.gwas.set$trait %in% df.enrichment.over$id)) - hitInSample


write.csv(df.enrichment.over, paste("/Shared/Everyone/Michael_Thomas/GWAS/df.enrichment.over.",as.character(Sys.time()), ".csv", sep = ""))
write.csv(df.enrichment.under, paste("/Shared/Everyone/Michael_Thomas/GWAS/df.enrichment.under.",as.character(Sys.time()), ".csv", sep = ""))

df.enrichment

v.traits <- unique(df.geneSet.gwas.set$trait)

df.enrichment <- data.frame(id = character(length(v.traits)), p.value = numeric(length(v.traits)), foldchange = numeric(length(v.traits)), stringsAsFactors = FALSE)

for(k in 1:length(v.traits)){
  
  hitInSample <- nrow(subset(df.geneSet.gwas.set, df.geneSet.gwas.set$trait == v.traits[k]))
  sampleSize <- nrow(df.geneSet.gwas.set)
  hitInPop <- nrow(subset(df.phenotype_gwas, df.phenotype_gwas$trait == v.traits[k]))
  popSize <- 3873 # nrow(df.phenotype_gwas)
  
  failInPop <- popSize - hitInPop
  
  foldchange <- (hitInSample/sampleSize)/(hitInPop/popSize)
  
  #tab <- matrix(c(hitInSample, hitInPop, sampleSize - hitInSample, popSize - hitInPop), nrow = 2)
  #print(paste("equal size comparison: FC:", foldchange," , p-value", fisher.test(tab)$p.value))
  
  df.enrichment$id[k] <- v.traits[k]
  #df.domain_enrichment$p.value[k] <- fisher.test(tab)$p.value
  df.enrichment$p.value[k] <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
  df.enrichment$foldchange[k] <-  as.numeric(foldchange)
  
}

df.enrichment <- subset(df.enrichment, df.enrichment$p.value < 0.5)
df.enrichment <- subset(df.enrichment, df.enrichment$foldchange != Inf)

#write.csv(df.enrichment,  "/Shared/Everyone/Michael_Thomas/df.gwas_traitEnrichment_not_in_gene.csv")


df.geneSet.gwas.set


df.geneSet.gwas.set

write.table(df.geneSet.gwas.set[,c(2,3,3)], "/Shared/Everyone/Michael_Thomas/df.gwas_notInGene_2000bpSnpProximity_1211.bed", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
write.csv(df.geneSet.gwas.set, "/Shared/Everyone/Michael_Thomas/df.gwas_notInGene_2000bpSnpProximity_1211.csv")

# ernichment test over the resulting 


#df.geneSet.gwas. <- subset(df.geneSet.gwas, !is.na(df.geneSet.gwas$gwas))



#v.hs.template.bg <- v.hs.template
plot(seq(1:1000), v.hs.template, type = "l")
lines(seq(1:1000),v.hs.template.bg, lwd=2, col="red")


v.percentage.center <- numeric(500)

for(j in 1:500){
  v.percentage.center[j]  <- sum(v.hs.template[(500-j):(500+j)])/ sum(v.hs.template) * 100
}

plot(seq(1:500), v.percentage.center, type = "l")


v.hs.template

# subset to snps
df.geneSet.snp <- subset(df.geneSet.snp, !is.na(df.geneSet.snp$snp))

# subset to gwas
df.geneSet.gwas <- subset(df.geneSet.gwas, !is.na(df.geneSet.gwas$gwas))

length(intersect(df.geneSet.gwas$id, df.geneSet.snp$id))


df.gwas_snps <- subset(df.geneSet.gwas, !is.na(df.geneSet.gwas$snp) & !is.na(df.geneSet.gwas$snp))
gwas <- sapply(df.gwas_snps$trait, function(m) strsplit(m, " "))
gwas <- (unlist(gwas))


write.table(df.geneSet.snp, "/Shared/Everyone/Michael_Thomas/traitTable.txt", sep = "\t")


#test <- subset(df.geneSet.snp, !is.na(df.geneSet.snp$snp) & !is.na(df.geneSet.snp$snp))

test <- df.geneSet.snp
test["tair"] <- NA
test["annotation"] <- NA

for(j in 1:nrow(test)){
  
  annot.j <- subset(df.annotation, df.annotation$V2 == test$id[j])[1,]
  test$tair[j] <- annot.j$V11
  test$annotation[j] <- annot.j$V13
  
}

names(df.rnaseq)[1] <- "id"
test <- merge(x = test, y = df.rnaseq, by = "id", all = TRUE)

test <- subset(test, !is.na(test$snp))

test["diffExp"] <- ifelse(test$pvalue < 0.05, TRUE, FALSE)
test["mode"] <- ifelse(test$log2FoldChange >= 0.5 & test$pvalue < 0.05, "up", ifelse(test$log2FoldChange <= - 0.5 & test$pvalue < 0.05, "down", ""))


write.csv(test, "/Shared/Everyone/Michael_Thomas/ASEReadCounter/genes_with_binding_peaks_in_promoters_with_rnaseq_1205.csv")


tb.all <- table(df.phenotype_gwas$trait)
tb.set <- table(gwas)

th.pvalue <- 0.1


df.enrichment <- data.frame(id = character(length(tb.set)), p.value = numeric(length(tb.set)), foldchange = numeric(length(tb.set)), stringsAsFactors = FALSE)

for(k in 1:length(tb.set)){
  
  hitInSample <- tb.set[k]
  sampleSize <-  sum(tb.set) 
  hitInPop <- tb.all[names(tb.set)][k]
  popSize <- sum(tb.all[names(tb.set)])
  
  failInPop <- popSize - hitInPop
  foldchange <- (hitInSample/sampleSize)/(hitInPop/popSize)
  
  #tab <- matrix(c(hitInSample, hitInPop, sampleSize - hitInSample, popSize - hitInPop), nrow = 2)
  #print(paste("equal size comparison: FC:", foldchange," , p-value", fisher.test(tab)$p.value))
  
  df.enrichment$id[k] <- as.character(names(tb.set)[k])
  #df.domain_enrichment$p.value[k] <- fisher.test(tab)$p.value
  df.enrichment$p.value[k] <- phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
  df.enrichment$foldchange[k] <-  as.numeric(foldchange)
  
}

df.enrichment <- subset(df.enrichment, df.enrichment$p.value < th.pvalue)



# enrichment test over traits
table(gwas)

table(df.phenotype_gwas$trait)


# df.geneSet2.snp <- subset(df.geneSet2, !is.na(df.geneSet2$gwas))

gwas <- sapply(df.geneSet2$gwas, function(m) strsplit(m, ","))
gwas <- unique(unlist(gwas))
postTotal.significant.set <- subset(postTotal.significant, !postTotal.significant$position %in% gwas)  


#   df.geneSet <- subset(df.geneSet, !is.na(df.geneSet$snp))
#   df.geneSet2 <- subset(df.geneSet2, !is.na(df.geneSet2$snp))
#   
#   snp <- sapply(df.geneSet2$snp, function(m) strsplit(m, ","))
#   snp <- unique(unlist(snp))
#   postTotal.significant.set <- subset(postTotal.significant, !postTotal.significant$position %in% snp)

#write.table(postTotal.significant.set, "/Shared/Everyone/Michael_Thomas/postTotal.significant_nongenes.bed", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)





## bed format 
postTotal.significant.set <- postTotal.significant.set[,c(1,2,2)]
postTotal.significant.set$chr <- paste("chr",postTotal.significant.set$chr, sep = "")
# 
# 
write.table(postTotal.significant.set, "/Shared/Everyone/Michael_Thomas/postTotal.significant_nongenes.bed", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)





i.set <- which(
  
  
  v.mapping <- df.report_gwas.i$V13
  names(v.mapping) <- df.report_gwas.i$V8
  
  df.phenotype_gwas.i["pos.AGPv3"] <- v.mapping[as.character(df.phenotype_gwas.i$pos)]
  
  
  i.gwas <- which(df.phenotype_gwas.i$pos.AGPv3 %in% postTotal.significant.i$position)
  
  print(df.phenotype_gwas.i$trait[i.gwas])
  
  
  test <- rbind(test, df.phenotype_gwas.i[i.gwas,])
  
  #df.report_gwas.i$V13[i.pos]
  }



## bed format 
test
test <- test[,c(2,7,7)]
test$chr <- paste("chr",test$chr, sep = "")
# 
# 
write.table(test, "/Shared/Everyone/Michael_Thomas/GWAS/df.snp_gwas_nonsignificant.bed", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)


df.report_gwas.set <- df.report_gwas 
df.report_gwas.set <- df.report_gwas.set[,c(2,4,4)]
df.report_gwas.set$V5 <- paste("chr",df.report_gwas.set$V5, sep = "")
write.table(df.report_gwas.set, "/Shared/Everyone/Michael_Thomas/df.snp_gwas_all.bed", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

