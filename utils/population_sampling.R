create_background_distribution_with_similar_gene_partitioning <- function(l.bQTL_gene_partitioning = l.bQTL_gene_partitioning, multiplyer = 10, v.partitions = v.partitions, b.duplicateRemoval = FALSE, seed.randomGenerator = 1234){
  
  set.seed(seed.randomGenerator)
  
  df.bQTL_gene_partitioning = l.bQTL_gene_partitioning[[1]]
  
  v.partitions.duplicates <- v.partitions[c(1,2,8)]
  
  if(b.duplicateRemoval){
    # identify duplicate partitions and remove => enter into general pipeline
    for(i in 1:2){
      i.set <- which(apply(l.bQTL_gene_partitioning[[i]][,v.partitions.duplicates], 1, function(m) length(which(m == "yes"))) == 2)
      i.set <- (!1:nrow(l.bQTL_gene_partitioning[[i]]) %in% i.set)
      l.bQTL_gene_partitioning[[i]] <- l.bQTL_gene_partitioning[[i]][i.set, ]
    }
  }

  df.bQTL_gene_partitioning_bg_equivalent <- c()
  
  for(i in 1:n.chromosomes){
    
    df.bQTL_gene_partitioning.asb.i <- subset(l.bQTL_gene_partitioning[[1]],  l.bQTL_gene_partitioning[[1]]$contig == i)
    df.bQTL_gene_partitioning_bg.i  <- subset(l.bQTL_gene_partitioning[[2]],  l.bQTL_gene_partitioning[[2]]$contig == i)

    v.distribution <- apply(df.bQTL_gene_partitioning.asb.i[,v.partitions], 2, table)["yes",]
    v.distribution <- v.distribution * multiplyer
    
    df.bQTL_gene_partitioning_bg_equivalent.i <- c()
    
    for(j in 1:length(v.partitions)){
      if(v.partitions[j] != "gene"){
        df.bQTL_gene_partitioning_bg.i.set <- df.bQTL_gene_partitioning_bg.i[which(as.character(df.bQTL_gene_partitioning_bg.i[,v.partitions[j]]) == "yes"),]
        i.set <- sample(nrow(df.bQTL_gene_partitioning_bg.i.set), v.distribution[j])
        df.bQTL_gene_partitioning_bg_equivalent.i <- rbind(df.bQTL_gene_partitioning_bg_equivalent.i, df.bQTL_gene_partitioning_bg.i.set[i.set,])
      }
    }
    
    df.bQTL_gene_partitioning_bg_equivalent <- rbind(df.bQTL_gene_partitioning_bg_equivalent, df.bQTL_gene_partitioning_bg_equivalent.i)
    
  }
  
  v.distribution <- apply(l.bQTL_gene_partitioning[[1]][,v.partitions], 2, table)["yes",]
  v.distribution.bg_equivalent <- apply(df.bQTL_gene_partitioning_bg_equivalent[,v.partitions], 2, table)["yes",]

  df.distribution <- data.frame(ASBs = v.distribution, bgSNPs = v.distribution.bg_equivalent, multiplyer = (v.distribution.bg_equivalent / v.distribution))

  l.bQTL_gene_partitioning[[2]] <- df.bQTL_gene_partitioning_bg_equivalent
  
  
  return(list(l.bQTL_gene_partitioning = l.bQTL_gene_partitioning, df.distribution = df.distribution))
  
}



sample_nonsignificant_background_snps_in_peaks <- function(l.postTotal, seed.randomGenerator = 1234){
  
  set.seed(seed.randomGenerator)
  
  message("sampling non significant background snps in peaks")
  # independent of peak and linkage 
  
  #n.total <- nrow(l.postTotal[[1]]) * n.bgSnp.multiplier
  # l.postTotal[[2]] <- c() # data.frame(contig = character(), position = numeric())
  
  for(i in 1:n.chromosomes){
    
    if(TRUE){
      
      df.snp.pos <- readRDS(paste("tmp/df.snp.pos_",i, ".rds")) # needs to be
      df.asb.i <- subset(l.postTotal[[1]], l.postTotal[[1]]$contig == i)
      
      # remove heterozygote snps
      if(FALSE){
        idx.snps <- which(lapply(df.snp.pos@elementMetadata$ALT, length) == 1)  
      }else{
        vec.snp.bases <- CharacterList(df.snp.pos@elementMetadata$ALT)
        vec.snp.bases = unstrsplit(vec.snp.bases, sep = ",")
        vec.snp.bases <- gsub("\\,.*", "", vec.snp.bases)
        df.snp.pos@elementMetadata$ALT <- vec.snp.bases
      }
      
      # A - generic SNPs
      if(FALSE){
        vec.snp.pos <- as.numeric(df.snp.pos@ranges@start[idx.snps])  
        vec.snp.bases <- as.character(unlist(df.snp.pos@elementMetadata$ALT[idx.snps]))
      }else{
        vec.snp.pos <- as.numeric(df.snp.pos@ranges@start)  
        vec.snp.bases <- as.character(unlist(df.snp.pos@elementMetadata$ALT))
      }
      
      vec.snp.bases <- ifelse(vec.snp.bases == "<DEL>", "N",vec.snp.bases) # both represented separately
      vec.snp.bases <- ifelse(vec.snp.bases == "<INS>", "N",vec.snp.bases)
      # remove empty strings # 
      idx.snp.exceptions <- which(vec.snp.bases == "")
      if(length(idx.snp.exceptions) > 0){
        vec.snp.pos <- vec.snp.pos[-idx.snp.exceptions]
        vec.snp.bases <- vec.snp.bases[-idx.snp.exceptions]
      }
      
      n.set <- table(l.postTotal[[1]]$contig)[[i]] * n.bgSnp.multiplier
      
      # sampling without bATL 
      df.tmp <- data.frame(contig = rep(i, n.set), position = vec.snp.pos[sample(length(vec.snp.pos), n.set, replace = FALSE)])
      df.tmp <- subset(df.tmp, !df.tmp$position %in% df.asb.i$position)
      l.postTotal[[2]] <- rbind(l.postTotal[[2]], df.tmp)
      
    }else{
      postTotal.i <- subset(postTotal, postTotal$contig == i)
      #postTotal.i <- subset(postTotal.i, postTotal.i$`p-value (corrected)` >= 0.95)
      postTotal.i <- subset(postTotal.i, postTotal.i$`p-value (corrected)` >= 0.5)
      l.postTotal[[2]] <- rbind(l.postTotal[[2]], postTotal.i)
    }
  }
  
  return(l.postTotal=l.postTotal)
  
}


