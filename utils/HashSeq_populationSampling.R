
sample_nonsignificant_background_snps_in_peaks <- function(l.postTotal, seed.randomGenerator = 1234){
    
  set.seed(seed.randomGenerator)
  
  message("sampling non significant background snps in peaks")
  # independent of peak and linkage 
  
  #n.total <- nrow(l.postTotal[[1]]) * n.bgSnp.multiplier
  # l.postTotal[[2]] <- c() # data.frame(contig = character(), position = numeric())
  
  for(i in 1:n.chromosomes){
    
    if(TRUE){
      
      df.snp.pos <- readRDS(paste(folder_tmp, "df.snp.pos_",i, ".rds", sep ="/")) # needs to be
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


