rm(list=ls()) # clear workspace 


setwd("D:/junkDNA.ai/Projects/HASCHSEQ/HaschSeqBasePrototype")


source("Hashseq_Parameters.R")

source("utils/general.R")
install_and_load_libraries()

# define input datasets 
message("load and preprocess basic datasets...")
source("HachSeq_datasets.R")


if(b.firstRun){
  
  message("Build mutated genome from input SNPs")
  
  for(i in 1:n.chromosomes){
    
    print(paste("processing chromosome ", i))
    
    filename <- paste(path.snptables,v.mo17_files[i],sep="")
    
    # assemble the entire genome (all 12 chromosomes)
    vcf <- readVcf(filename,"BSgenome.Zmays.NCBI.AGPv4")
    
    df <- snpSummary(vcf)
    df.snp.pos <- rowRanges(vcf)
    
    if(b.doInitialHomozygousFilter){
      idx.homozygous <- which(df$g11 == 1) # remove all non homozygous
      df.snp.pos <- df.snp.pos[idx.homozygous,]
    }
    
    # subset to first alternative 
    saveRDS(df.snp.pos, paste("tmp/df.snp.pos_",i, ".rds"))
    
  }
  
  df.snp_positions <- c()
  l.genome.mutant <- vector(mode = "list", length = n.chromosomes)
  # welche haben laenge 
  for(i in 1:n.chromosomes){
    
    print(paste("processing chromosome ", i))
    
    df.snp.pos <- readRDS(paste("tmp/df.snp.pos_",i, ".rds"))
    
    # n.snps <- n.snps + length(df.snp.pos)
    
    # remove heterozygote snps
    # idx.snps <- which(lapply(df.snp.pos@elementMetadata$ALT, length) == 1)
    
    # A - generic SNPs
    #vec.snp.pos <- as.numeric(df.snp.pos@ranges@start[idx.snps])
    vec.snp.pos <- as.numeric(df.snp.pos@ranges@start)
    # vec.snp.bases <- as.character(unlist(df.snp.pos@elementMetadata$ALT)) # 
    
    # vec.snp.bases <- as.character(unlist(lapply(df.snp.pos@elementMetadata$ALT, function(m) m[1])))
    vec.ref.bases <- as.character(df.snp.pos@elementMetadata$REF)
    
    vec.snp.bases <- CharacterList(df.snp.pos@elementMetadata$ALT)
    vec.snp.bases = unstrsplit(vec.snp.bases, sep = ",")
    vec.snp.bases <- gsub("\\,.*", "", vec.snp.bases)
    
    #vec.snp.bases <- as.character(unlist(df.snp.pos@elementMetadata$ALT[idx.snps]))
    
    vec.snp.bases <- ifelse(vec.snp.bases == "<DEL>", "N",vec.snp.bases) # both represented separately
    vec.snp.bases <- ifelse(vec.snp.bases == "<INS>", "N",vec.snp.bases)
    
    # remove empty strings # 
    idx.snp.exceptions <- which(vec.snp.bases == "")
    
    if(length(idx.snp.exceptions) > 0){
      vec.snp.pos <- vec.snp.pos[-idx.snp.exceptions]
      vec.snp.bases <- vec.snp.bases[-idx.snp.exceptions]
      vec.ref.bases <- vec.ref.bases[-idx.snp.exceptions]
    }
    
    df.snp_positions.i <- data.frame(ID = paste(i ,"_", vec.snp.pos, sep ="") , chromosome = rep(paste(i ,"_maternal", sep ="") , length(vec.snp.pos)) , 
                                     position = vec.snp.pos, strand =  rep(1, length(vec.snp.pos)),
                                     Ref_SNP = paste(vec.ref.bases ,"/", vec.snp.bases, sep =""))
    
    df.snp_positions <- rbind(df.snp_positions, df.snp_positions.i)
    
    genome.reference <- DNAString(genome[[i]])
    genome.mutant    <- replaceLetterAt(genome.reference, vec.snp.pos, vec.snp.bases) # vergleich vor austausch
    
    l.genome.mutant[[i]] <- genome.mutant
    
    #  v.sequences <- unlist(lapply(l.sequences, function(m) {m[[1]]}))
    # if(FALSE){
    #   Sequences = DNAStringSet(genome.mutant)
    #   #Sequences <- read.DNAStringSet(FastaFile, "fasta")
    #   names(Sequences) <- as.character(paste("chr",i,sep = ""))
    #   writeXStringSet(Sequences, paste("output/mo17Genome/chr",i,".fasta", sep = ""), format="fasta")
    # }
    
  }
  
  write.table(df.snp_positions,  paste("data/df.snp_positions.csv", sep = ""), col.names = FALSE, row.names = FALSE, sep ="\t")
  saveRDS(l.genome.mutant, "tmp/genome_mutant.rds")
  
}else{
  df.snp_positions <- read.table(paste("data/df.snp_positions.csv", sep = ""),  sep ="\t")
  l.genome.mutant <- readRDS("tmp/genome_mutant.rds")
}




# binding peaks
message("loading binding peak datasets...")

if(!b.load_filtered_binding_peaks){
  
  strt<-Sys.time()
  df.peaks <- read.table(path.bindingPeaks, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
  df.peaks["seqnames"] <- as.numeric(unlist(lapply(strsplit(df.peaks$Position, ":"), function(m) m[[1]])))
  df.peaks["posPeak"] <- as.numeric(unlist(lapply(strsplit(df.peaks$Position, ":"), function(m) m[[2]])))
  df.peaks["start"] <- df.peaks$posPeak - s.half_window_size
  df.peaks["end"] <- df.peaks$posPeak + s.half_window_size
  
  l.df.peaks.all <- vector(mode = "list", length = length(l.path.bindingPeaksAll))
  for(i in 1:length(l.path.bindingPeaksAll)){
    l.df.peaks.all[[i]] <- read.table(l.path.bindingPeaksAll[[i]], header = F, sep = "\t", stringsAsFactors = FALSE, fill = TRUE, quote = "")  
    l.df.peaks.all[[i]] <- l.df.peaks.all[[i]][2:nrow(l.df.peaks.all[[i]]),]
  }
  
  # 6 von 6 
  l.df.peaks.all[[1]]["left"] <- l.df.peaks.all[[1]]$V2 - s.half_window_size.input
  l.df.peaks.all[[1]]["right"] <- l.df.peaks.all[[1]]$V3 + s.half_window_size.input
  
  if(FALSE){
    strt<-Sys.time() 
    #cl<-makeCluster(min(n.chromosomes, n.cpus))
    #registerDoParallel(cl)
    # l.peaks.filtered <-  foreach(i = 1:n.chromosomes, .packages=c("seqinr", "VariantAnnotation", "Biostrings")) %dopar% { 
    df.peaks.filtered <- c()
    for(i in 1:n.chromosomes){
      l.peaks.chromosome <- vector(mode = "list", length = 6)
      l.idx.blacklist <- vector(mode = "list", length = 6)
      for(l in 1:length(l.path.bindingPeaksAll)){
        l.peaks.chromosome[[l]] <- subset(l.df.peaks.all[[l]], l.df.peaks.all[[l]]$V1 == paste("chr", i, sep = ""))
        l.idx.blacklist[[l]] <- numeric(nrow(l.peaks.chromosome[[l]]))
      }
      v.idx.keep <- numeric(nrow(l.peaks.chromosome[[1]]))
      pb <- txtProgressBar(min = 0, max = nrow(l.peaks.chromosome[[1]]), style = 3)
      for(j in 1:nrow(l.peaks.chromosome[[1]])){
        setTxtProgressBar(pb, j)
        n.found <- 0
        for(l in 2:length(l.peaks.chromosome)){
          b.found = FALSE
          for(k in 1:nrow(l.peaks.chromosome[[l]])){
            
            if(!k %in% l.idx.blacklist[[l]]){        
              
              if( l.peaks.chromosome[[1]]$V1[j] == l.peaks.chromosome[[l]]$V1[k]){
                # give me a middle 
                if(l.peaks.chromosome[[l]]$V3[k] - 100 > l.peaks.chromosome[[1]]$left[j] & l.peaks.chromosome[[l]]$V3[k] - 100 < l.peaks.chromosome[[1]]$right[j]){
                  n.found = n.found + 1
                  b.found = TRUE
                  l.idx.blacklist[[l]][k] <- k
                  break
                }
              }
            }
          }
          if(b.found == FALSE){
            break
          }
        }
        if(n.found == 5){
          v.idx.keep[j] <- j  
        }
      }
      close(pb)
      
      v.idx.keep <- v.idx.keep[v.idx.keep != 0]
      l.peaks.chromosome[[1]] <- l.peaks.chromosome[[1]][v.idx.keep,]
      
      df.peaks.filtered <- rbind(df.peaks.filtered, l.peaks.chromosome[[1]])
      # l.peaks.chromosome[[1]]
    }
    #stopCluster(cl)
    print(Sys.time()-strt)
    
  }else{
    
    l.peaks.filtered <- readRDS("tmp/l.peaks.filtered.rds")  
    
    df.peaks.filtered <- c()
    for(i in 1:n.chromosomes){
      df.peaks.filtered <- rbind(df.peaks.filtered, l.peaks.filtered[[i]])
    }
    # 
    df.peaks.final <- c()
    for(i in 1:n.chromosomes){
      
      df.peaks.chromosome <- subset(df.peaks, df.peaks$seqnames == i)
      df.peaks.filtered.chromosome <- subset(df.peaks.filtered, df.peaks.filtered$V1 == paste("chr", i, sep = ""))
      
      v.idx.keep <- numeric(nrow(df.peaks.chromosome))
      pb <- txtProgressBar(min = 0, max = nrow(df.peaks.chromosome), style = 3)
      for(j in 1:nrow(df.peaks.chromosome)){
        setTxtProgressBar(pb, j)
        for(k in 1:nrow(df.peaks.filtered.chromosome)){
          dist <- abs(df.peaks.chromosome$posPeak[j]  -  df.peaks.filtered.chromosome$V3[k] - 100)  
          if(dist < 500){
            v.idx.keep[j] <- j
            break
          }
        }
      }
      close(pb)
      
      v.idx.keep <- v.idx.keep[v.idx.keep>0]
      df.peaks.chromosome <- df.peaks.chromosome[v.idx.keep, ]
      df.peaks.final <- rbind(df.peaks.final, df.peaks.chromosome)
    }
    
    if(FALSE){
      saveRDS(df.peaks.final, "tmp/df.peaks.final.rds")  
    }
    
  }
}else{
  df.peaks.final <- readRDS("tmp/df.peaks.final.rds") # TABLE S5 
}

# further peak processing
if(TRUE){
  
  df.peaks <- df.peaks.final
  df.peaks["start"] <- df.peaks$posPeak - s.half_window_size
  df.peaks["end"] <- df.peaks$posPeak + s.half_window_size
  
  message("number of binding peaks: ", nrow(df.peaks))
  
  message("loading binding QTL datasets...")
  strt<-Sys.time()
  postTotal = read.csv(path.bQTLs, header=T, stringsAsFactors = FALSE)
  df.bQTLsInput = read.csv(path.bQTLsInput, header=T, stringsAsFactors = FALSE)
  print(Sys.time() - strt)
  #postTotal = read.table(pd, sep="\t", header=T, stringsAsFactors = FALSE)
  
  df.bQTLsInput <- subset(df.bQTLsInput, df.bQTLsInput$refCount > 0 & df.bQTLsInput$altCount > 0)
  df.bQTLsInput["id"] <- paste(df.bQTLsInput$contig, "_", df.bQTLsInput$position, sep = "")
  postTotal["id"] <- paste(postTotal$contig, "_", postTotal$position, sep = "")
  
  # store all snps 
  df.bQTLsOutputformat <- df.bQTLsInput[,c("id","contig", "position", "refAllele", "altAllele")]
  
  # input numbers
  n.bQTLs_readDepth50 <- nrow(postTotal)
  n.bQTLs_readDepth5 <- nrow(df.bQTLsInput)
  
  print(n.bQTLs_readDepth50)
  print(n.bQTLs_readDepth5)
  
  # re-insert ASBs with at least minReadsPerAllele reads in ref and alt (based on the input with presence in at least 5 in 6 repeats)
  
  # blacklist areas 
  minReadsPerAllele = minReadDepth * 0.1 
  postTotalSafe <- subset(postTotal, postTotal$refCount >= minReadsPerAllele & postTotal$altCount >= minReadsPerAllele)
  postTotalCheck <- subset(postTotal, !postTotal$id %in% postTotalSafe$id)
  postTotalCheck <- subset(postTotalCheck, postTotalCheck$id %in% df.bQTLsInput$id)
  postTotal <- rbind(postTotalSafe, postTotalCheck)
  
  # if(b.save_intermediate_results){
  #   write.csv(postTotal, "paper_tmp_files/bQTL_after_input.csv")
  # }
  # 
  
  # filtering
  n.bQTLs_readDepth50.both_alleles_10Percent_bias <- nrow(postTotal)
  n.bQTLs_readDepth5.both_alleles_nonZero_cound <- nrow(df.bQTLsInput)
  
  print(n.bQTLs_readDepth50.both_alleles_10Percent_bias)
  print(n.bQTLs_readDepth5.both_alleles_nonZero_cound)
  
  message("input bQTL after min reads per allele filter: ", nrow(postTotal))
  
  v.postFreq <- postTotal$refCount / postTotal$totalCount
  v.postFreq <- v.postFreq[which(v.postFreq != 0)]
  v.postFreq <- v.postFreq[which(v.postFreq != 1)]
  prob.bias <- median(v.postFreq)
  
  message("bias probability: ", prob.bias)
  
  #p = 0.5 # 0.551948 # expectancy
  p = 1 - prob.bias # 0.5120976 # post frequency - ref / all
  test <- sapply(1:nrow(postTotal), function(i) binom.test(as.integer(postTotal$altCount[i]), as.integer(postTotal$totalCount[i]), p = p)$p.value) #  analytic.pv(preProps[i],preVars[i],postProps[i],depths[i]))
  fdr <- p.adjust(test, "bonferroni")
  postTotal["p-value (corrected)"] <- fdr
  
  postTotal <- subset(postTotal, postTotal$totalCount >= minReadDepth)
  
  postTotal["POSTfreq"] <- postTotal$refCount / postTotal$totalCount
  postTotal["POSTallele"] <- ifelse(postTotal$altCount > postTotal$refCount, postTotal$altAllele, postTotal$refAllele)
  
  # postTotal <- subset(postTotal, postTotal$POSTfreq < 1 & postTotal$POSTfreq > 0)
  median(postTotal$POSTfreq)
  
  message("make figure per paper")
  df.bQTLsInput["POSTfreq"] <- df.bQTLsInput$refCount / df.bQTLsInput$totalCount
  # df.bQTLsInput <- subset(df.bQTLsInput, df.bQTLsInput$totalCount > 5)
  
  # perform chromosome based - blacklisting
  
  #message("input bQTL: ", nrow(postTotal))
  #bQTL_scatterplot(postTotal=postTotal)
  
  ## manual artifact removal
  message("remove manually labelled artifact regions - further reducing significant ")
  # make figure 
  # bQTL_scatterplot_chr(postTotal=postTotal, chr = 3)
  chr = 3
  cut.left <-  186240476
  cut.right <- 210080072
  print(abs(cut.left - cut.right) / 1e6)
  
  idx.cut <- which(df.bQTLsInput$contig == chr & df.bQTLsInput$position < cut.right & df.bQTLsInput$position > cut.left)
  df.bQTLsInput <- df.bQTLsInput[!seq(1:nrow(df.bQTLsInput)) %in%  idx.cut,]
  
  chr = 10
  cut.left <- 2965337
  cut.right <- 4124958
  print(abs(cut.left - cut.right) / 1e6)
  
  idx.cut <- which(df.bQTLsInput$contig == chr & df.bQTLsInput$position < cut.right & df.bQTLsInput$position > cut.left)
  df.bQTLsInput <- df.bQTLsInput[!seq(1:nrow(df.bQTLsInput)) %in%  idx.cut,]
  
  # double check
  # bQTL_scatterplot(postTotal=df.bQTLsInput, chr = 10)
  
  postTotal <- subset(postTotal, postTotal$id %in% df.bQTLsInput$id)
  postTotal <- postTotal[,!names(postTotal) %in% "id"]
  
  message("input bQTL after min reads per allele filter: ", nrow(postTotal))
  
  # if(b.save_intermediate_results){
  #   write.csv(postTotal, "paper_tmp_files/bQTL_after_input_blacklisting.csv")
  # }
}


l.postTotal <- vector(mode = "list", length = 2)
l.postTotal[[1]] <- subset(postTotal, postTotal$`p-value (corrected)` < th.p.bQTL) # significant bQTL
l.postTotal[[1]] <- subset(l.postTotal[[1]], l.postTotal[[1]]$POSTfreq < 1 & l.postTotal[[1]]$POSTfreq > 0)


if(b.save_intermediate_results){
  write.table(l.postTotal[[1]], paste(folder_output, "/ASBs.txt", sep = ""),  row.names = FALSE, quote = FALSE, sep ="\t")
}

message("bQTL (significance p < ", th.p.bQTL, ", ", nrow(l.postTotal[[1]]), ")")
message("bQTLs after linkage ", nrow(l.postTotal[[1]]))
hist(l.postTotal[[1]]$POSTfreq, breaks = 100, main = "Allelic bias of bQTL", xlab = "allelic bias")


if(b.peak_analysis){
  
  message("Estimate significant bQTLs in binding peak ranges...")
  
  l.selection <- vector(mode = "list", length = 2)
  l.sigSnps_in_peaks <- vector(mode = "list", length = 2)
  
  df.non_ASB_peaks <- c()
  
  for(s in 1:1){
    
    n.total <- 0
    
    postTotal.significant <- l.postTotal[[s]]
    
    strt<-Sys.time() 
    cl<-makeCluster(min(n.chromosomes, n.cpus))
    registerDoParallel(cl)
    
    l.res <- foreach(i = 1:n.chromosomes, .packages = c("seqinr", "Biostrings", "VariantAnnotation")) %dopar% {     
      
      sigSnps_in_peaks <- 0
      df.selection <- c()
      
      df.snp.pos <- readRDS(paste("tmp/df.snp.pos_",i, ".rds"))
      
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
      
      n.total <- n.total + length(vec.snp.pos)
      
      # remove empty strings # 
      idx.snp.exceptions <- which(vec.snp.bases == "")
      
      if(length(idx.snp.exceptions) > 0){
        vec.snp.pos <- vec.snp.pos[-idx.snp.exceptions]
        vec.snp.bases <- vec.snp.bases[-idx.snp.exceptions]
      }
      
      genome.reference <- DNAString(genome[[i]])
      genome.mutant    <- replaceLetterAt(genome.reference, vec.snp.pos, vec.snp.bases) # vergleich vor austausch
      
      # binding peaks 
      if(s == 1){
        tf_target_bind.sset <- subset(df.peaks, df.peaks$seqnames == i)
      }else if(s == 2){
        tf_target_bind.sset <- subset(df.non_ASB_peaks, df.non_ASB_peaks$seqnames == i)
      }
      
      peak <- tf_target_bind.sset
      
      #   for(p in 1:nrow(peak)){ # double check for all
      #     idx.snps <- which(vec.snp.pos >= peak$start[p] & vec.snp.pos <= peak$end[p])
      #     snp_in_peaks <- snp_in_peaks + length(idx.snps)  
      #   }
      
      postTotal.significant.i <- subset(postTotal.significant ,postTotal.significant$contig == i)
      
      idx.non_asb_peaks <- numeric()
      
      pb <- txtProgressBar(min = 0, max = nrow(peak), style = 3)
      
      for(p in 1:nrow(peak)){ # double check for all
        
        setTxtProgressBar(pb, p)
        
        idx.snps <- which(postTotal.significant.i$position >= peak$start[p] & postTotal.significant.i$position <= peak$end[p])
        
        if(length(idx.snps) > 0){
          df.selection <- rbind(df.selection ,postTotal.significant.i[idx.snps,])
          sigSnps_in_peaks <- sigSnps_in_peaks + length(idx.snps)  
        }else{
          idx.non_asb_peaks <- c(idx.non_asb_peaks, p)
        }
        
      }
      
      close(pb)
      
      l.res <- list(df.selection = df.selection, sigSnps_in_peaks = sigSnps_in_peaks, df.non_ASB_peaks = peak[idx.non_asb_peaks,])
      
    }
    
    stopCluster(cl)
    print(Sys.time()-strt)
    
    l.selection[[s]] <- (l.res[[1]]$df.selection)
    l.sigSnps_in_peaks[[s]] <- (l.res[[1]]$sigSnps_in_peak)
    
    df.non_ASB_peaks <- l.res[[1]]$df.non_ASB_peaks
    
    for(i in 2:n.chromosomes){
      l.selection[[s]] <- rbind(l.selection[[s]], (l.res[[i]]$df.selection))
      l.sigSnps_in_peaks[[s]] <- l.sigSnps_in_peaks[[s]] + (l.res[[i]]$sigSnps_in_peak)
      df.non_ASB_peaks <- rbind(df.non_ASB_peaks, l.res[[i]]$df.non_ASB_peaks)
    }
    
  }
  
  message("...finished")
  message("ASBs" , nrow(l.selection[[1]]))
  
  if(b.save_intermediate_results){
    write.table(l.selection[[1]], paste(folder_output, "/bQTL_after_input_blacklisting_pvalue_filter_in_peaks.txt", sep = ""), quote = FALSE, row.names = FALSE, sep ="\t")
  }
  
  # message("Foldchange of strongly bound bQTLs in peaks vs random snps, FC:", fc, " / p-value: ", p.val)
  
}
n_snps_peak_filtered = nrow(l.selection[[1]])

if(b.disequilibriumLinkageHandling){
  
  # run linkage 
  message("Processing linkage disequilibrium of significant bQTLs...")
  
  postTotal.significant <- l.selection[[1]] # 
  
  # initial linkage disequilibrium filtering 
  postTotal.significant.tmp <- postTotal.significant
  postTotal.significant <- postTotal.significant[-(1:nrow(postTotal.significant)),]
  
  n.flips <- 0
  
  strt<-Sys.time() 
  cl<-makeCluster(min(n.chromosomes, n.cpus))
  registerDoParallel(cl)
  
  l.postTotal.significant <- foreach(i = 1:n.chromosomes) %dopar% {     
    
    postTotal.significant.i <- subset(postTotal.significant.tmp, postTotal.significant.tmp$contig == i)
    postTotal.significant.i <- postTotal.significant.i[order(postTotal.significant.i$position),]
    
    finished <- FALSE
    df.keep <- postTotal.significant[-(1:nrow(postTotal.significant)),]
    
    while(!finished){
      
      # estimate distances of neighboring SNPs
      dist <- postTotal.significant.i$position[2:nrow(postTotal.significant.i)] - postTotal.significant.i$position[1:(nrow(postTotal.significant.i) - 1)]
      
      if(!any(dist < s.disequilibriumDistance)){
        finished <- TRUE
      }
      
      for(j in 1:length(dist)){
        
        if(dist[j] < s.disequilibriumDistance){
          
          # pairwise evaluation
          set <- postTotal.significant.i[c(j,(j+1)),]
          
          # using p-values to guide selection 
          if(s == 1) # in case of significant
            i.max <- which(set$`p-value (corrected)` == max(set$`p-value (corrected)`))
          else
            i.max <- which(set$`p-value (corrected)` == min(set$`p-value (corrected)`))
          
          # different directions - keep both
          if(set$POSTfreq[1] < 0.5 & set$POSTfreq[2] > 0.5 | set$POSTfreq[1] > 0.5 & set$POSTfreq[2] < 0.5){
            
            # keep both
            n.flips <- n.flips + 1
            df.keep <- rbind(df.keep, set)
            
          }
          
          # if s = 1 : remove SNP within proximity with higher p-value
          idx <- ifelse(i.max == 1, j, j + 1)
          postTotal.significant.i <- postTotal.significant.i[-idx,]
          
          break
          
        }
      }
    }
    
    postTotal.significant.i <- rbind(postTotal.significant.i, df.keep)
    postTotal.significant.i <- unique(postTotal.significant.i)
    
    # postTotal.significant <- rbind(postTotal.significant, postTotal.significant.i)
    postTotal.significant.i
    
  }
  
  stopCluster(cl)
  print(Sys.time()-strt)
  
  # filtered significant bQTLs
  postTotal.significant <- c()
  for(i in 1:n.chromosomes){
    postTotal.significant <- rbind(postTotal.significant, l.postTotal.significant[[i]])
  }
  
  l.selection[[1]] <- postTotal.significant
  
  hist(l.selection[[1]]$POSTfreq, breaks = 100) # histogram of significant peaks
  
  if(b.save_intermediate_results){
    write.table(l.selection[[1]], paste(folder_output, "/bQTL_after_input_blacklisting_pvalue_filter_in_peaks_after_LinkageDisequilibrium.txt", sep = ""),  row.names = FALSE, quote = FALSE, sep ="\t")
  }
}
n_snps_linkage_filtered = nrow(l.selection[[1]])

bQTL_scatterplot(postTotal=l.selection[[1]]) # fig 2 c



message("create population")

if(!b.load_snp_to_gene_partitioning_and_background_sampling){
  
  source("HashSeq_populationSampling.R")
  l.postTotal = sample_nonsignificant_background_snps_in_peaks(l.postTotal)
  
  
  # run gene partitioning 
  # run the background and gene partitioning 
  source("HachSeq_partitioning.R")
  
  # l.bQTL_gene_partitioning <- readRDS(paste("tmp/l.bQTL_gene_partitioning_", timeStamp, ".rds", sep = ""))
  l.bQTL_gene_partitioning = perform_SNP_to_gene_partitioning(df.gene_annotation=df.gene_annotation,
                                                              v.genePartitions=v.genePartitions)
  
  
  
  l.res <- create_background_distribution_with_similar_gene_partitioning(l.bQTL_gene_partitioning = l.bQTL_gene_partitioning, multiplyer = s.multiplyer, 
                                                                         v.partitions = v.partitions, b.duplicateRemoval = b.duplicateRemoval, seed.randomGenerator = 1234)
  
  l.bQTL_gene_partitioning.sampling <- l.res$l.bQTL_gene_partitioning
  df.distribution <- l.res$df.distribution
  
  saveRDS(l.res,"tmp/l.bQTL_gene_partitioning_withGeneDistances_backgroundSampled.rds") 
  
  if(b.save_intermediate_results){
    
    #saveRDS(l.bQTL_gene_partitioning.sampling, paste("tmp/l.bQTL_gene_partitioning_no_duplicates_bg_sampled.", timeStamp, ".rds", sep = ""))  
    write.csv(df.distribution, paste("output/df.ASB_and_bgSnp_gene_partitions.csv", sep = ""))  
    write.csv(l.bQTL_gene_partitioning.sampling[[1]], paste("output/df.ASB_gene_partitioning.csv", sep = ""))  
    write.csv(l.bQTL_gene_partitioning.sampling[[2]], paste("output/df.bgSNP_gene_partitioning.csv", sep = ""))  
    #df.bQTL_gene_partitioning_bg_equivalent <- read.csv("output/df.bQTL_gene_partitioning_bg_equivalent.csv")
  }
  
}else{
  
  l.bQTL_gene_partitioning <- readRDS(paste("tmp/l.bQTL_gene_partitioning_withGeneDistances_backgroundSampled.rds", sep = ""))$l.bQTL_gene_partitioning
  
  write.table(l.bQTL_gene_partitioning[[1]], paste(folder_output, "manuscript/supplement/S8.txt", sep = ""),  row.names = FALSE, quote = FALSE, sep ="\t")
  
  # distribution of 3297 ASBs  - Figure 2 d
  v.distribution <- apply(l.bQTL_gene_partitioning[[1]][,v.partitions], 2, table)["yes",] 
   
  
}

message("motif analysis")

source("HachSeq_motifs.R")
predefined_motif_analysis(l.bQTL_gene_partitioning, motifs)


l.motif_analysis <- readRDS(paste("tmp/l.motif_analysis.rds", sep = ""))
l.nucleotideInPeaks <- readRDS(paste("tmp/l.nucleotideInPeaks.rds", sep = ""))

l.motif_analysis.postprocessed = readRDS(paste("tmp/l.motif_analysis.postprocessed.rds", sep = ""))



# GWAS 



message("17463 high confidence ZmBZR1 binding peaks") # TABLE S2
if(FALSE){
  
  l.df.ChipSeqB73 <- vector(mode = "list", length = length(l.path.ChipSeqB73))
  for(i in 1:length(l.path.ChipSeqB73)){
    l.df.ChipSeqB73[[i]] <- read.table(l.path.ChipSeqB73[[i]], header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE, quote = "")  
    l.df.ChipSeqB73[[i]]["chr"] <- as.numeric(gsub("\\:.*", "", l.df.ChipSeqB73[[i]]$Position))
    l.df.ChipSeqB73[[i]]["pos"] <- as.numeric(gsub(".*:", "", l.df.ChipSeqB73[[i]]$Position))
    l.df.ChipSeqB73[[i]]["left"] <- l.df.ChipSeqB73[[i]]$pos - s.half_window_size.input
    l.df.ChipSeqB73[[i]]["right"] <- l.df.ChipSeqB73[[i]]$pos + s.half_window_size.input
  }
  
  strt<-Sys.time() 
  #cl<-makeCluster(min(n.chromosomes, n.cpus))
  #registerDoParallel(cl)
  # l.peaks.filtered <-  foreach(i = 1:n.chromosomes, .packages=c("seqinr", "VariantAnnotation", "Biostrings")) %dopar% { 
  df.ChipSeq.filtered <- c()
  for(i in 1:n.chromosomes){
    l.peaks.chromosome <- vector(mode = "list", length = 3)
    l.idx.blacklist <- vector(mode = "list", length = 3)
    for(l in 1:length(l.df.ChipSeqB73)){
      l.peaks.chromosome[[l]] <- subset(l.df.ChipSeqB73[[l]], l.df.ChipSeqB73[[l]]$chr == i)
      l.idx.blacklist[[l]] <- numeric(nrow(l.peaks.chromosome[[l]]))
    }
    v.idx.keep <- numeric(nrow(l.peaks.chromosome[[1]]))
    pb <- txtProgressBar(min = 0, max = nrow(l.peaks.chromosome[[1]]), style = 3)
    for(j in 1:nrow(l.peaks.chromosome[[1]])){
      setTxtProgressBar(pb, j)
      n.found <- 0
      for(l in 2:length(l.peaks.chromosome)){
        b.found = FALSE
        for(k in 1:nrow(l.peaks.chromosome[[l]])){
          if(!k %in% l.idx.blacklist[[l]]){        
            if(l.peaks.chromosome[[1]]$chr[j] == l.peaks.chromosome[[l]]$chr[k]){
              if(l.peaks.chromosome[[l]]$pos[k] > l.peaks.chromosome[[1]]$left[j] & l.peaks.chromosome[[l]]$pos[k] < l.peaks.chromosome[[1]]$right[j]){
                n.found = n.found + 1
                b.found = TRUE
                l.idx.blacklist[[l]][k] <- k
                break
              }
            }
          }
        }
        if(b.found == FALSE){
          break
        }
      }
      if(n.found == 2){
        v.idx.keep[j] <- j  
      }
    }
    close(pb)
    v.idx.keep <- v.idx.keep[v.idx.keep != 0]
    l.peaks.chromosome[[1]] <- l.peaks.chromosome[[1]][v.idx.keep,]
    df.ChipSeq.filtered <- rbind(df.ChipSeq.filtered, l.peaks.chromosome[[1]])
    # l.peaks.chromosome[[1]]
  }
  #stopCluster(cl)
  print(Sys.time()-strt)

  saveRDS(df.ChipSeq.filtered, "tmp/df.ChipSeq.filtered.rds")
}else{
  df.ChipSeq.filtered <- readRDS("tmp/df.ChipSeq.filtered.rds")
}




message("add gene features")
### different input datasets 

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




#### 
source("utils/genomic_location.R")

test = add_genomic_location(postTotal.significant, 
                                 pos_peak,
                                 chr_Nr,
                                 df.gene_annotation,
                                 v.genePartitions = c("gene", "five_prime_UTR",  "CDS", "three_prime_UTR", "exon"),
                                 n.cpus = 1)

# S4 - genes of 17431 high confidence ZmBZR1 binding peaks

# check for BR regulated (up and down)
length(unique(test$gene.ID)) # 6371
v.distribution <- apply(test[,v.partitions], 2, table)["yes",] 
v.distribution / sum(v.distribution)

length(intersect(df.rnaseq.down_regulated$X, unique(test$gene.ID))) # 580
length(intersect(df.rnaseq.up_regulated$X, unique(test$gene.ID))) # 469
                                 

df.bQTL_RNAseq


                                 
# else{
#   message("no denovo motif and in peak analysis")
#   if(b.BRRE_GBOX.only){
#     l.motif_analysis <- readRDS(paste("tmp/l.motif_analysis.BRRE_GBOX_only_",timeStamp, ".rds", sep = ""))
#     l.nucleotideInPeaks <- readRDS(paste("tmp/l.nucleotideInPeaks.BRRE_GBOX_only_",timeStamp, ".rds", sep = ""))
#   }else{
#     l.motif_analysis <- readRDS(paste("tmp/l.motif_analysis_",timeStamp, ".rds", sep = ""))
#     l.nucleotideInPeaks <- readRDS(paste("tmp/l.nucleotideInPeaks_",timeStamp, ".rds", sep = ""))
#   }
# }

write.csv(df.bQTL_gene_partitioning, "output/df.bQTL_gene_partitioning_peaks.csv", row.names = FALSE)


### orthologs ... 

df.bQTL_gene_partitioning.subset <- subset(df.bQTL_gene_partitioning, df.bQTL_gene_partitioning$nonGenic == "no")
length(unique(df.bQTL_gene_partitioning.subset$gene.ID))

length(unique(df.bQTL_gene_partitioning.subset$Arabidopsis_ortholog))


df.ChipSeq.gene_partitioning.subset <- subset(df.ChipSeq.gene_partitioning, !is.na(df.ChipSeq.gene_partitioning$Arabidopsis_ortholog))
df.ChipSeq.gene_partitioning.subset$Arabidopsis_ortholog <- gsub("\\..*","", df.ChipSeq.gene_partitioning.subset$Arabidopsis_ortholog)

# v.genePartitions = c("promoter_5kb", "promoter_1kb","gene", "five_prime_UTR", "intron", "CDS", "three_prime_UTR", "exon", "post_gene_1kb")
table(df.ChipSeq_gene_partitioning$non_genic)
test = subset(df.ChipSeq_gene_partitioning, df.ChipSeq_gene_partitioning$non_genic == "no")
length(unique(test$gene.ID))

test2 = add_features_to_gene_partitioning(test)
length(unique(test2$Arabidopsis_ortholog))

test2 = test2[,c("gene.ID", "Arabidopsis_ortholog")]
test3 = subset(test2, !is.na(test2$Arabidopsis_ortholog))

write.table(test3, "S4.txt", sep ="\t", row.names = F, quote = F)

test = read.table("manuscript/supplement/Table_S4.txt", sep ="\t", header = T)






message("Figure 2 - Relative comparative partition analysis ")

df.partitions["percentage_significant"] <- df.partitions[,2] / df.partitions[10,2]
df.partitions["percentage_non_significant"] <- df.partitions[,3] / df.partitions[10,3]

write.csv(df.partitions, paste("output/df.partition_plus_", timeStamp, ".csv", sep = ""))

