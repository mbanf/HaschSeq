
df.ASBs_methylation <- readRDS("df.ASBs_methylation.rds")

df.ASBs = l.bQTL_gene_partitioning[[1]]
df.bpSNPs = l.bQTL_gene_partitioning[[2]]


# TODO: get a background selection on the motif by chance occurrence ... 
asbs_explaining_motifs <- function(df.ASBs, genome, l.genome.mutant, offset = 5, core_motif = "CGTG", n.cpus = 2){
  
  brre.fw = c("CGTGCG", "CGTGTG", # original motifs
              "CGTGAG", "CGTGGG", 
              "CGTGCA", "CGTGCC", "CGTGCT", 
              "CGTGTA", "CGTGTC", "CGTGTT")
  
  brre.bw = c("CGCACG","CACACG", # original motifs
              "CCCACG", "CTCACG",
              "AGCACG", "GGCACG", "TGCACG", 
              "AACACG", "GACACG", "TACACG")
  
  ebox.combinations = c("CACGTG") 
  
  # ebox.combinations = c("CAAATG", "CAACTG", "CAAGTG", "CAATTG",
  #                       "CACATG", "CACCTG", "CACGTG", "CACTTG",
  #                       "CAGATG", "CAGCTG", "CAGGTG", "CAGTTG",
  #                       "CATATG", "CATCTG", "CATGTG", "CATTTG")

  motifs = unique(c(brre.fw, brre.bw, ebox.combinations))
                        
  message("Performing predefined motif analysis...")
  
  # allow for double motifs in ASBs
  postTotal.significant <- df.ASBs #  df.bpSNPs #df.ASBs # df.bpSNPs #  df.ASBs # df.bpSNPs #  df.ASBs
  
  strt<-Sys.time() 
  
  df.motif_analysis <- postTotal.significant[-(1:nrow(postTotal.significant)),]
  df.motif_analysis["motif.ref"] <- c()
  df.motif_analysis["motif.mutant"] <- c()

  # most basic way to estimate the numbers => no dataset generation 
  
  
  strt <- Sys.time()
  cl <- makeCluster(n.cpus)
  registerDoParallel(cl)
  
  l.SNPs.motifs <- foreach(i = 1:n.chromosomes, .packages=c("seqinr", "VariantAnnotation", "Biostrings")) %dopar% {

      print(paste("processing chromosome ", i))
      
      genome.reference <- DNAString(genome[[i]])
      genome.mutant    <- l.genome.mutant[[i]]
      
      postTotal.significant.i <- subset(postTotal.significant, postTotal.significant$contig == i) 
      postTotal.significant.i["has_motif"] = FALSE
      pb <- txtProgressBar(min = 0, max = nrow(postTotal.significant.i), style = 3)
      for(k in 1:nrow(postTotal.significant.i)){
        setTxtProgressBar(pb, k)
      
        pos = postTotal.significant.i$position[k]
        offset = 5
        seq.reference = subseq(genome.reference, pos - offset, pos + offset)
        seq.mutant    = subseq(genome.mutant, pos - offset, pos + offset)
        
        # perfect hits
        n.total = 0
        for(m in 1:length(motifs)){ # DNAStringSet(matches) => check if these have the core motif
          n.reference = length(start(matchPattern(motifs[m], seq.reference, fixed = TRUE, max.mismatch = 0)))
          n.mutant = length(start(matchPattern(motifs[m], seq.mutant, fixed = TRUE, max.mismatch = 0)))
          n.total = n.total + n.reference + n.mutant
        }

        if(n.total > 0){
          postTotal.significant.i[k, "has_motif"] <- TRUE
        }
        
      }
      close(pb)
    postTotal.significant.i
  }
  stopCluster(cl)
  print(Sys.time() - strt)  
      
    
  df.ASBs_w_motifs <- c()
  for(i in 1:n.chromosomes){
    df.ASBs_w_motifs <- rbind(df.ASBs_w_motifs, l.SNPs.motifs[[i]])
  }
  
  df.ASBs_methylation[,"has_motif"] <- df.ASBs_w_motifs$has_motif
  df.ASBs_methylation[,"overlap"] <- apply(df.ASBs_methylation[,c("is_methylated", "has_motif")], 1, all)
  
  
  table(df.ASBs_methylation$overlap)
  
  test <- subset(df.ASBs_w_motifs, df.ASBs_w_motifs$POSTfreq >= 0.85 | df.ASBs_w_motifs$POSTfreq <= 0.15)
  table(df.ASBs_w_motifs$has_motif)
      table(test$has_motif)
      
      
      
      
    
      df.motif_analysis <- postTotal.significant[-(1:nrow(postTotal.significant)),]
      
      df.motif_analysis["motif.ref"] <- c()
      df.motif_analysis["motif.mutant"] <- c()
      df.motif_analysis["unique"] <- c()
      df.motif_analysis["isMotifInPeak"] <- c()
      df.motif_analysis["pos.motif"] <- c()
      
      
      
      for(m in 1:length(motifs)){
        
        
        
        
        # returns the starting position of the found motif
        v.reference <- as.numeric(start(matchPattern(motifs[m], genome.reference, fixed = TRUE)))
        v.mutant <- as.numeric(start(matchPattern(motifs[m], genome.mutant, fixed = TRUE)))
        
        vec.intersect <- intersect(v.reference, v.mutant)
        vec.reference.unique <- v.reference[!v.reference %in% v.mutant] # B73
        vec.mutant.unique <- v.mutant[!v.mutant %in% v.reference] # MO17
        
        #postTotal.significant.i[which(postTotal.significant.i$position %in% vec.mutant.unique),
        
        ## analyze peak based surroundings to 
        for(j in 1:length(vec.intersect)){
          
          #cat("Processing... ", round(j/length(vec.intersect) * 100, digits = 2) , "%", "\r"); flush.console()  
          
          # is the motif (not just the ASB) in peak
          idx.j <- which(tf_target_bind.sset$start <= vec.intersect[j] & tf_target_bind.sset$end >= vec.intersect[j] + s.motif_offset) # found a mutual motif in the peak -> take this peak
          
          #idx.j <- which(postTotal.significant.i$position >= vec.intersect[j] & postTotal.significant.i$position <= vec.intersect[j] + 5) # 
          
          if(length(idx.j) > 0){
            
            peak <- tf_target_bind.sset[idx.j,] 
            
            for(p in 1:nrow(peak)){ # double check for all
              
              idx.snps <- which(postTotal.significant.i$position >= peak$start[p] & postTotal.significant.i$position <= peak$end[p])
              
              if(length(idx.snps) > 0){
                
                for(k in 1:length(idx.snps)){
                  if(postTotal.significant.i$position[idx.snps[k]] - vec.intersect[j] < 0){
                    dist <- abs(postTotal.significant.i$position[idx.snps[k]] - vec.intersect[j])
                  }else{
                    dist <- abs(postTotal.significant.i$position[idx.snps[k]] - (vec.intersect[j] + s.motif_offset))
                  }
                  
                  # break
                  tb.nucleotides <- table(postTotal.significant.i[idx.snps[k],]$POSTallele)
                  
                  newrow <- data.frame(A = 0, G = 0, C = 0, T = 0, motif = motifs[m])
                  newrow[,names(tb.nucleotides)] <- dist #as.numeric(tb.nucleotides)
                  
                  df.nucleotideInPeaks <- rbind(df.nucleotideInPeaks, newrow)
                  
                }
              }
            }
          }
        }
        
        # analyze effect unique motifs
        for(j in 1:length(vec.reference.unique)){
          
          cat("Processing... ", round(j/length(vec.reference.unique) * 100, digits = 2) , "%", "\r"); flush.console()  
          
          idx.j <- which(postTotal.significant.i$position >= vec.reference.unique[j] & postTotal.significant.i$position <= vec.reference.unique[j] + s.motif_offset) 
          
          if(length(idx.j) > 0){
            
            v.seq.mutant <- subseq(genome.mutant, start=vec.reference.unique[j], end=vec.reference.unique[j] + s.motif_offset)  
            v.seq.reference <- subseq(genome.reference, start=vec.reference.unique[j], end=vec.reference.unique[j] + s.motif_offset)  
            
            newrow <- postTotal.significant.i[idx.j,] 
            
            newrow["motif.ref"] <- as.character(v.seq.reference)
            newrow["motif.mutant"] <- as.character(v.seq.mutant)
            newrow["unique"] <- "reference"
            
            newrow["pos.motif"] <- vec.reference.unique[j]
            
            in.peak <- any(tf_target_bind.sset$start <= vec.reference.unique[j] & tf_target_bind.sset$end >= vec.reference.unique[j] + s.motif_offset) # found a mutual motif in the peak -> take this peak
            
            if(in.peak)
              newrow["isMotifInPeak"] <- "yes"
            else
              newrow["isMotifInPeak"] <- "no"
            
            df.motif_analysis <- rbind(df.motif_analysis, newrow)
            
          }
          
        }
        
        for(j in 1:length(vec.mutant.unique)){
          
          cat("Processing... ", round(j/length(vec.mutant.unique) * 100, digits = 2) , "%", "\r"); flush.console()  
          
          idx.j <- which(postTotal.significant.i$position >= vec.mutant.unique[j] & postTotal.significant.i$position <= vec.mutant.unique[j] + s.motif_offset)  
          
          if(length(idx.j) > 0){
            
            # this is the original motif
            v.seq.mutant <- subseq(genome.mutant, start=vec.mutant.unique[j], end=vec.mutant.unique[j] + s.motif_offset)  
            
            # this is a position corresponding de novo motif
            v.seq.reference <- subseq(genome.reference, start=vec.mutant.unique[j], end=vec.mutant.unique[j] + s.motif_offset)  
            
            newrow <- postTotal.significant.i[idx.j,] 
            
            newrow["motif.ref"] <- as.character(v.seq.reference)
            newrow["motif.mutant"] <- as.character(v.seq.mutant)
            newrow["unique"] <- "mutant"
            
            newrow["pos.motif"] <- vec.mutant.unique[j]
            
            # motifs (not snps) in peak
            in.peak <- any(tf_target_bind.sset$start <= vec.mutant.unique[j] & tf_target_bind.sset$end >= vec.mutant.unique[j] + s.motif_offset) # found a mutual motif in the peak -> take this peak
            
            if(in.peak)
              newrow["isMotifInPeak"] <- "yes"
            else
              newrow["isMotifInPeak"] <- "no"
            
            df.motif_analysis <- rbind(df.motif_analysis, newrow) # use these motifs 
            
          }
        }
      }
      
      l.sets <- vector(mode = "list", length = 2)
      
      l.sets[[1]] <- df.nucleotideInPeaks
      l.sets[[2]] <- df.motif_analysis
      
      l.sets
    }
    
  }
  print(Sys.time()-strt)
  
  df.nucleotideInPeaks <- c()
  df.motif_analysis <- c()
  
  for(i in 1:n.chromosomes){
    df.nucleotideInPeaks <- rbind(df.nucleotideInPeaks, l.sets[[i]][[1]])
    df.motif_analysis <- rbind(df.motif_analysis, l.sets[[i]][[2]])
  }
  
  l.motif_analysis[[s]] <- df.motif_analysis
  l.nucleotideInPeaks[[s]] <- df.nucleotideInPeaks
  
  
  
  saveRDS(l.motif_analysis, paste("tmp/l.motif_analysis.rds", sep = ""))
  saveRDS(l.nucleotideInPeaks, paste("tmp/l.nucleotideInPeaks.rds", sep = ""))
  
  write.table(l.motif_analysis[[1]], paste(folder_output, "/motifs/S0_1.txt", sep = ""),  row.names = FALSE, quote = FALSE, sep ="\t")
  
  message("...finished")
  
}