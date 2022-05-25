

predefined_motif_analysis <- function(l.bQTL_gene_partitioning, motifs){
  
  message("Performing predefined motif analysis...")
  
  l.motif_analysis <- vector(mode = "list", length = 2)
  l.nucleotideInPeaks <- vector(mode = "list", length = 2)
  
  for(s in 1:2){ # significant and non-significant
    
    # work only on asb partitioning data (remove duplicates)
    postTotal.significant <- l.bQTL_gene_partitioning[[s]]
    
    strt<-Sys.time()  
    cl<-makeCluster(min(n.chromosomes, n.cpus))
    registerDoParallel(cl)
    
    l.sets <-  foreach(i = 1:n.chromosomes, .packages=c("seqinr", "VariantAnnotation", "Biostrings")) %dopar% { 
      
      df.nucleotideInPeaks <- data.frame(A = numeric(), G = numeric(), C = numeric(), T = numeric(), motif = character())
      
      df.motif_analysis <- postTotal.significant[-(1:nrow(postTotal.significant)),]
      df.motif_analysis["motif.ref"] <- c()
      df.motif_analysis["motif.mutant"] <- c()
      df.motif_analysis["unique"] <- c()
      df.motif_analysis["isMotifInPeak"] <- c()
      df.motif_analysis["pos.motif"] <- c()
      
      print(paste("processing chromosome ", i))
      
      genome.reference <- DNAString(genome[[i]])
      genome.mutant    <- l.genome.mutant[[i]]
      
      # binding peaks 
      tf_target_bind.sset <- subset(df.peaks, df.peaks$seqnames == i)
      
      # bQTL
      postTotal.significant.i <- subset(postTotal.significant, postTotal.significant$contig == i) 
      
      for(m in 1:length(motifs)){
        
        print(paste("processing motif ", m))
        
        s.motif_offset = v.motif_offset[m]
        
        # returns the starting position of the found motif
        v.reference <- as.numeric(start(matchPattern(motifs[m], genome.reference, fixed = TRUE)))
        v.mutant <- as.numeric(start(matchPattern(motifs[m], genome.mutant, fixed = TRUE)))
        
        vec.intersect <- intersect(v.reference, v.mutant)
        vec.reference.unique <- v.reference[!v.reference %in% v.mutant] # B73
        vec.mutant.unique <- v.mutant[!v.mutant %in% v.reference] # MO17
        
        ## analyze peak based surroundings to 
        for(j in 1:length(vec.intersect)){
          
          # is the motif (not just the ASB) in peak
          idx.j <- which(tf_target_bind.sset$start <= vec.intersect[j] & tf_target_bind.sset$end >= vec.intersect[j] + s.motif_offset) # found a mutual motif in the peak -> take this peak
          
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
    
    stopCluster(cl)
    print(Sys.time()-strt)
    
    df.nucleotideInPeaks <- c()
    df.motif_analysis <- c()
    
    for(i in 1:n.chromosomes){
      df.nucleotideInPeaks <- rbind(df.nucleotideInPeaks, l.sets[[i]][[1]])
      df.motif_analysis <- rbind(df.motif_analysis, l.sets[[i]][[2]])
    }
    
    l.motif_analysis[[s]] <- df.motif_analysis
    l.nucleotideInPeaks[[s]] <- df.nucleotideInPeaks
    
  }
  
  saveRDS(l.motif_analysis, paste(folder_tmp, "l.motif_analysis.rds", sep = "/"))
  saveRDS(l.nucleotideInPeaks, paste(folder_tmp, "l.nucleotideInPeaks.rds", sep = "/"))
  
  message("...finished")
}


#TODO: function - l.motif_analysis, l.nucleotideInPeaks, motifs, v.sets <- c(5,10,15,20, 25)


if(b.predefined_motif_analysis & b.peak_analysis){
  
  # To identify ASBs which may be explained by motif variation (Fig. 2h), we extracted the
  # +/- 5 bp of the high affinity BZR1 bound allele surrounding ASBs. Using R we scanned
  # those 11 bp fragments for canonical BRRE (CGTG[T/C]G, C[G/A]CACG), allowing a
  # single base pair mismatch outside the core motif (CGTG), or G-box (CACGTG) motifs
  # and determined ASBs where the SNP changed a BRRE or G-box motif into an altered
  # (non BRRE or G-box) motif. 
  
  message("Performing predefined motif and peak analysis...")
  
  #df.motif_analysis.significant <- df.motif_analysis
  #df.nucleotideInPeaks.significant <- df.nucleotideInPeaks
  l.motif_analysis.postprocessed <- vector(mode = "list", length = 2)
  
  for(s in 1:2){
    
    df.motif_analysis <- l.motif_analysis[[s]]
    df.nucleotideInPeaks <- l.nucleotideInPeaks[[s]]
    
    # nucleotide environments 
    v.sets <- c(5,10,15,20, 25)
    
    for(k in 1:length(v.sets)){
      for(m in 1:length(motifs)){
        
        df.nucleotideInPeaks.m <- subset(df.nucleotideInPeaks, df.nucleotideInPeaks$motif ==  motifs[m])[,1:4]
        df.nucleotideInPeaks.m[df.nucleotideInPeaks.m > v.sets[k]] <- 0
        df.nucleotideInPeaks.m <- df.nucleotideInPeaks.m[which(rowSums(df.nucleotideInPeaks.m) > 0),]
        
        if(nrow(df.nucleotideInPeaks.m) > 0)
          res <- colSums(df.nucleotideInPeaks.m) / sum(df.nucleotideInPeaks.m)
      }
    }  
    
    df.motif_analysis.postprocessed <- subset(df.motif_analysis, !(df.motif_analysis$motif.ref %in% motifs & df.motif_analysis$motif.mutant %in% motifs))
    
    if(s == 1)
      df.motif_analysis.postprocessed <- subset(df.motif_analysis.postprocessed, df.motif_analysis.postprocessed$altCount >= 5 & df.motif_analysis.postprocessed$refCount >= 5)
    
    if(s == 2)
      df.motif_analysis.postprocessed["POSTfreq"] = 0.5
    
    v.sets <- c("reference", "mutant")
    
    # position_beforeMotifInDataframe <- 29 
    # 
    # if(s == 1){
    #   for(k in 1:2){ # 
    #     
    #     print(v.sets[k]) # bargraph
    #     
    #     if(k == 1){
    #       test <- subset(df.motif_analysis.postprocessed, df.motif_analysis.postprocessed$unique == "reference")
    #     }else{
    #       test <- subset(df.motif_analysis.postprocessed, df.motif_analysis.postprocessed$unique == "mutant")
    #     }
    #     
    #     if(nrow(test) > 0){
    #       for(l in 1:length(unique(test[,position_beforeMotifInDataframe+k]))){
    #         
    #         print(unique(test[,position_beforeMotifInDataframe+k])[l])
    #         
    #         test.l <- subset(test,test[,position_beforeMotifInDataframe+k] == unique(test[,position_beforeMotifInDataframe+k])[l])
    #         a <- test.l$POSTfreq
    #         print(length(a[a < 0.5]))
    #         print(length(a[a > 0.5]))
    #         print("")
    #         if(k == 1){
    #           test.l <- subset(test.l, test.l$POSTfreq > 0.5)
    #         }else{
    #           test.l <- subset(test.l, test.l$POSTfreq < 0.5)
    #         }
    #         snp.pos <- (test.l$position - test.l$pos.motif) + 1
    #         #test.l$POSTallele
    #         
    #         print(table(snp.pos))
    #         
    #         #       for(j in 1:6){
    #         #         print(median(test.l$POSTfreq[which(snp.pos == j)]))
    #         #       }
    #         
    #       }
    #     }
    #   }
    # }
    
    df.motif_analysis.postprocessed["direction"] <- NA
    
    for(k in 1:nrow(df.motif_analysis.postprocessed)){
      
      if(df.motif_analysis.postprocessed$unique[k] == "reference" & df.motif_analysis.postprocessed$POSTfreq[k] > 0.5){
        df.motif_analysis.postprocessed$direction[k] <- "yes"
      }else if(df.motif_analysis.postprocessed$unique[k] == "mutant" & df.motif_analysis.postprocessed$POSTfreq[k] < 0.5){
        df.motif_analysis.postprocessed$direction[k] <- "yes"
      }else{
        df.motif_analysis.postprocessed$direction[k] <- "no"
      }
    }
    
    l.motif_analysis.postprocessed[[s]] <- df.motif_analysis.postprocessed
    
  }
  
  saveRDS(l.motif_analysis.postprocessed, paste(folder_tmp, "l.motif_analysis.postprocessed.rds", sep = "/"))
  
  df.motif_directionality <- data.frame(motif = motifs,  yes = numeric(length(motifs)),  no = numeric(length(motifs)))
  
  for(m in 1:length(motifs)){
    
    test <- l.motif_analysis.postprocessed[[1]] 
    tmp <- subset(test, test$motif.mutant ==  motifs[m] & test$unique == "mutant")
    tmp <- rbind(tmp, subset(test, test$motif.ref ==  motifs[m] & test$unique == "reference"))
    
    df.motif_directionality$yes[m] <- table(tmp$direction)["yes"]
    df.motif_directionality$no[m] <- table(tmp$direction)["no"]
  }
  
  
  df.motif_nucleotidePositionPartitions <- data.frame(motif = character(),
                                                      type = character(),
                                                      inPeak = character(),
                                                      nucleotidePos_1 = numeric(),
                                                      nucleotidePos_2 = numeric(),
                                                      nucleotidePos_3 = numeric(),
                                                      nucleotidePos_4 = numeric(),
                                                      nucleotidePos_5 = numeric(),
                                                      nucleotidePos_6 = numeric())
  
  
  v.inPeak <- which(df.motif_analysis.postprocessed$isMotifInPeak == "no")
  v.direction <- which(df.motif_analysis.postprocessed$direction == "no")
  
  df.motif_analysis.postprocessed <- l.motif_analysis.postprocessed[[1]]
  
  v.inPeak <- which(df.motif_analysis.postprocessed$isMotifInPeak == "yes")
  v.direction <- which(df.motif_analysis.postprocessed$direction == "yes")
  
  length(intersect(v.inPeak,v.direction ))
  length(v.inPeak)
  length(v.direction)
  
  set_sign <- c("ASBs", "bgSNPs")
  set_peak <- c("yes", "no")
  
  peak <- c("yes", "no")
  
  for(s in 1:2){
    
    for(p in 1:2){
      
      df.motif_analysis.postprocessed <- l.motif_analysis.postprocessed[[s]]
      df.motif_analysis.postprocessed <- subset(df.motif_analysis.postprocessed, df.motif_analysis.postprocessed$isMotifInPeak == peak[p])
      
      for(l in 1:length(motifs)){
        
        test.l <- subset(df.motif_analysis.postprocessed, df.motif_analysis.postprocessed$motif.ref == motifs[l] | df.motif_analysis.postprocessed$motif.mutant == motifs[l])
        snp.pos <- (test.l$position - test.l$pos.motif) + 1
        
        # directionality test 
        if(FALSE){
          test.ld <- subset(test.l, test.l$direction == "yes")
          snp.pos <- (test.ld$position - test.ld$pos.motif) + 1
        }
        
        snp.pos = table(snp.pos)
        nucleotidePos = c(0,0,0,0,0,0)
        
        for(i in 1:length(nucleotidePos)){
          if(as.character(i) %in% names(snp.pos)){
            nucleotidePos[i] <- snp.pos[as.character(i)]
          }
        }
        
        
        df.motif_nucleotidePositionPartitions <- rbind(df.motif_nucleotidePositionPartitions, 
                                                       data.frame(motif = motifs[l],
                                                                  type = set_sign[s],
                                                                  isMotifInPeak = set_peak[p],
                                                                  nucleotidePos_1 = nucleotidePos[1],
                                                                  nucleotidePos_2 = nucleotidePos[2],
                                                                  nucleotidePos_3 = nucleotidePos[3],
                                                                  nucleotidePos_4 = nucleotidePos[4],
                                                                  nucleotidePos_5 = nucleotidePos[5],
                                                                  nucleotidePos_6 = nucleotidePos[6]))
        
      }
      
    }
    
  }
  
  message("...finished")
  
  # return 
  # df.motif_directionality
  # df.motif_nucleotidePositionPartitions
  
  
}else{
  message("no predefined motif and in peak analysis")
}



# +/- 5 AND 50 bp (excluding asbs in motifs)
# l.bQTL_gene_partitioning <- readRDS(paste("tmp/l.bQTL_gene_partitioning_withGeneDistances_backgroundSampled_", timeStamp, ".rds", sep = ""))

message("check motif snp position significance")


#### Fig 2 g - motif position analysis 
motif_snp_position_significance <- function(l.bQTL_gene_partitioning,
                                            l.postTotal,
                                            df.peaks,
                                            n.chromosomes,
                                            genome,
                                            l.genome.mutant){
  
  df.motifPositionAnalysis <- c() 
  
  for(i in 1:n.chromosomes){
    
    postTotal.combined <- c()
    
    df.nucleotideInPeaks <- data.frame(A = numeric(), G = numeric(), C = numeric(), T = numeric(), motif = character())
    
    if(TRUE){
      postTotal.significant <- subset(l.bQTL_gene_partitioning[[1]], l.bQTL_gene_partitioning[[1]]$contig == i)  
      #postTotal.chromosome <- subset(df.gene_function_partitioning, df.gene_function_partitioning$contig == i)  
    }else{
      postTotal.significant <- subset(l.postTotal[[1]], l.postTotal[[1]]$contig == i)  
    }
    
    postTotal.significant["ref_motif"] <- NA
    postTotal.significant["alt_motif"] <- NA
    postTotal.significant["distance_snp_to_peak"] <- NA
    
    postTotal.significant["canonical_in_alt"] <- NA
    postTotal.significant["canonical_in_ref"] <- NA
    postTotal.significant["canonical_to_canonical"] <- NA
    
    postTotal.significant["pos_1"] <- NA
    postTotal.significant["pos_2"] <- NA
    postTotal.significant["pos_3"] <- NA
    postTotal.significant["pos_4"] <- NA
    postTotal.significant["pos_5"] <- NA
    postTotal.significant["pos_6"] <- NA
    
    
    genome.reference <- DNAString(genome[[i]])
    genome.mutant    <- l.genome.mutant[[i]]
    
    # binding peaks 
    tf_target_bind.sset <- subset(df.peaks, df.peaks$seqnames == i)
    
    tf_target_bind.sset["start_window"] <- tf_target_bind.sset$posPeak - 20 
    tf_target_bind.sset["end_window"] <-  tf_target_bind.sset$posPeak + 20
    
    # bQTL
    postTotal.chromosome <- subset(postTotal.significant, postTotal.significant$contig == i) 
    
    for(m in 1:length(motifs)){ 
      
      s.motif_offset = v.motif_offset[m]
      
      print(paste("processing motif ", m))
      
      v.reference <- as.numeric(start(matchPattern(motifs[m], genome.reference, fixed = TRUE)))
      v.mutant <- as.numeric(start(matchPattern(motifs[m], genome.mutant, fixed = TRUE)))
      
      # present in both species
      vec.intersect <- intersect(v.reference, v.mutant)
      vec.reference.unique <- v.reference[!v.reference %in% v.mutant] # B73
      vec.mutant.unique <- v.mutant[!v.mutant %in% v.reference] # MO17
      
      postTotal.chromosome.expanded <- c()
      
      # analyze every peak 
      for(l in 1:nrow(postTotal.chromosome)){
        
        # new motif 
        if(any(vec.mutant.unique <= postTotal.chromosome$position[l] & postTotal.chromosome$position[l] <= vec.mutant.unique + s.motif_offset)){
          b.inPeak <- any(postTotal.chromosome$position[l] >= tf_target_bind.sset$start_window & postTotal.chromosome$position[l] <= tf_target_bind.sset$end_window)
          if(b.inPeak){
            idx.Peak <- which(postTotal.chromosome$position[l] >= tf_target_bind.sset$start_window & postTotal.chromosome$position[l] <= tf_target_bind.sset$end_window)
            
            # snp in motif in peak
            idx.pos_motif <- which(vec.mutant.unique <= postTotal.chromosome$position[l] & postTotal.chromosome$position[l] <= vec.mutant.unique + s.motif_offset)
            idx.pos_motif_snp <- abs(vec.mutant.unique[idx.pos_motif] - postTotal.chromosome$position[l]) + 1 # double check 
            
            for(k in 1:length(idx.pos_motif_snp)){
              
              motif.reference <- as.character(subseq(genome.reference, start = vec.mutant.unique[idx.pos_motif[k]], end = vec.mutant.unique[idx.pos_motif[k]] + s.motif_offset))    
              
              if(idx.pos_motif_snp[k] == 1){
                postTotal.chromosome$pos_1[l] = postTotal.chromosome$POSTfreq[l]
              }else if(idx.pos_motif_snp[k] == 2){
                postTotal.chromosome$pos_2[l] = postTotal.chromosome$POSTfreq[l]
              }else if(idx.pos_motif_snp[k] == 3){
                postTotal.chromosome$pos_3[l] = postTotal.chromosome$POSTfreq[l]
              }else if(idx.pos_motif_snp[k] == 4){
                postTotal.chromosome$pos_4[l] = postTotal.chromosome$POSTfreq[l]
              }else if(idx.pos_motif_snp[k] == 5){
                postTotal.chromosome$pos_5[l] = postTotal.chromosome$POSTfreq[l]
              }else if(idx.pos_motif_snp[k] == 6){
                postTotal.chromosome$pos_6[l] = postTotal.chromosome$POSTfreq[l]
              }
              
              postTotal.chromosome$distance_snp_to_peak[l] <- min(abs(postTotal.chromosome$position[l] - tf_target_bind.sset$posPeak[idx.Peak]))
              
              postTotal.chromosome$ref_motif[l] <- motif.reference
              postTotal.chromosome$alt_motif[l] <- motifs[m]
              
              
              if(motif.reference %in% motifs){
                postTotal.chromosome$canonical_to_canonical[l] <- 1
              }else{
                postTotal.chromosome$canonical_in_alt[l] <- 1
                postTotal.chromosome$canonical_in_ref[l] <- 0
              }
              
              postTotal.chromosome.expanded <- rbind(postTotal.chromosome.expanded, postTotal.chromosome[l,])
            }
          }
        }
        
        # destroyed motif 
        if(any(vec.reference.unique <= postTotal.chromosome$position[l] & postTotal.chromosome$position[l] <= vec.reference.unique + s.motif_offset)){
          
          b.inPeak <- any(postTotal.chromosome$position[l] >= tf_target_bind.sset$start_window & postTotal.chromosome$position[l] <= tf_target_bind.sset$end_window)
          
          if(b.inPeak){
            
            idx.Peak <- which(postTotal.chromosome$position[l] >= tf_target_bind.sset$start_window & postTotal.chromosome$position[l] <= tf_target_bind.sset$end_window)
            
            # snp in motif in peak
            idx.pos_motif <- which(vec.reference.unique <= postTotal.chromosome$position[l] & postTotal.chromosome$position[l] <= vec.reference.unique + s.motif_offset)
            idx.pos_motif_snp <- abs(vec.reference.unique[idx.pos_motif] - postTotal.chromosome$position[l]) + 1 # double check 
            
            for(k in 1:length(idx.pos_motif_snp)){
              
              motif.mutant <- as.character(subseq(genome.mutant, start = vec.reference.unique[idx.pos_motif[k]], end = vec.reference.unique[idx.pos_motif[k]] + s.motif_offset))    
              
              if(idx.pos_motif_snp[k] == 1){
                postTotal.chromosome$pos_1[l] = postTotal.chromosome$POSTfreq[l]
              }else if(idx.pos_motif_snp[k] == 2){
                postTotal.chromosome$pos_2[l] = postTotal.chromosome$POSTfreq[l]
              }else if(idx.pos_motif_snp[k] == 3){
                postTotal.chromosome$pos_3[l] = postTotal.chromosome$POSTfreq[l]
              }else if(idx.pos_motif_snp[k] == 4){
                postTotal.chromosome$pos_4[l] = postTotal.chromosome$POSTfreq[l]
              }else if(idx.pos_motif_snp[k] == 5){
                postTotal.chromosome$pos_5[l] = postTotal.chromosome$POSTfreq[l]
              }else if(idx.pos_motif_snp[k] == 6){
                postTotal.chromosome$pos_6[l] = postTotal.chromosome$POSTfreq[l]
              }
              
              postTotal.chromosome$distance_snp_to_peak[l] <- min(abs(postTotal.chromosome$position[l] - tf_target_bind.sset$posPeak[idx.Peak]))
              
              postTotal.chromosome$ref_motif[l] <- motifs[m]
              postTotal.chromosome$alt_motif[l] <- motif.mutant
              
              if(motif.mutant %in% motifs){
                postTotal.chromosome$canonical_to_canonical[l] <- 1
              }else{
                postTotal.chromosome$canonical_in_alt[l] <- 0
                postTotal.chromosome$canonical_in_ref[l] <- 1
              }
              
              postTotal.chromosome.expanded <- rbind(postTotal.chromosome.expanded, postTotal.chromosome[l,])
              
            }
          }
        }
        
      }
      
      postTotal.combined <- rbind(postTotal.combined, postTotal.chromosome.expanded)
      
    }
    
    df.motifPositionAnalysis <- rbind(df.motifPositionAnalysis, postTotal.combined)
    
  }
  
  write.table(df.motifPositionAnalysis, paste(folder_output, "/motifs/S0_5.txt", sep = ""),  row.names = FALSE, quote = FALSE, sep ="\t")
  
  
}

# stopCluster(cl)
# print(Sys.time()-strt)
# 
# 


write.csv(df.motifPositionAnalysis, paste(folder_output, "/df.motifPositionAnalysis_BRRE_GBOX_ABI_TCP.csv", sep =""))



df.motifPositionAnalysis
l.postTotal.chromosomes


# l.motif_analysis



df.motif_analysis <- l.motif_analysis[[1]]

df.motif_analysis.postprocessed <- l.motif_analysis.postprocessed[[1]]

df.test <- subset(df.motif_analysis.postprocessed, df.motif_analysis.postprocessed$unique == "reference")
hist(df.test$POSTfreq, breaks = 100)

hist(subset(df.motif_analysis.postprocessed, df.motif_analysis.postprocessed$unique == "reference")$POSTfreq, breaks = 100) #  & df.motif_analysis$motif.ref %in% motifs.sel)$POSTfreq,breaks = 50)
hist(subset(df.motif_analysis.postprocessed, df.motif_analysis.postprocessed$unique == "mutant")$POSTfreq, breaks = 100) # & df.motif_analysis$motif.mutant %in% motifs.sel)$POSTfreq,breaks = 50)

hist(subset(df.motif_analysis, df.motif_analysis$unique == "reference")$POSTfreq, breaks = 100)#  & df.motif_analysis$motif.ref %in% motifs.sel)$POSTfreq,breaks = 50)
hist(subset(df.motif_analysis, df.motif_analysis$unique == "mutant")$POSTfreq, breaks = 100) # & df.motif_analysis$motif.mutant %in% motifs.sel)$POSTfreq,breaks = 50)

# only to be performed for significant
print("perform de novo motif analysis")



denovo_motif_discovery <- function(postTotal.significant = l.selection[[1]], th.seq = 5){
  
  #postTotal.significant <- l.postTotal[[1]]http://heliumfloats.com/helium_logo.png
  
  # postTotal.significant <- l.selection[[1]]
  
  strt<-Sys.time() 
  cl<-makeCluster(n.cpus)
  registerDoParallel(cl)
  
  l.sequences <- foreach(i = 1:n.chromosomes, .packages=c("seqinr", "VariantAnnotation", "Biostrings")) %dopar% {   
    
    #   df.nucleotideInPeaks <- data.frame(A = numeric(), G = numeric(), C = numeric(), T = numeric(), motif = character())
    #   
    #   df.motif_analysis <- postTotal.significant[-(1:nrow(postTotal.significant)),]
    #   df.motif_analysis["motif.ref"] <- c()
    #   df.motif_analysis["motif.mutant"] <- c()
    #   df.motif_analysis["unique"] <- c()
    #   df.motif_analysis["isInPeak"] <- c()
    #   df.motif_analysis["pos.motif"] <- c()
    
    # for(i in 1:10){
    
    print(paste("processing chromosome ", i))
    
    genome.reference <- DNAString(genome[[i]])
    genome.mutant    <- l.genome.mutant[[i]]
    
    
    # binding peaks 
    tf_target_bind.sset <- subset(df.peaks, df.peaks$seqnames == i)
    
    # bQTL
    postTotal.significant.i <- subset(postTotal.significant, postTotal.significant$contig == i)
    
    
    df.snp_surrounding <- as.data.frame(matrix("", nrow = 0, ncol = 12))
    v.seqs <- character(nrow(postTotal.significant.i))
    
    for(j in 1:nrow(postTotal.significant.i)){
      
      cat("Processing... ", round(j/nrow(postTotal.significant.i) * 100, digits = 2) , "%", "\r"); flush.console()  
      
      i.snp <- postTotal.significant.i$position[j]
      
      if(postTotal.significant.i$POSTfreq[j] < 0.5){
        
        v.seq.mutant <- as.character(subseq(genome.mutant, start=i.snp - th.seq, end=i.snp + th.seq))    
        
        v.seqs[j] <- v.seq.mutant
        
        # for denovo motif
        v.seq.mutant <- (strsplit(v.seq.mutant, ""))
        tmp <- cbind(as.data.frame(t(v.seq.mutant[[1]])), "mo17")
        names(tmp)[12] <- "V42"
        df.snp_surrounding <- rbind(df.snp_surrounding, tmp)
        
      }else{
        
        v.seq.reference <- as.character(subseq(genome.reference, start=i.snp - th.seq, end=i.snp + th.seq))    
        
        v.seqs[j] <- v.seq.reference
        
        # for denovo motif
        v.seq.reference <- (strsplit(v.seq.reference, ""))
        tmp <- cbind(as.data.frame(t(v.seq.reference[[1]])), "B73")
        names(tmp)[12] <- "V42"
        df.snp_surrounding <- rbind(df.snp_surrounding, tmp)
      }
    }
    
    list(v.seqs=v.seqs, df.snp_surrounding=df.snp_surrounding)
  }
  
  stopCluster(cl)
  print(Sys.time()-strt)
  
  
}












#####################



evaluating_BZR1_core_motif <- function(df.ASBs, core_motif = "CGTG"){
  
  core_motif = "CGTG"
  s.motif_offset = 5
  
  df.ASBs.core_motifs = c()
  
  
  for(i in 1:n.chromosomes){
    
    df.ASBs.i = subset(df.ASBs, df.ASBs$contig == i)
    
    df.ASBs.i["B73_sequence"] = ""
    df.ASBs.i["MO17_sequence"] = ""
    df.ASBs.i["B73_pos_core_motif"] = NA
    df.ASBs.i["MO17_pos_core_motif"] = NA
    
    
    print(paste("processing chromosome ", i))
    
    genome.reference <- DNAString(genome[[i]])
    genome.mutant    <- l.genome.mutant[[i]]
    
    for(j in 1:nrow(df.ASBs.i)){
      
      pos.ASB.j = df.ASBs.i$position[j]
      
      motif.reference <- as.character(subseq(genome.reference, start = pos.ASB.j - s.motif_offset, end = pos.ASB.j + s.motif_offset))    
      motif.mutant <- as.character(subseq(genome.mutant, start = pos.ASB.j - s.motif_offset, end = pos.ASB.j + s.motif_offset))    
      
      df.ASBs.i$B73_sequence[j] = motif.reference
      df.ASBs.i$MO17_sequence[j] = motif.mutant
      
      res = unlist(gregexpr(pattern = core_motif, motif.reference))
      if(res != -1){
        res = res - s.motif_offset - 1
        df.ASBs.i$B73_pos_core_motif[j] = paste(res, collapse = ",")
      }
      
      res = unlist(gregexpr(pattern = core_motif, motif.mutant))
      if(res != -1){
        res = res - s.motif_offset - 1
        df.ASBs.i$MO17_pos_core_motif[j] = paste(res, collapse = ",")
      }
      
      
    }
    
    df.ASBs.core_motifs = rbind(df.ASBs.tmp, df.ASBs.i)
    
  }
  
  
  write.table(df.ASBs.core_motifs, paste(folder_output, "/SX3.txt", sep = ""), quote = FALSE, row.names = FALSE, sep ="\t")
  # write.csv2(df.ASBs.tmp, "A:/junkDNA.ai/df.ASBs.csv",row.names = F)
  
  
  ## NAS
  
  
}
























