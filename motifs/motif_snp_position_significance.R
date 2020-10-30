motif_snp_position_significance <- function(){
  
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