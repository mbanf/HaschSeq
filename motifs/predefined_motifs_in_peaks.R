predefined_motifs_in_peaks = function(){
  
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
        #print( motifs[m])  
        
        df.nucleotideInPeaks.m <- subset(df.nucleotideInPeaks, df.nucleotideInPeaks$motif ==  motifs[m])[,1:4]
        df.nucleotideInPeaks.m[df.nucleotideInPeaks.m > v.sets[k]] <- 0
        
        df.nucleotideInPeaks.m <- df.nucleotideInPeaks.m[which(rowSums(df.nucleotideInPeaks.m) > 0),]
        
        #print(colSums(df.nucleotideInPeaks.m) / sum(df.nucleotideInPeaks.m))
        if(nrow(df.nucleotideInPeaks.m) > 0)
          res <- colSums(df.nucleotideInPeaks.m) / sum(df.nucleotideInPeaks.m)
        #print(res[2] + res[3])
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
}