
linkage_disequlilibrium <- function(df.bQTLs,
                                    df.peaks,
                                    mode = "peak_distance",
                                    s.disequilibriumDistance = 75,
                                    n.chromosomes = 10){

  df.bQTLs.peak_distance <- df.bQTLs[-(1:nrow(df.bQTLs)),]
  
  for(i in 1:n.chromosomes){    
    
    chr <- paste("B73-chr", i, sep = "")
    df.peaks.B73.i <- subset(df.peaks, df.peaks$seqnames == chr)
    
    df.bQTLs.B73.i <- subset(df.bQTLs, df.bQTLs$`B73-chr` == chr)
    df.bQTLs.B73.i <- subset(df.bQTLs.B73.i, df.bQTLs.B73.i$POSTfreq > 0.5)
    df.bQTLs.B73.i <- df.bQTLs.B73.i[order(df.bQTLs.B73.i$`B73-pos`), ]
    
    while(TRUE){
      dist <- df.bQTLs.B73.i$`B73-pos`[2:nrow(df.bQTLs.B73.i)] - df.bQTLs.B73.i$`B73-pos`[1:(nrow(df.bQTLs.B73.i) - 1)]
      idx.dist <- which(dist < s.disequilibriumDistance)
      if(length(idx.dist) == 0){
        break
      }
      idx.df <- idx.dist[1]
      set <- df.bQTLs.B73.i[c(idx.df,(idx.df+1)),]
      if(mode == "peak_distance"){
        i.max <- ifelse(min(abs(set$`B73-pos`[1] - df.peaks.B73.i$center)) < min(abs(set$`B73-pos`[2] - df.peaks.B73.i$center)), idx.df + 1, idx.df)
      }else{
        i.max <- ifelse(set$`p-value (corrected)`[1] < set$`p-value (corrected)`[2], idx.df + 1, idx.df)
      }
      df.bQTLs.B73.i <- df.bQTLs.B73.i[-i.max,]
    }
    
    df.bQTLs.peak_distance <- rbind(df.bQTLs.peak_distance, df.bQTLs.B73.i)
    
    #### Mo17 ###
    
    chr <- paste("Mo17-chr", i, sep = "")
    df.bQTLs.Mo17.i <- subset(df.bQTLs, df.bQTLs$`Mo17-chr` == chr)
    df.bQTLs.Mo17.i <- subset(df.bQTLs.Mo17.i, df.bQTLs.Mo17.i$POSTfreq < 0.5)
    df.bQTLs.Mo17.i <- df.bQTLs.Mo17.i[order(df.bQTLs.Mo17.i$`Mo17-pos`), ]
    
    df.peaks.Mo17.i <- subset(df.peaks, df.peaks$seqnames == chr)
    
    while(TRUE){
      
      dist <- df.bQTLs.Mo17.i$`Mo17-pos`[2:nrow(df.bQTLs.Mo17.i)] - df.bQTLs.Mo17.i$`Mo17-pos`[1:(nrow(df.bQTLs.Mo17.i) - 1)]
      idx.dist <- which(dist < s.disequilibriumDistance)
      
      if(length(idx.dist) == 0){
        break
      }
      
      idx.df <-idx.dist[1]
      set <- df.bQTLs.Mo17.i[c(idx.df,(idx.df+1)),]
      if(mode == "peak_distance"){
        i.max <- ifelse(min(abs(set$`Mo17-pos`[1] - df.peaks.Mo17.i$center)) < min(abs(set$`Mo17-pos`[2] - df.peaks.Mo17.i$center)), idx.df + 1, idx.df)
      }else{
        i.max <- ifelse(set$`p-value (corrected)`[1] < set$`p-value (corrected)`[2], idx.df + 1, idx.df)
      }
      df.bQTLs.Mo17.i <- df.bQTLs.Mo17.i[-i.max,]
    }
    df.bQTLs.peak_distance <- rbind(df.bQTLs.peak_distance, df.bQTLs.Mo17.i)
  }
  
  ## sanity check 
  for(i in 1:n.chromosomes){    
    
    chr <- paste("B73-chr", i, sep = "")
    df.bQTLs.B73.i <- subset(df.bQTLs.peak_distance, df.bQTLs.peak_distance$`B73-chr` == chr)
    df.bQTLs.B73.i <- subset(df.bQTLs.B73.i, df.bQTLs.B73.i$POSTfreq > 0.5)
    df.bQTLs.B73.i <- df.bQTLs.B73.i[order(df.bQTLs.B73.i$`B73-pos`), ]
    
    dist <- df.bQTLs.B73.i$`B73-pos`[2:nrow(df.bQTLs.B73.i)] - df.bQTLs.B73.i$`B73-pos`[1:(nrow(df.bQTLs.B73.i) - 1)]
    if(length(which(dist < s.disequilibriumDistance)) > 0){
      print("warning - still linkage")
    }
    
    chr <- paste("Mo17-chr", i, sep = "")
    df.bQTLs.Mo17.i <- subset(df.bQTLs.peak_distance, df.bQTLs.peak_distance$`Mo17-chr` == chr)
    df.bQTLs.Mo17.i <- subset(df.bQTLs.Mo17.i, df.bQTLs.Mo17.i$POSTfreq < 0.5)
    df.bQTLs.Mo17.i <- df.bQTLs.Mo17.i[order(df.bQTLs.Mo17.i$`Mo17-pos`), ]
    
    dist <- df.bQTLs.Mo17.i$`Mo17-pos`[2:nrow(df.bQTLs.Mo17.i)] - df.bQTLs.Mo17.i$`Mo17-pos`[1:(nrow(df.bQTLs.Mo17.i) - 1)]
    if(length(which(dist < s.disequilibriumDistance)) > 0){
      print("warning - still linkage")
    }
    
  }
  
  df.bQTLs.peak_distance
}
