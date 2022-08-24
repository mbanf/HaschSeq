load_peaks <- function(path.peaks){

  df.peaks <- read.table(path.peaks, header = FALSE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE, quote = "")
  # df.peaks <- df.peaks[2:nrow(df.peaks),] # TODO: do we still need this?
  names(df.peaks) <- c("seqnames", "start", "end", "id", "val")
  #df.peaks$seqnames <- str_replace(df.peaks$seqnames, "chr", "")
  df.peaks["center"] <- as.numeric(unlist(lapply(strsplit(df.peaks$id, ":"), function(m) m[[2]])))
  
  return(df.peaks)
}


load_filtered_peaks <- function(folder){
  return(readRDS(paste(folder, "df.filtered_peaks.rds", sep ="/")))
}


bQTLs_in_peaks <- function(df.bQTLs, df.peaks, n.chromosomes = 10){
  # TODO: DRY
  message("Estimate significant bQTLs in binding peak ranges...")
  
  df.selection <- c()
  
  for(i in 1:n.chromosomes){    
    
    chr <- paste("B73-chr", i, sep = "")
    df.peaks.B73.i <- subset(df.peaks, df.peaks$seqnames == chr)
    df.bQTLs.B73.i <- subset(df.bQTLs, df.bQTLs$`B73-chr` == chr)
    df.bQTLs.B73.i <- subset(df.bQTLs.B73.i, df.bQTLs.B73.i$POSTfreq > 0.5)
    
    message(chr)
    for(p in 1:nrow(df.peaks.B73.i)){
      idx.snps <- which(df.peaks.B73.i$start[p] <= df.bQTLs.B73.i$`B73-pos` & df.bQTLs.B73.i$`B73-pos` <= df.peaks.B73.i$end[p])
      if(length(idx.snps) > 0){
        df.selection <- rbind(df.selection, df.bQTLs.B73.i[idx.snps,])
      }
    }
    
    chr <- paste("Mo17-chr", i, sep = "")
    df.peaks.Mo17.i <- subset(df.peaks, df.peaks$seqnames == chr)
    df.bQTLs.Mo17.i <- subset(df.bQTLs, df.bQTLs$`Mo17-chr` == chr)
    df.bQTLs.Mo17.i <- subset(df.bQTLs.Mo17.i, df.bQTLs.Mo17.i$POSTfreq < 0.5)
    
    message(chr)
    for(p in 1:nrow(df.peaks.Mo17.i)){
      idx.snps <- which(df.peaks.Mo17.i$start[p] <= df.bQTLs.Mo17.i$`Mo17-pos` & df.bQTLs.Mo17.i$`Mo17-pos` <= df.peaks.Mo17.i$end[p])
      if(length(idx.snps) > 0){
        df.selection <- rbind(df.selection, df.bQTLs.Mo17.i[idx.snps,])
      }
    }
   
  }

  return(unique(df.selection))
}




merge_peaks <- function(path.peaks,
                         v.paths.peaks.replicates,
                         i.ref=1,
                         th.min_number = 5,
                         folder = "tmp/"
){
  
  l.df.peaks <- vector(mode = "list", length = length(v.paths.peaks.replicates))
  for(i in 1:length(v.paths.peaks.replicates)){
    l.df.peaks[[i]] <- load_peaks(v.paths.peaks.replicates[[i]])
  }
  
  print(table(l.df.peaks[[i.ref]]$seqnames))
  v.chromosomes <- unique(l.df.peaks[[i.ref]]$seqnames)
  
  message("Merge peak replicate files, using file ", v.paths.peaks.replicates[i.ref], " as reference, with min number of replicates to find peak found in ref. file being ", th.min_number)
  
  strt<-Sys.time() 
  
  df.merged_peaks <- c()
  
  for(i in 1:length(v.chromosomes)){
    
    message(v.chromosomes[i])
    
    l.df.peaks.i <- vector(mode = "list", length = length(l.df.peaks))
    l.idx.blacklist <- vector(mode = "list", length = length(l.df.peaks))
    for(l in 1:length(l.df.peaks)){
      l.df.peaks.i[[l]] <- subset(l.df.peaks[[l]], l.df.peaks[[l]]$seqnames == v.chromosomes[i])
      l.idx.blacklist[[l]] <- numeric(nrow(l.df.peaks.i[[l]]))
    }
    
    v.idx.keep <- c()
    pb <- txtProgressBar(min = 0, max = nrow(l.df.peaks.i[[i.ref]]), style = 3)
    for(j in 1:nrow(l.df.peaks.i[[i.ref]])){
      setTxtProgressBar(pb, j)
      n.found <- 0
      for(l in 2:length(l.df.peaks.i)){
        b.found = FALSE
        for(k in 1:nrow(l.df.peaks.i[[l]])){
          if(!k %in% l.idx.blacklist[[l]]){        
            if(l.df.peaks.i[[i.ref]]$start[j] < l.df.peaks.i[[l]]$center[k] & l.df.peaks.i[[l]]$center[k] < l.df.peaks.i[[i.ref]]$end[j]){
              n.found = n.found + 1
              b.found = TRUE
              l.idx.blacklist[[l]][k] <- k
              break
            }
          }
        }
        if(b.found == FALSE){
          break
        }
      }
      if(n.found == th.min_number){
        v.idx.keep <- c(v.idx.keep, j)  
      }
    }
    close(pb) 
    
    df.merged_peaks <- rbind(df.merged_peaks, l.df.peaks.i[[i.ref]][v.idx.keep,])
  }
  
  print(Sys.time()-strt)
  
  saveRDS(df.merged_peaks, paste(folder, "df.merged_peaks.rds", sep ="/"))  
  write.csv(df.merged_peaks, paste(folder, "df.merged_peaks.csv", sep ="/"), row.names = F, quote = F)  
  
}



filter_peaks_w_merged_replicaties <- function(path.peaks,
                                       th.dist = 500,
                                       folder = "tmp/"
                                        ){
  
  
  df.peaks <- load_peaks(path.peaks)
  
  print(table(df.peaks$seqnames))
  v.chromosomes <- unique(df.peaks$seqnames)
  
  df.merged_peaks <- readRDS(paste(folder, "df.merged_peaks.rds", sep ="/"))  
    
  message("filter original peaks by merged input read replicates - exclude all with center distance > ", th.dist)
  
  df.filtered_peaks <- c()
  strt<-Sys.time() 
  
  for(i in 1:length(v.chromosomes)){
  
    message(v.chromosomes[i])
    
    df.peaks.i <- subset(df.peaks, df.peaks$seqnames == v.chromosomes[i])
    df.merged_peaks.i <- subset(df.merged_peaks, df.merged_peaks$seqnames == v.chromosomes[i])
    
    v.idx.keep <- c()
    pb <- txtProgressBar(min = 0, max = nrow(df.peaks.i), style = 3)
    for(j in 1:nrow(df.peaks.i)){
      setTxtProgressBar(pb, j)
      for(k in 1:nrow(df.merged_peaks.i)){
        dist <- abs(df.peaks.i$center[j]  -  df.merged_peaks.i$center)  
        if(dist < th.dist){
          v.idx.keep <- c(v.idx.keep, j)
          break
        }
      }
    }
    close(pb)
    
    df.filtered_peaks <- rbind(df.filtered_peaks, df.peaks.i[v.idx.keep, ])
  }
  
  print(Sys.time()-strt)
  
  write.csv(df.filtered_peaks, paste(folder, "df.filtered_peaks.csv", sep ="/"), row.names = F, quote = F)  
  saveRDS(df.filtered_peaks, paste(folder, "df.filtered_peaks.rds", sep ="/"))  

}
