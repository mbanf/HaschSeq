# Quality control of the resulting peaks - needs to be present in 5 out of 6 replicates

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
  }
  print(Sys.time()-strt)
  
}else{
  
  l.peaks.filtered <- readRDS(paste(folder_tmp,"l.peaks.filtered.rds", sep="/"))  
  
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
    saveRDS(df.peaks.final, paste(folder_tmp, "df.peaks.final.rds", sep ="/"))  
  }
  
}