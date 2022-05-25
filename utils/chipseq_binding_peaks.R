chipseq_binding_peaks <- function(l.path.ChipSeqB73, s.half_window_size.input = 201, n.chromosomes = 10){

  l.df.ChipSeqB73 <- vector(mode = "list", length = length(l.path.ChipSeqB73))
  for(i in 1:length(l.path.ChipSeqB73)){
    l.df.ChipSeqB73[[i]] <- read.table(l.path.ChipSeqB73[[i]], header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE, quote = "")  
    l.df.ChipSeqB73[[i]]["chr"] <- as.numeric(gsub("\\:.*", "", l.df.ChipSeqB73[[i]]$Position))
    l.df.ChipSeqB73[[i]]["pos"] <- as.numeric(gsub(".*:", "", l.df.ChipSeqB73[[i]]$Position))
    l.df.ChipSeqB73[[i]]["left"] <- l.df.ChipSeqB73[[i]]$pos - s.half_window_size.input
    l.df.ChipSeqB73[[i]]["right"] <- l.df.ChipSeqB73[[i]]$pos + s.half_window_size.input
  }
  
  strt<-Sys.time() 
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
  }
  print(Sys.time()-strt)
  
  saveRDS(df.ChipSeq.filtered, paste(folder_tmp,"df.ChipSeq.filtered.rds", sep ="/"))
}