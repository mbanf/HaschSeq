asb_distances_to_average_gene = function(){
  
  v.paths.bindingPeaks <- c("datasets_paper/GEM/AGTCAA/q2_3fold/q2_3fold.GEM_events.txt",
                            "datasets_paper/GEM/ATGTCA/q2_3fold/q2_3fold.GEM_events.txt",
                            "datasets_paper/GEM/CCGTCC/q2_3fold/q2_3fold.GEM_events.txt")
  
  # needs to be replaced by df.geneAnnotation - replace annotations
  df.gene_only_annotation <- subset(df.gene_annotation, df.gene_annotation$partition == "gene")
  
  
  message("peak distances to TSS, gene, stop")
  
  ## novel B73 datasets - updated 06.2018 
  
  
  ## old - check if needed 
  dist_to_TSS <- 5000
  dist_to_stop <- 5000
  
  l.gene_distances <- vector(mode= "list", length = length(v.paths.bindingPeaks))
  
  for(s in 1:length(v.paths.bindingPeaks)){
    
    print(s)
    
    df.peaks <- read.table(v.paths.bindingPeaks[s], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    df.peaks["seqnames"] <- as.numeric(unlist(lapply(strsplit(df.peaks$Position, ":"), function(m) m[[1]])))
    df.peaks["posPeak"] <- as.numeric(unlist(lapply(strsplit(df.peaks$Position, ":"), function(m) m[[2]])))
    df.peaks["start"] <- df.peaks$posPeak - 75
    df.peaks["end"] <- df.peaks$posPeak + 75
    
    df.snp_gene_distances <- c()
    
    for(i in 1:n.chromosomes){
      
      print(i)
      
      df.peaks.i <- subset(df.peaks, df.peaks$seqnames == i)  
      df.gene_only_annotation.i <- subset(df.gene_only_annotation, df.gene_only_annotation$chr == i)
      
      
      df.peaks.i.pos <- subset(df.peaks.i, df.peaks.i$Strand %in% c("+","*"))
      
      df.peaks.i.pos["distance_to_TSS"] <- 1e20
      df.peaks.i.pos["in_gene"] <- "no"
      df.peaks.i.pos["distance_to_stop"] <- 1e20
      
      for(k in 1:nrow(df.peaks.i.pos)){
        
        cat("Processing... ", round(k/nrow(df.peaks.i.pos) * 100, digits = 3) , "%", "\r"); flush.console()  
        
        ### TEST 1 - distance to TSS
        
        df.gene_only_annotation.i.pos <- subset(df.gene_only_annotation.i, df.gene_only_annotation.i$strand == "+")
        start <- df.gene_only_annotation.i.pos$pos.start - dist_to_TSS
        end <- df.gene_only_annotation.i.pos$pos.stop # TSS
        
        i.set <- which(start < df.peaks.i.pos[,"posPeak"][k])
        
        start <- start[i.set]
        end <- end[i.set]
        
        i.set <- which(df.peaks.i.pos[,"posPeak"][k] < end)
        
        dist.pos <- 1e20
        if(length(i.set) > 0){
          start <- start[i.set]
          end <- end[i.set]
          dist.pos <- min(abs(end - df.peaks.i.pos[,"posPeak"][k]))
        }
        
        df.gene_only_annotation.i.neg <- subset(df.gene_only_annotation.i, df.gene_only_annotation.i$strand == "-")
        start <- df.gene_only_annotation.i.neg$pos.stop # TSS
        end <- df.gene_only_annotation.i.neg$pos.stop + dist_to_TSS
        
        i.set <- which(start < df.peaks.i.pos[,"posPeak"][k])
        
        start <- start[i.set]
        end <- end[i.set]
        
        i.set <- which(df.peaks.i.pos[,"posPeak"][k] < end)
        
        if(length(i.set) > 0){
          start <- start[i.set]
          end <- end[i.set]
          df.peaks.i.pos$distance_to_TSS[k] <- min(dist.pos, min(abs(start - df.peaks.i.pos[,"posPeak"][k])))
        }else{
          df.peaks.i.pos$distance_to_TSS[k] <- dist.pos
        }
        
        
        #### TEST 2 - in gene 
        
        df.gff.i.pos <- subset(df.gff.i, df.gff.i$V7 == "+")
        start <- df.gff.i.pos$V4
        end <- df.gff.i.pos$V5
        
        i.set <- which(start < df.test[,"posPeak"][k])
        
        start <- start[i.set]
        end <- end[i.set]
        
        i.set <- which(df.test[,"posPeak"][k] < end)
        
        if(length(i.set) > 0){
          df.test$in_gene[k] <- "yes"
        }
        
        df.gff.i.neg <- subset(df.gff.i, df.gff.i$V7 == "-")
        start <- df.gff.i.neg$V4
        end <- df.gff.i.neg$V5
        
        
        i.set <- which(start < df.test[,"posPeak"][k])
        
        start <- start[i.set]
        end <- end[i.set]
        
        i.set <- which(df.test[,"posPeak"][k] < end)
        
        if(length(i.set) > 0){
          df.test$in_gene[k] <- "yes"
        }
        
        
        # Test 3 - behind gene distance 
        
        df.gff.i.pos <- subset(df.gff.i, df.gff.i$V7 == "+")
        start <- df.gff.i.pos$V5 # stop
        end <- df.gff.i.pos$V5 + dist_to_stop
        
        i.set <- which(start < df.test[,"posPeak"][k])
        
        start <- start[i.set]
        end <- end[i.set]
        
        i.set <- which(df.test[,"posPeak"][k] < end)
        
        
        dist.pos <- 1e20
        if(length(i.set) > 0){
          start <- start[i.set]
          end <- end[i.set]
          dist.pos <- min(abs(start - df.test[,"posPeak"][k]))
        }
        
        ###
        
        df.gff.i.neg <- subset(df.gff.i, df.gff.i$V7 == "-")
        start <- df.gff.i.neg$V4 - dist_to_stop
        end <- df.gff.i.neg$V4
        
        
        i.set <- which(start < df.test[,"posPeak"][k])
        
        start <- start[i.set]
        end <- end[i.set]
        
        i.set <- which(df.test[,"posPeak"][k] < end)
        
        if(length(i.set) > 0){
          start <- start[i.set]
          end <- end[i.set]
          df.test$distance_to_stop[k] <- min(dist.pos, min(abs(end - df.test[,"posPeak"][k])))
        }else{
          df.test$distance_to_stop[k] <- dist.pos
        }
        
      }
      
      ###
      
      for(k in 1:nrow(df.test)){
        
        cat("Processing... ", round(k/nrow(df.test) * 100, digits = 3) , "%", "\r"); flush.console()  
        
        ### TEST 1 - distance to TSS
        
        df.gene_only_annotation.i.pos <- subset(df.gene_only_annotation.i, df.gene_only_annotation.i$strand == "+")
        start <- df.gene_only_annotation.i.pos$pos.start - dist_to_TSS
        end <- df.gene_only_annotation.i.pos$pos.stop # TSS
        
        i.set <- which(start < df.test[,"posPeak"][k])
        
        start <- start[i.set]
        end <- end[i.set]
        
        i.set <- which(df.test[,"posPeak"][k] < end)
        
        dist.pos <- 1e20
        if(length(i.set) > 0){
          start <- start[i.set]
          end <- end[i.set]
          dist.pos <- min(abs(end - df.test[,"posPeak"][k]))
        }
        
        df.gene_only_annotation.i.neg <- subset(df.gene_only_annotation.i, df.gene_only_annotation.i$strand == "-")
        start <- df.gene_only_annotation.i.neg$pos.stop # TSS
        end <- df.gene_only_annotation.i.neg$pos.stop + dist_to_TSS
        
        i.set <- which(start < df.test[,"posPeak"][k])
        
        start <- start[i.set]
        end <- end[i.set]
        
        i.set <- which(df.test[,"posPeak"][k] < end)
        
        if(length(i.set) > 0){
          start <- start[i.set]
          end <- end[i.set]
          df.test$distance_to_TSS[k] <- min(dist.pos, min(abs(start - df.test[,"posPeak"][k])))
        }else{
          df.test$distance_to_TSS[k] <- dist.pos
        }
        
        
        #### TEST 2 - in gene 
        
        df.gff.i.pos <- subset(df.gff.i, df.gff.i$V7 == "+")
        start <- df.gff.i.pos$V4
        end <- df.gff.i.pos$V5
        
        i.set <- which(start < df.test[,"posPeak"][k])
        
        start <- start[i.set]
        end <- end[i.set]
        
        i.set <- which(df.test[,"posPeak"][k] < end)
        
        if(length(i.set) > 0){
          df.test$in_gene[k] <- "yes"
        }
        
        df.gff.i.neg <- subset(df.gff.i, df.gff.i$V7 == "-")
        start <- df.gff.i.neg$V4
        end <- df.gff.i.neg$V5
        
        
        i.set <- which(start < df.test[,"posPeak"][k])
        
        start <- start[i.set]
        end <- end[i.set]
        
        i.set <- which(df.test[,"posPeak"][k] < end)
        
        if(length(i.set) > 0){
          df.test$in_gene[k] <- "yes"
        }
        
        
        # Test 3 - behind gene distance 
        
        df.gff.i.pos <- subset(df.gff.i, df.gff.i$V7 == "+")
        start <- df.gff.i.pos$V5 # stop
        end <- df.gff.i.pos$V5 + dist_to_stop
        
        i.set <- which(start < df.test[,"posPeak"][k])
        
        start <- start[i.set]
        end <- end[i.set]
        
        i.set <- which(df.test[,"posPeak"][k] < end)
        
        
        dist.pos <- 1e20
        if(length(i.set) > 0){
          start <- start[i.set]
          end <- end[i.set]
          dist.pos <- min(abs(start - df.test[,"posPeak"][k]))
        }
        
        ###
        
        df.gff.i.neg <- subset(df.gff.i, df.gff.i$V7 == "-")
        start <- df.gff.i.neg$V4 - dist_to_stop
        end <- df.gff.i.neg$V4
        
        
        i.set <- which(start < df.test[,"posPeak"][k])
        
        start <- start[i.set]
        end <- end[i.set]
        
        i.set <- which(df.test[,"posPeak"][k] < end)
        
        if(length(i.set) > 0){
          start <- start[i.set]
          end <- end[i.set]
          df.test$distance_to_stop[k] <- min(dist.pos, min(abs(end - df.test[,"posPeak"][k])))
        }else{
          df.test$distance_to_stop[k] <- dist.pos
        }
        
      }
      
      
      
      ###
      
      df.snp_gene_distances <- rbind(df.snp_gene_distances, df.test)
      
    }
    
    l.gene_distances[[s]] <- df.snp_gene_distances
    
  }
  
  saveRDS(l.gene_distances, paste("tmp/l.gene_distances_", timeStamp, ".rds", sep = ""))    
}