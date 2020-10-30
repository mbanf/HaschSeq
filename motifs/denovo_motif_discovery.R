
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



