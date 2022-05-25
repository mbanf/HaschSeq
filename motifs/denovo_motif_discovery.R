
denovo_motif_discovery <- function(postTotal.significant = l.selection[[1]], df.peaks, genome, l.genome.mutant, th.seq = 5){

  # For local motif enrichment analysis (Fig. 2E), we extracted +/- 50 bp of the high affinity
  # BZR1 bound allele surrounding ASBs. The MEME CentriMo suite (v. 5.2.0) was used to
  # determine the local distribution along the 101 bp fragments for the canonical BZR1 motifs
  # BRRE (CGTG[T/C]G, C[G/A]CACG) or G-box (CACGTG) and a control motif, with SBP
  # ("GTACGG", "CCGTAC")14, with a similar GC content
  
  # for external extraction 
  
  strt<-Sys.time() 
  cl<-makeCluster(n.cpus)
  registerDoParallel(cl)
  
  l.sequences <- foreach(i = 1:n.chromosomes, .packages=c("seqinr", "VariantAnnotation", "Biostrings")) %dopar% {   
    
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
  
  # store for external analysis
  
}



