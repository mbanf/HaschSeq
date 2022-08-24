extract_sequences_affinity_based <- function(df.bQTLs,
                                            genome_sequences,
                                            s.env_bps = 100,
                                            n.chromosomes = 10
){
  
  df.bQTLs["id"] <- seq(1, nrow(df.bQTLs))
  si <- seqinfo(genome_sequences)
  message("extract sequences surrounding ASBs from diploid genome")
  
  v.sequences <- character(length = nrow(df.bQTLs))
  
  for(i in 1:n.chromosomes){ 
    print(i)
    
    chr_name <- paste("B73-chr", i, sep = "")
    chr_seq <- getSeq(genome_sequences, GRanges(chr_name, IRanges(1, seqlengths(si)[chr_name])))
    
    df.bQTLs.B73.i <- subset(df.bQTLs, df.bQTLs$`B73-chr` == chr_name)
    df.bQTLs.B73.i <- subset(df.bQTLs.B73.i, df.bQTLs.B73.i$POSTfreq > 0.5)
    
    for(j in 1:nrow(df.bQTLs.B73.i)){
      id <- df.bQTLs.B73.i$id[j]
      v.sequences[id] <- subseq(chr_seq, start=df.bQTLs.B73.i$`B73-pos`[j] - s.env_bps, end=df.bQTLs.B73.i$`B73-pos`[j] + s.env_bps)  
    }
    
    chr_name <- paste("Mo17-chr", i, sep = "")
    chr_seq <- getSeq(genome_sequences, GRanges(chr_name, IRanges(1, seqlengths(si)[chr_name])))
    
    df.bQTLs.Mo17.i <- subset(df.bQTLs, df.bQTLs$`Mo17-chr` == chr_name)
    df.bQTLs.Mo17.i <- subset(df.bQTLs.Mo17.i, df.bQTLs.Mo17.i$POSTfreq < 0.5)
    
    for(j in 1:nrow(df.bQTLs.Mo17.i)){
      id <- df.bQTLs.Mo17.i$id[j]
      v.sequences[id] <- subseq(chr_seq, start=df.bQTLs.Mo17.i$`Mo17-pos`[j] - s.env_bps, end=df.bQTLs.Mo17.i$`Mo17-pos`[j] + s.env_bps)  
    }
  }
  
  descs <-  paste("asb", df.bQTLs$id, sep = "")
  asb_environment <- DNAStringSet(v.sequences)
  names(asb_environment) <- descs
  
  asb_environment
}


extract_sequences <- function(df.bQTLs,
                               genome_sequences,
                               s.env_bps = 5,
                               n.chromosomes = 10
){
  
  df.bQTLs["id"] <- seq(1, nrow(df.bQTLs))
  descs <-  paste("asb", df.bQTLs$id, sep = "")
  
  si <- seqinfo(genome_sequences)
  message("extract sequences surrounding ASBs from diploid genome")
  
  v.sequences <- character(length = nrow(df.bQTLs))
  
  for(i in 1:n.chromosomes){ 

    chr_name <- paste("B73-chr", i, sep = "")
    chr_seq <- getSeq(genome_sequences, GRanges(chr_name, IRanges(1, seqlengths(si)[chr_name])))
    df.bQTLs.B73.i <- subset(df.bQTLs, df.bQTLs$`B73-chr` == chr_name)

    for(j in 1:nrow(df.bQTLs.B73.i)){
      id <- df.bQTLs.B73.i$id[j]
      v.sequences[id] <- subseq(chr_seq, start=df.bQTLs.B73.i$`B73-pos`[j] - s.env_bps, end=df.bQTLs.B73.i$`B73-pos`[j] + s.env_bps)  
    }
    
  }
  
  
  asb_environment.B73 <- DNAStringSet(v.sequences)
  names(asb_environment.B73) <- descs
  
  v.sequences <- character(length = nrow(df.bQTLs))
  for(i in 1:n.chromosomes){ 
    
    chr_name <- paste("Mo17-chr", i, sep = "")
    chr_seq <- getSeq(genome_sequences, GRanges(chr_name, IRanges(1, seqlengths(si)[chr_name])))
    df.bQTLs.Mo17.i <- subset(df.bQTLs, df.bQTLs$`Mo17-chr` == chr_name)

    for(j in 1:nrow(df.bQTLs.Mo17.i)){
      id <- df.bQTLs.Mo17.i$id[j]
      v.sequences[id] <- subseq(chr_seq, start=df.bQTLs.Mo17.i$`Mo17-pos`[j] - s.env_bps, end=df.bQTLs.Mo17.i$`Mo17-pos`[j] + s.env_bps)  
    }
  }
  
  asb_environment.Mo17 <- DNAStringSet(v.sequences)
  names(asb_environment.Mo17) <- descs
  
  return(list(asb_environment.B73=asb_environment.B73, asb_environment.Mo17=asb_environment.Mo17))
  
}


predefined_motif_analysis <- function(df.bQTLs, 
                                      v.motifs, 

                                      s.env_bps = 5,
                                      n.chromosomes = 10,
                     
                                      path.genomic_sequences = "data/revision/ref_B73Mo17.fasta"
                                      ){
  
  # To identify ASBs which may be explained by motif variation (Fig. 2h), we extracted the
  # +/- 5 bp of the high affinity BZR1 bound allele surrounding ASBs. Using R we scanned
  # those 11 bp fragments for canonical BRRE (CGTG[T/C]G, C[G/A]CACG), allowing a
  # single base pair mismatch outside the core motif (CGTG), or G-box (CACGTG) motifs
  # and determined ASBs where the SNP changed a BRRE or G-box motif into an altered
  # (non BRRE or G-box) motif. 
  
  message("Extracting +/- 5 bp of the high affinity BZR1 bound allele surrounding ASBs and scanning for canonical BRRE (CGTG[T/C]G, C[G/A]CACG), \n allowing a single base pair mismatch outside the core motif (CGTG), or G-box (CACGTG) motifs")
  genome_sequences = readDNAStringSet(path.genomic_sequences)

  res.asb_environment <- extract_sequences(df.bQTLs,
                                            genome_sequences,
                                            s.env_bps,
                                            n.chromosomes)

  rm(genome_sequences)
  
  
  # Fraction of ASBs in motifs (and nucleotide position)
  m.motif_variation <- matrix(0, nrow = length(v.motifs), ncol = s.env_bps + 1) 
  rownames(m.motif_variation) <- v.motifs
  
  m.motif_asbs <- matrix(0, nrow = length(v.motifs), ncol = 6) 
  rownames(m.motif_asbs) <- v.motifs
  colnames(m.motif_asbs) <- c("B73", "Mo17", "dir. incorrect", "dir. correct", "ambiguous (both)", "total (non-ambiguous)")
  
  df.bQTLs[,v.motifs] <- F
  s.asb.pos <- s.env_bps + 1
  
  for(motif in v.motifs){
  
    print(motif)
    
    v.matches.B73 <- vmatchPattern(motif, res.asb_environment$asb_environment.B73, fixed = TRUE)
    v.matches.Mo17 <- vmatchPattern(motif, res.asb_environment$asb_environment.Mo17, fixed = TRUE)
    
    for(i in 1:nrow(df.bQTLs)){
    
      n.B73 <- length(which(as.numeric(lapply(v.matches.B73[i], length)) > 0))
      n.Mo17 <- length(which(as.numeric(lapply(v.matches.Mo17[i], length)) > 0))
      
      if(n.B73 == 0 & n.Mo17 == 0){
        next
      }
      
      if(n.B73 >= 1 | n.Mo17 >= 1)
        df.bQTLs[i, motif] <- T
      
      if(n.B73 >= 1 & n.Mo17 >= 1){
        m.motif_asbs[motif, 5] <- m.motif_asbs[motif, 5] + 1
        next
      }
      
      m.motif_asbs[motif, 6] <- m.motif_asbs[motif, 6] + 1

      if(df.bQTLs$POSTfreq[i] > 0.5){
        if(n.B73 >= 1 & n.Mo17 == 0){
          # only for correct direction 
          v.pos <- s.asb.pos - start(v.matches.B73[[i]]) + 1
          m.motif_variation[motif,v.pos] <- m.motif_variation[motif,v.pos] + 1
          
          m.motif_asbs[motif, 4] <- m.motif_asbs[motif, 4] + 1
          m.motif_asbs[motif, 1] <- m.motif_asbs[motif, 1] + 1  
        }
        if(n.B73 == 0 & n.Mo17 >= 1){
          m.motif_asbs[motif, 3] <- m.motif_asbs[motif, 3] + 1
          m.motif_asbs[motif, 2] <- m.motif_asbs[motif, 2] + 1
        }
      }
    
      if(df.bQTLs$POSTfreq[i] < 0.5){
        if(n.B73 == 0 & n.Mo17 >= 1){
          # only for correct direction 
          v.pos <- s.asb.pos - start(v.matches.B73[[i]]) + 1
          m.motif_variation[motif,v.pos] <- m.motif_variation[motif,v.pos] + 1
          
          m.motif_asbs[motif, 4] <- m.motif_asbs[motif, 4] + 1
          m.motif_asbs[motif, 2] <- m.motif_asbs[motif, 2] + 1
        }
        if(n.B73 >= 1 & n.Mo17 == 0){
          m.motif_asbs[motif, 3] <- m.motif_asbs[motif, 3] + 1
          m.motif_asbs[motif, 1] <- m.motif_asbs[motif, 1] + 1
        }
      }
    }
    
  }
  

  write.csv(as.data.frame(m.motif_variation), "motif_base_6143_ASBs_variation.csv")
  print(m.motif_variation)
  
  m.motif_variation <- round(m.motif_variation / rowSums(m.motif_variation) * 100, 2) 
  print(m.motif_variation)
  
  write.csv(as.data.frame(m.motif_variation), "motif_base_6143_ASBs_variation_percentage.csv")
  
  df.motif_asbs <- as.data.frame(m.motif_asbs)
  df.motif_asbs["percentage dir. correct"] <- round(df.motif_asbs$`dir. correct` / df.motif_asbs$`total (non-ambiguous)`* 100, 2) 
  df.motif_asbs["motif-seq"] <- v.motifs
  df.motif_asbs["motif-name"] <- names(v.motifs)
  print(df.motif_asbs)
  
  write.csv(df.motif_asbs, "motif_directionality_6143_ASBs.csv")
  
  write.csv(df.bQTLs, "6143_ASBs_w_motif_annotation.csv", quote = F, row.names = F)
  
  #saveRDS(df.bQTLs, paste(folder_tmp, "317094_bgSNPs_w_genomic_location.rds", sep = "/"))

  
}