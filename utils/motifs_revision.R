predefined_motif_analysis <- function(df.bQTLs, 
                                      v.motifs, 

                                      s.env_bps = 5,
                                      n.chromosomes = 10,
                     
                                      path.genomic_sequences = "data/revision/ref_B73Mo17.fasta",
                                      ){
  
  # To identify ASBs which may be explained by motif variation (Fig. 2h), we extracted the
  # +/- 5 bp of the high affinity BZR1 bound allele surrounding ASBs. Using R we scanned
  # those 11 bp fragments for canonical BRRE (CGTG[T/C]G, C[G/A]CACG), allowing a
  # single base pair mismatch outside the core motif (CGTG), or G-box (CACGTG) motifs
  # and determined ASBs where the SNP changed a BRRE or G-box motif into an altered
  # (non BRRE or G-box) motif. 
  
  message("Extracting +/- 5 bp of the high affinity BZR1 bound allele surrounding ASBs and scanning for canonical BRRE (CGTG[T/C]G, C[G/A]CACG), \n allowing a single base pair mismatch outside the core motif (CGTG), or G-box (CACGTG) motifs")
  
  # TODO: invert to identify the sequences on the opposite strand to see creation and deletion of motifs?
  # TODO: directionality - ASB to one side, by motif has been removed - wrong directinalioty 
  
  asb_environment <- extract_sequences(df.bQTLs,
                                       genome_sequences = readDNAStringSet(path.genomic_sequences),
                                       s.env_bps = s.env_bps,
                                       n.chromosomes = 10)

  
  # Fraction of ASBs in motifs (and nucleotide position)
  m.motif_variation <- matrix(0, nrow = length(v.motifs), ncol = s.env_bps + 1) 
  rownames(m.motif_variation) <- v.motifs
  
  df.bQTLs[,v.motifs] <- F
  s.asb.pos <- s.env_bps + 1
  
  for(motif in v.motifs){
  
    v.matches <- vmatchPattern(motif, asb_environment, fixed = TRUE)
    idx.bQTLs <- which(as.numeric(lapply(v.matches, length)) > 0)
    df.bQTLs[idx.bQTLs, motif] <- T
    
    for(i in idx.bQTLs){
      v.pos <- s.asb.pos - start(v.matches[[i]]) + 1
      m.motif_variation[motif,v.pos] <- m.motif_variation[motif,v.pos] + 1
      
      # check directionality 
      df.bQTLs$POSTfreq[i]
      
      
    }
    
  }
  
  m.motif_variation <- round(m.motif_variation / nrow(df.bQTLs) * 100, 2) 
  print(m.motif_variation)
  
  
  
  m.motif_asbs <- matrix(0, nrow = length(v.motifs), ncol = 2) 
  rownames(m.motif_asbs) <- v.motifs
  colnames(m.motif_asbs) <- c("no", "yes")
  for(motif in v.motifs){
    m.motif_asbs[motif, ] <- as.numeric(table(df.bQTLs[,motif]))
  }
  print(m.motif_asbs)
  
  n.asb_motifs <- length(which(as.numeric(rowSums(df.bQTLs[,v.motifs])) > 0))
  print(n.asb_motifs)
  
  
  # directionality - motif created resulting in higher binding affinity - per motif analysis 
  for(motif in v.motifs){
    
  }
  
  # change of position but still 
  
  
  
  
  
}