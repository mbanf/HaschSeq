

evaluating_BZR1_core_motif <- function(df.ASBs, core_motif = "CGTG"){
  
  core_motif = "CGTG"
  s.motif_offset = 5
  
  df.ASBs.core_motifs = c()
  
  
  for(i in 1:n.chromosomes){
    
    df.ASBs.i = subset(df.ASBs, df.ASBs$contig == i)
    
    df.ASBs.i["B73_sequence"] = ""
    df.ASBs.i["MO17_sequence"] = ""
    df.ASBs.i["B73_pos_core_motif"] = NA
    df.ASBs.i["MO17_pos_core_motif"] = NA
    
    
    print(paste("processing chromosome ", i))
    
    genome.reference <- DNAString(genome[[i]])
    genome.mutant    <- l.genome.mutant[[i]]
    
    for(j in 1:nrow(df.ASBs.i)){
      
      pos.ASB.j = df.ASBs.i$position[j]
      
      motif.reference <- as.character(subseq(genome.reference, start = pos.ASB.j - s.motif_offset, end = pos.ASB.j + s.motif_offset))    
      motif.mutant <- as.character(subseq(genome.mutant, start = pos.ASB.j - s.motif_offset, end = pos.ASB.j + s.motif_offset))    
      
      df.ASBs.i$B73_sequence[j] = motif.reference
      df.ASBs.i$MO17_sequence[j] = motif.mutant
      
      res = unlist(gregexpr(pattern = core_motif, motif.reference))
      if(res != -1){
        res = res - s.motif_offset - 1
        df.ASBs.i$B73_pos_core_motif[j] = paste(res, collapse = ",")
      }
      
      res = unlist(gregexpr(pattern = core_motif, motif.mutant))
      if(res != -1){
        res = res - s.motif_offset - 1
        df.ASBs.i$MO17_pos_core_motif[j] = paste(res, collapse = ",")
      }
      
      
    }
    
    df.ASBs.core_motifs = rbind(df.ASBs.tmp, df.ASBs.i)
    
  }
  
  
  write.table(df.ASBs.core_motifs, paste(folder_output, "/SX3.txt", sep = ""), quote = FALSE, row.names = FALSE, sep ="\t")
  # write.csv2(df.ASBs.tmp, "A:/junkDNA.ai/df.ASBs.csv",row.names = F)
  
  
  ## NAS
  
  
}



