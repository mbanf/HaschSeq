directionality_and_snp_distribution <- function(l.motif_analysis.postprocessed){
  
  message("figure 2 - directionality of ASB")

  
  df.motif_directionality <- data.frame(motif = motifs,  yes = numeric(length(motifs)),  no = numeric(length(motifs)))
  
  for(m in 1:length(motifs)){
    test <- l.motif_analysis.postprocessed[[1]] 
    tmp <- subset(test, test$motif.mutant ==  motifs[m] & test$unique == "mutant")
    tmp <- rbind(tmp, subset(test, test$motif.ref ==  motifs[m] & test$unique == "reference"))
    
    df.motif_directionality$yes[m] <- table(tmp$direction)["yes"]
    df.motif_directionality$no[m] <- table(tmp$direction)["no"]
  }
  
  saveRDS(df.motif_directionality, paste(folder_tmp, "df.motif_directionality.rds", sep = "/"))
  
  
  message("figure 2 - motif analysis - snp distribution")
  
  df.motif_nucleotidePositionPartitions <- data.frame(motif = character(),
                                                      type = character(),
                                                      inPeak = character(),
                                                      nucleotidePos_1 = numeric(),
                                                      nucleotidePos_2 = numeric(),
                                                      nucleotidePos_3 = numeric(),
                                                      nucleotidePos_4 = numeric(),
                                                      nucleotidePos_5 = numeric(),
                                                      nucleotidePos_6 = numeric())
  
  
  v.inPeak <- which(df.motif_analysis.postprocessed$isMotifInPeak == "no")
  v.direction <- which(df.motif_analysis.postprocessed$direction == "no")
  
  df.motif_analysis.postprocessed <- l.motif_analysis.postprocessed[[1]]
  
  v.inPeak <- which(df.motif_analysis.postprocessed$isMotifInPeak == "yes")
  v.direction <- which(df.motif_analysis.postprocessed$direction == "yes")
  
  length(intersect(v.inPeak,v.direction ))
  length(v.inPeak)
  length(v.direction)
  
  set_sign <- c("ASBs", "bgSNPs")
  set_peak <- c("yes", "no")
  
  peak <- c("yes", "no")
  
  for(s in 1:2){
    for(p in 1:2){

      df.motif_analysis.postprocessed <- l.motif_analysis.postprocessed[[s]]
      df.motif_analysis.postprocessed <- subset(df.motif_analysis.postprocessed, df.motif_analysis.postprocessed$isMotifInPeak == peak[p])
      
      for(l in 1:length(motifs)){
        
        test.l <- subset(df.motif_analysis.postprocessed, df.motif_analysis.postprocessed$motif.ref == motifs[l] | df.motif_analysis.postprocessed$motif.mutant == motifs[l])
        snp.pos <- (test.l$position - test.l$pos.motif) + 1
        
        # directionality test 
        if(FALSE){
          test.ld <- subset(test.l, test.l$direction == "yes")
          snp.pos <- (test.ld$position - test.ld$pos.motif) + 1
        }
        
        
        snp.pos = table(snp.pos)
        
        nucleotidePos = c(0,0,0,0,0,0)
        
        for(i in 1:length(nucleotidePos)){
          
          if(as.character(i) %in% names(snp.pos)){
            
            nucleotidePos[i] <- snp.pos[as.character(i)]
            
          }
        }
        
        
        
        # print(table(snp.pos)) # all
        df.motif_nucleotidePositionPartitions <- rbind(df.motif_nucleotidePositionPartitions, 
                                                       data.frame(motif = motifs[l],
                                                                  type = set_sign[s],
                                                                  isMotifInPeak = set_peak[p],
                                                                  nucleotidePos_1 = nucleotidePos[1],
                                                                  nucleotidePos_2 = nucleotidePos[2],
                                                                  nucleotidePos_3 = nucleotidePos[3],
                                                                  nucleotidePos_4 = nucleotidePos[4],
                                                                  nucleotidePos_5 = nucleotidePos[5],
                                                                  nucleotidePos_6 = nucleotidePos[6]))
        
      }
      
    }
    
  }
  
  message("...finished")
  
  saveRDS(df.motif_nucleotidePositionPartitions, paste(folder_tmp, "df.motif_nucleotidePositionPartitions.rds", sep = "/"))

}