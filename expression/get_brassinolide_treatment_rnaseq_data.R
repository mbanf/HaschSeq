get_brassinolide_treatment_rnaseq_data <- function(path.rnaseq.down_regulated = "data/expression/BR_repressed_adjpvalue_0.05&l2FC_0.5.txt",
                                                   path.rnaseq.up_regulated = "data/expression/BR_induced_adjpvalue_0.05&l2FC_0.5.txt"){

  df.rnaseq.down_regulated <- read.table(path.rnaseq.down_regulated, header = TRUE, sep = "\t", fill = TRUE)
  df.rnaseq.up_regulated <- read.table(path.rnaseq.up_regulated, header = TRUE, sep = "\t", fill = TRUE)
  
  df.bQTL_RNAseq <- df.rnaseq.up_regulated[,c(1,3)]
  df.bQTL_RNAseq <- rbind(df.bQTL_RNAseq, df.rnaseq.down_regulated[,c(1,3)])
  
  df.bQTL_RNAseq["mode"] <- "down"
  df.bQTL_RNAseq$mode[1:nrow(df.rnaseq.up_regulated)] <- "up"
  names(df.bQTL_RNAseq)[1] <- "gene.ID"

  df.bQTL_RNAseq  
}