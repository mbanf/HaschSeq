
br_vs_non_br_regulated <- function(){
  
  
  df.ZM_Ath_overlaps <- read.csv("datasets_paper/ZmvsAth_1_to_1.csv")
  df.ZM_Ath_overlaps["p.val"] <- 1
  df.ZM_Ath_overlaps["fc"] <- 0
  
  # overla analysis 
  for(i in 1:nrow(df.ZM_Ath_overlaps)){
    
    ## domain genes in target genes
    hitInSample = n_A_B = df.ZM_Ath_overlaps$Overlap[i]
    sampleSize = n_A = df.ZM_Ath_overlaps$unique_Ath[i] +  hitInSample
    ## domain genes in genome
    hitInPop = n_B =  df.ZM_Ath_overlaps$unique_Zm[i] + hitInSample
    popSize = n_C =  df.ZM_Ath_overlaps$total[1]
    
    failInPop = n_C-n_B
    
    fc <- (n_A_B / n_A) / (n_B / n_C)
    pval = fisher.test(matrix(c(hitInSample, hitInPop-hitInSample, sampleSize-hitInSample, failInPop-sampleSize +hitInSample), 2, 2), alternative='greater')$p.value; 
    
    df.ZM_Ath_overlaps$fc[i] <- fc
    df.ZM_Ath_overlaps$p.val[i] <- pval
    
  }
  
  ####
  
  v.categories <- c("BR targets but not regulated", "BR targets and regulated", "no BR targets but regulated", "no BR targets and not regulated")
  
  m.heatmap.fc <- matrix(0, nrow = 4, ncol = 4, dimnames = list(v.categories, v.categories))
  m.heatmap.pval <- matrix(1, nrow = 4, ncol = 4, dimnames = list(v.categories, v.categories))
  
  m.heatmap.fc[1:16] <- df.ZM_Ath_overlaps$fc
  m.heatmap.fc <- t(m.heatmap.fc)
  
  m.heatmap.pval[1:16] <- df.ZM_Ath_overlaps$p.val
  m.heatmap.pval <- t(m.heatmap.pval)
  
  library(pheatmap)
  p <- pheatmap((m.heatmap.fc), cluster_rows = FALSE,cluster_cols = FALSE, labels_row = NULL, labels_col = NULL)
  
  
  
}

