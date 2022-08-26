
br_vs_non_br_regulated <- function(filename = "data/arabidopsis_overlap_genelists/ZmvsAth_1_to_1.txt",
                                   folder_tmp ="tmp/"){
  
  
  df.ZM_Ath_overlaps <- read.table(filename, stringsAsFactors = F, fill = T, header = T)
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

  write.csv(paste(folder_tmp, df.ZM_Ath_overlaps, "df.ZM_Ath_overlaps.csv", sep ="/"))
  
  ####
  
  v.categories <- c("non-regulated BR targets", "regulated BR targets", "regulated non-BR targets", "non-regulated non-BR targets")
  
  m.heatmap.fc <- matrix(0, nrow = 4, ncol = 4, dimnames = list(v.categories, v.categories))
  m.heatmap.pval <- matrix(1, nrow = 4, ncol = 4, dimnames = list(v.categories, v.categories))
  
  m.heatmap.fc[1:16] <- df.ZM_Ath_overlaps$fc
  m.heatmap.fc <- t(m.heatmap.fc)
  
  m.heatmap.pval[1:16] <- df.ZM_Ath_overlaps$p.val
  m.heatmap.pval <- t(m.heatmap.pval)
  
  library(pheatmap)
  p <- pheatmap((m.heatmap.fc), cluster_rows = FALSE,cluster_cols = FALSE, labels_row = NULL, labels_col = NULL)
  
}



arabidopsis_chipseq_venn_overlap <- function(df.ChipSeq.gene_partitioning){
  
  message("paper only - overlap analysis")
  
  library(VennDiagram)
  
  df.At_BZR1_ChipSeqTargets <- read.table("data/ArabidopsisScriptsAndDatasets/At_BZR1_ChipSeqTargets.txt", header = TRUE, sep ="\t", quote = "", stringsAsFactors = FALSE)
  names(df.At_BZR1_ChipSeqTargets) <- "locus"
  
  # ath chIP chip - light
  df.At_BZR1_ChIPchipTargets <- read.table("data/ArabidopsisScriptsAndDatasets/At_BZR1_ChIPchip.txt", header = TRUE, sep ="\t", quote = "", stringsAsFactors = FALSE)
  
  df.ChipSeq.gene_partitioning.subset <- subset(df.ChipSeq.gene_partitioning, !is.na(df.ChipSeq.gene_partitioning$Arabidopsis_ortholog))
  df.ChipSeq.gene_partitioning.subset$Arabidopsis_ortholog <- gsub("\\..*","", df.ChipSeq.gene_partitioning.subset$Arabidopsis_ortholog)
  
  ######## 
  
  a1 <- length(unique(df.At_BZR1_ChIPchipTargets$locus))
  a2 <- length(unique(df.At_BZR1_ChipSeqTargets$locus))
  a3 <- length(unique(df.ChipSeq.gene_partitioning.subset$Arabidopsis_ortholog)) 
  
  a12 <- length(intersect(unique(df.At_BZR1_ChIPchipTargets$locus), unique(df.At_BZR1_ChipSeqTargets$locus)))
  a23 <- length(intersect(unique(df.At_BZR1_ChipSeqTargets$locus), unique(df.ChipSeq.gene_partitioning.subset$Arabidopsis_ortholog)))
  a13 <- length(intersect(unique(df.At_BZR1_ChIPchipTargets$locus), unique(df.ChipSeq.gene_partitioning.subset$Arabidopsis_ortholog)))
  a123 <- length(intersect(intersect(unique(df.At_BZR1_ChIPchipTargets$locus), unique(df.At_BZR1_ChipSeqTargets$locus)), unique(df.ChipSeq.gene_partitioning.subset$Arabidopsis_ortholog)))
  
  
  lst.sets <- vector(mode = "list", length = 7)
  lst.sets[[1]] <- intersect(intersect(unique(df.At_BZR1_ChIPchipTargets$locus), unique(df.At_BZR1_ChipSeqTargets$locus)), unique(df.ChipSeq.gene_partitioning.subset$Arabidopsis_ortholog))
  
  lst.sets[[2]] <- intersect(unique(df.At_BZR1_ChIPchipTargets$locus), unique(df.At_BZR1_ChipSeqTargets$locus)) 
  lst.sets[[3]] <- intersect(unique(df.At_BZR1_ChipSeqTargets$locus), unique(df.ChipSeq.gene_partitioning.subset$Arabidopsis_ortholog))
  lst.sets[[4]] <- intersect(unique(df.At_BZR1_ChIPchipTargets$locus), unique(df.ChipSeq.gene_partitioning.subset$Arabidopsis_ortholog))
  lst.sets[[2]] <- lst.sets[[2]][!lst.sets[[2]] %in% lst.sets[[1]]]
  lst.sets[[3]] <- lst.sets[[3]][!lst.sets[[3]] %in% lst.sets[[1]]]
  lst.sets[[4]] <- lst.sets[[4]][!lst.sets[[4]] %in% lst.sets[[1]]]
  
  lst.sets[[5]] <- unique(df.At_BZR1_ChIPchipTargets$locus)
  lst.sets[[6]] <- unique(df.At_BZR1_ChipSeqTargets$locus)
  lst.sets[[7]] <- unique(df.ChipSeq.gene_partitioning.subset$Arabidopsis_ortholog)
  
  lst.sets[[5]] <- lst.sets[[5]][!lst.sets[[5]] %in% c(unique(df.At_BZR1_ChipSeqTargets$locus), unique(df.ChipSeq.gene_partitioning.subset$Arabidopsis_ortholog))]
  lst.sets[[6]] <- lst.sets[[6]][!lst.sets[[6]] %in% c(unique(df.At_BZR1_ChIPchipTargets$locus), unique(df.ChipSeq.gene_partitioning.subset$Arabidopsis_ortholog))]
  lst.sets[[7]] <- lst.sets[[7]][!lst.sets[[7]] %in% c(unique(df.At_BZR1_ChIPchipTargets$locus), unique(df.At_BZR1_ChipSeqTargets$locus))]
  
  names(lst.sets) <- c("a123_466", "a12_1038","a23_525","a13_420", "a1_1487","a2_2272","a3_2784")
  
  
  # write genes 
  a123_466 <- lst.sets[[1]]
  a12_1038 <- lst.sets[[2]]
  a23_525 <- lst.sets[[3]]
  a13_420 <- lst.sets[[4]]
  a1_1487 <- lst.sets[[5]]
  a2_2272 <- lst.sets[[6]]
  a3_2784 <- lst.sets[[7]]
  
  # library(xlsx)
  # save.xlsx("GeneSets.xlsx", a123_466, a12_1038, a23_525, a13_420, a1_1487, a2_2272, a3_2784)
  
  
  ## remark: do venn diagramm
  grid.newpage()
  draw.triple.venn(area1 = a1, area2 = a2, area3 = a3, n12 = a12, n23 = a23, n13 = a13, cex = 1.5,  col = 1, cat.cex = 1.5, SetNames=c( "A", "B","A", "B","A", "B","A"),
                   n123 = a123, category = c("AtBZR1 ChipChip", "AtBZR1 ChipSeq", "ZmBZR1 Seq"), lty = 1, 
                   fill = c("blue", "red", "green"))
  
  
}


