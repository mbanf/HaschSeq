bQTL_on_ArabidopsisHomolog_geneexpression_evaluation <- function(){
  
  l.df.overlap_zm_vs_ath = vector(mode="list", length = 4)
  l.df.overlap_zm_vs_ath[[1]] <- read.csv("data/arabidopsis_overlap_genelists/AtBT_AtBR_exlusive_to_arab.txt", stringsAsFactors = FALSE, header = F, skip = 1, sep = "\t")
  l.df.overlap_zm_vs_ath[[2]] <- read.csv("data/arabidopsis_overlap_genelists/ZmBT_ZmBR_AtBT_AtBR_overlap.txt", stringsAsFactors = FALSE,header = F, skip = 1, sep = "\t")
  l.df.overlap_zm_vs_ath[[3]] <- read.csv("data/arabidopsis_overlap_genelists/ZmBT_ZmBR_exlusive_to_maize.txt", stringsAsFactors = FALSE, header = F, skip = 1,sep = "\t")
  l.df.overlap_zm_vs_ath[[4]] <- read.csv("data/arabidopsis_overlap_genelists/Zm_noBT_BR_vsAt_noBT_BR_overlap.txt", stringsAsFactors = FALSE, header = F, skip = 1, sep = "\t")
  
  
  # population
  path.gene_orthologs <- "data/GeneAnnotation/Phytozome13_anno_1_1.txt"
  df.gene_orthologs = read.table(path.gene_orthologs, header=T, stringsAsFactors = FALSE, sep = "\t", fill = T)
  
  
  path.rnaseq.down_regulated <- "data/expression//BR_repressed_adjpvalue_0.05&l2FC_0.5.txt"
  path.rnaseq.up_regulated <- "data/expression/BR_induced_adjpvalue_0.05&l2FC_0.5.txt"
  
  df.rnaseq.down_regulated <- read.table(path.rnaseq.down_regulated, header = TRUE, sep = "\t", fill = TRUE)
  df.rnaseq.up_regulated <- read.table(path.rnaseq.up_regulated, header = TRUE, sep = "\t", fill = TRUE)
  
  df.bQTL_RNAseq <- df.rnaseq.up_regulated[,c(1,3)]
  df.bQTL_RNAseq <- rbind(df.bQTL_RNAseq, df.rnaseq.down_regulated[,c(1,3)])
  
  names(df.bQTL_RNAseq)[1] <- "gene.ID"
  
  df.bQTL_RNAseq <- merge(df.bQTL_RNAseq, df.gene_orthologs, by = "gene.ID")
  df.expression_zmays = df.bQTL_RNAseq
  
  
  v.gns.overlap = l.df.overlap_zm_vs_ath[[2]][,1]
  
  #df.overlap_zm_vs_ath <- read.csv("data/arabidopsis_overlap_genelists/ZmnoBTBRvsAtnoBTBR.csv", stringsAsFactors = FALSE)
  #df.overlap_zm_vs_ath <- read.csv("data/arabidopsis_overlap_genelists/ZmBTBRvsAtBTBR.csv", stringsAsFactors = FALSE)
  
  df.expression_ath <- read.csv("data/arabidopsis_overlap_genelists/Ath_BR_elife_03031_fig3_data1_v2.txt", stringsAsFactors = FALSE, sep = "\t")
  
  
  # df.expression_zmays <- read.csv("data/arabidopsis_overlap_genelists/AllmaizegeneswithBR_RNAseq.txt", stringsAsFactors = FALSE, sep = "\t")
  #v.gns.overlap <- as.character(df.overlap_zm_vs_ath[,1])
  #v.gns.overlap <- v.gns.overlap[v.gns.overlap != ""]
  
  # expression arabidopsis 
  df.expression_zmays.overlap <- subset(df.expression_zmays, df.expression_zmays$Ara_gene.ID %in% v.gns.overlap)
  v.gns.overlap <- unique(df.expression_zmays.overlap$Ara_gene.ID)
  
  
  m.expr_directionality <- matrix(NA, nrow = length(v.gns.overlap), ncol = 3, dimnames = list(v.gns.overlap, c("Ath", "Zm1", "Zm2")))
  
  v.directionality = numeric(length(v.gns.overlap))
  
  for(i in 1:length(v.gns.overlap)){
    
    df.expression_zmays.overlap.i <- subset(df.expression_zmays.overlap, df.expression_zmays.overlap$Ara_gene.ID == v.gns.overlap[i])
    df.expression_ath.i <- subset(df.expression_ath, df.expression_ath$Gene_ID_Ath == v.gns.overlap[i])
    
    # if non duplicate 
    if(nrow(df.expression_zmays.overlap.i) == 1 & nrow(df.expression_ath.i) == 1){
      
      m.expr_directionality[i,1] <- as.numeric(df.expression_ath.i$log2FC_WT_.BL_vs_WT.BL)
      m.expr_directionality[i,2] <- as.numeric(df.expression_zmays.overlap.i$log2FoldChange)
      
      if(df.expression_zmays.overlap.i$log2FoldChange > 0 & df.expression_ath.i$log2FC_WT_.BL_vs_WT.BL > 0){
        v.directionality[i] <- 1
      }else if(df.expression_zmays.overlap.i$log2FoldChange < 0 & df.expression_ath.i$log2FC_WT_.BL_vs_WT.BL < 0){
        v.directionality[i] <- 1
      }else if(df.expression_zmays.overlap.i$log2FoldChange > 0 & df.expression_ath.i$log2FC_WT_.BL_vs_WT.BL < 0){
        v.directionality[i] <- - 1
      }else if(df.expression_zmays.overlap.i$log2FoldChange < 0 & df.expression_ath.i$log2FC_WT_.BL_vs_WT.BL > 0){
        v.directionality[i] <- - 1
      }
      
    }else{
      #print(df.expression_zmays.overlap.i$log2FoldChange.x)
      m.expr_directionality[i,1] <- as.numeric(df.expression_ath.i$log2FC_WT_.BL_vs_WT.BL)[1]
      m.expr_directionality[i,2] <- as.numeric(df.expression_zmays.overlap.i$log2FoldChange[1])
      m.expr_directionality[i,3] <- as.numeric(df.expression_zmays.overlap.i$log2FoldChange[2])
    }
    
  }
  
  
  m.expr_directionality[m.expr_directionality < -3.5] <- -3.5
  m.expr_directionality[m.expr_directionality > 3.5] <- 3.5
  
  library(gplots)
  library(pheatmap)
  library(RColorBrewer)
  #p <- pheatmap((m.expr_directionality), cluster_rows = TRUE,cluster_cols = FALSE, labels_row = NULL, labels_col = NULL, legend_labels = c("< -3", "-2", "-1", "0", "1", "2", "> 3.5"))
  m.expr_directionality.dual <- m.expr_directionality[which(m.expr_directionality[,3] == 0),]
  m.expr_directionality.dual <- m.expr_directionality.dual[,1:2]
  
  p <- pheatmap((m.expr_directionality.dual), cluster_rows = TRUE,cluster_cols = FALSE, color = greenred(100), legend_labels = c("< -3", "-2", "-1", "0", "1", "2", "> 3.5"))
  
  m.expr_directionality.triple <- m.expr_directionality[which(m.expr_directionality[,3] > 0),]
  #m.expr_directionality.dual <- m.expr_directionality.dual[,1:2]
  
  p <- pheatmap((m.expr_directionality), cluster_rows = TRUE,cluster_cols = FALSE,  color = greenred(100), legend_labels = c("< -3", "-2", "-1", "0", "1", "2", "> 3.5"))
  
  save_pheatmap_pdf(p,  paste(folder_output, "expression/F0.txt", sep = ""),  width=10, height=10)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #### 
  
  v.Arahomologs_ZmBZR1_targets_nonBRreg <- toupper(as.character(read.table("data/expression/Arahomologs_ZmBZR1_targets_nonBRreg.txt", sep = "\t", header = FALSE)[,1]))
  v.Arahomologs_ZmBZR1_targets_BRreg <- toupper(as.character(read.table("data/expression/Arahomologs_ZmBZR1_targets_BRreg.txt", sep = "\t", header = FALSE)[,1]))
  v.Arahomologs_non_ZmBZR1_targets_BRreg <- toupper(as.character(read.table("data/expression/Arahomologs_non_ZmBZR1_targets_BRreg.txt", sep = "\t", header = FALSE)[,1]))
  
  
  # df.bQTL_gene_partitioning <- read.csv("output/df.bQTL_gene_partitioning_B73chip.csv")
  df.bQTL_gene_partitioning$Arabidopsis_ortholog <- gsub("\\..*", "", df.bQTL_gene_partitioning$Arabidopsis_ortholog) # Arabidopsis 
  
  # v.Arabidopsis_ortholog <- as.character(read.table("data/expression/ZmBZR1_targets_Ara_homologs_virtualplant.txt", sep = "\t", header = FALSE)[,1])
  
  length(v.Arahomologs_ZmBZR1_targets_nonBRreg)
  length(v.Arahomologs_ZmBZR1_targets_BRreg)
  length(v.Arahomologs_non_ZmBZR1_targets_BRreg)
  
  
  length(intersect(v.Arahomologs_ZmBZR1_targets_nonBRreg, v.Arahomologs_non_ZmBZR1_targets_BRreg))
  length(intersect(v.Arahomologs_ZmBZR1_targets_nonBRreg, v.Arahomologs_ZmBZR1_targets_BRreg))
  length(v.Arahomologs_ZmBZR1_targets_nonBRreg)
  
  
  v.AT_BR_regulated <- as.character(read.table("data/expression/At_BR_regulated_genes.txt", sep = "\t", header = FALSE)[,1])
  
  length(intersect(v.AT_BR_regulated, df.At_BZR1_ChipSeqTargets$locus))
  length(v.AT_BR_regulated)
  
  
  length(intersect(df.At_BZR1_ChipSeqTargets$locus, v.Arahomologs_ZmBZR1_targets_nonBRreg))
  length(intersect(intersect(v.AT_BR_regulated, df.At_BZR1_ChipSeqTargets$locus), v.Arahomologs_ZmBZR1_targets_BRreg))
  length(intersect(v.AT_BR_regulated, v.Arahomologs_non_ZmBZR1_targets_BRreg))
  
  ##### older data 
  
  
  message("Venn diagram overlap analysis - Arabidopsis comparison")
  
  library(VennDiagram)
  library(xlsx)
  
  df.homologies <- read.table("datasets_paper/GeneAnnotation/Zmays_284_6a.annotation_info.txt", 
                              header = FALSE, sep ="\t", quote = "", stringsAsFactors = FALSE, fill = TRUE)
  names(df.homologies) <- c("zmays", "ath")
  
  ## read BR genes
  v.genes <- as.character(read.csv("BR_regulated_genes_talk.csv", header = FALSE, quote = "", stringsAsFactors = FALSE)[,1])
  
  df.homologies.subset <- subset(df.homologies, df.homologies$zmays %in% v.genes)
  df.homologies.subset <- subset(df.homologies.subset, df.homologies.subset$ath != "")
  df.homologies.subset["locus"] <- gsub("\\..*","", df.homologies.subset$ath)
  
  
  df.homologies.mapping <- subset(df.homologies.subset, df.homologies.subset$locus %in% df.At_BZR1_ChIPchipTargets$locus)
  nrow(df.homologies.mapping)
  
  df.homologies.mapping <- subset(df.homologies.subset, df.homologies.subset$locus %in% df.At_BZR1_ChipSeqTargets$locus)
  nrow(df.homologies.mapping)
  
  length(unique(intersect(intersect(df.At_BZR1_ChipSeqTargets$locus, df.At_BZR1_ChIPchipTargets$locus), df.homologies.subset$locus)))
  
  a1 <- length(unique(df.At_BZR1_ChIPchipTargets$locus))
  a2 <- length(unique(df.At_BZR1_ChipSeqTargets$locus))
  a3 <- length(unique(df.homologies.subset$locus)) 
  
  a12 <- length(intersect(unique(df.At_BZR1_ChIPchipTargets$locus), unique(df.At_BZR1_ChipSeqTargets$locus)))
  a23 <- length(intersect(unique(df.At_BZR1_ChipSeqTargets$locus), unique(df.homologies.subset$locus)))
  a13 <- length(intersect(unique(df.At_BZR1_ChIPchipTargets$locus), unique(df.homologies.subset$locus)))
  a123 <- length(intersect(intersect(unique(df.At_BZR1_ChIPchipTargets$locus), unique(df.At_BZR1_ChipSeqTargets$locus)), unique(df.homologies.subset$locus)))
  
  
  lst.sets <- vector(mode = "list", length = 7)
  lst.sets[[1]] <- intersect(intersect(unique(df.At_BZR1_ChIPchipTargets$locus), unique(df.At_BZR1_ChipSeqTargets$locus)), unique(df.homologies.subset$locus))
  
  lst.sets[[2]] <- intersect(unique(df.At_BZR1_ChIPchipTargets$locus), unique(df.At_BZR1_ChipSeqTargets$locus)) 
  lst.sets[[3]] <- intersect(unique(df.At_BZR1_ChipSeqTargets$locus), unique(df.homologies.subset$locus))
  lst.sets[[4]] <- intersect(unique(df.At_BZR1_ChIPchipTargets$locus), unique(df.homologies.subset$locus))
  lst.sets[[2]] <- lst.sets[[2]][!lst.sets[[2]] %in% lst.sets[[1]]]
  lst.sets[[3]] <- lst.sets[[3]][!lst.sets[[3]] %in% lst.sets[[1]]]
  lst.sets[[4]] <- lst.sets[[4]][!lst.sets[[4]] %in% lst.sets[[1]]]
  
  lst.sets[[5]] <- unique(df.At_BZR1_ChIPchipTargets$locus)
  lst.sets[[6]] <- unique(df.At_BZR1_ChipSeqTargets$locus)
  lst.sets[[7]] <- unique(df.homologies.subset$locus)
  
  lst.sets[[5]] <- lst.sets[[5]][!lst.sets[[5]] %in% c(unique(df.At_BZR1_ChipSeqTargets$locus), unique(df.homologies.subset$locus))]
  lst.sets[[6]] <- lst.sets[[6]][!lst.sets[[6]] %in% c(unique(df.At_BZR1_ChIPchipTargets$locus), unique(df.homologies.subset$locus))]
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
  save.xlsx("GeneSets.xlsx", a123_466, a12_1038, a23_525, a13_420, a1_1487, a2_2272, a3_2784)
  
  
  ## remark: do venn diagramm
  grid.newpage()
  draw.triple.venn(area1 = a1, area2 = a2, area3 = a3, n12 = a12, n23 = a23, n13 = a13, cex = 1.5,  col = 1, cat.cex = 1.5, SetNames=c( "A", "B","A", "B","A", "B","A"),
                   n123 = a123, category = c("AtBZR1 Chip", "AtBZR1 Seq", "ZmBZR1 Seq"), lty = 1, 
                   fill = c("blue", "red", "green"))
  
  
}