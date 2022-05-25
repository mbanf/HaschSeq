
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






BR_up_and_Down <- function(df.bQTL_gene_partitioning){
  #### TODO: create table S1 #### 
  message("List of BR responsive genes in Maize determined by RNA-seq - Table S1")
  
  # BES 1 
  # df.rnaseq.ATBES1_targets <- read.table(path.rnaseq.ATBES1_targets, header = FALSE, sep = "\t", fill = TRUE)
  # v.rnaseq.ATBES1_targets <- unique(c(as.character(df.rnaseq.ATBES1_targets[,1]), as.character(df.rnaseq.ATBES1_targets[,2]), as.character(df.rnaseq.ATBES1_targets[,3]),
  #                                     as.character(df.rnaseq.ATBES1_targets[,4]), as.character(df.rnaseq.ATBES1_targets[,5]), as.character(df.rnaseq.ATBES1_targets[,6]),
  #                                     as.character(df.rnaseq.ATBES1_targets[,7]), as.character(df.rnaseq.ATBES1_targets[,8])))
  # 
  # 
  # print(length(v.rnaseq.ATBES1_targets))
  
  # BR regulated 
  path.rnaseq.down_regulated <- "data/expression/BR_repressed_adjpvalue_0.05&l2FC_0.5.txt"
  path.rnaseq.up_regulated <- "data/expression/BR_induced_adjpvalue_0.05&l2FC_0.5.txt"
  path.rnaseq.ATBES1_targets <- "data/expression/AtBES1_targets.txt"
  
  df.rnaseq.down_regulated <- read.table(path.rnaseq.down_regulated, header = TRUE, sep = "\t", fill = TRUE)
  df.rnaseq.up_regulated <- read.table(path.rnaseq.up_regulated, header = TRUE, sep = "\t", fill = TRUE)
  
  df.bQTL_RNAseq <- df.rnaseq.up_regulated[,c(1,3)]
  df.bQTL_RNAseq <- rbind(df.bQTL_RNAseq, df.rnaseq.down_regulated[,c(1,3)])
  
  df.bQTL_RNAseq["mode"] <- "down"
  df.bQTL_RNAseq$mode[1:nrow(df.rnaseq.up_regulated)] <- "up"
  names(df.bQTL_RNAseq)[1] <- "gene.ID"
  #df.bQTL_RNAseq <- unique(df.bQTL_RNAseq[,c("gene.ID", "diffExp", "mode")])
  message("Number of BR responsive genes: ", nrow(df.bQTL_RNAseq) )
  write.table(df.bQTL_RNAseq, paste(folder_output, "/S1.txt", sep = ""), quote = FALSE, row.names = FALSE, sep = "\t")
  
  v.gns.total <- unique(df.bQTL_gene_partitioning$gene.ID)
  
  # todo: select only genes, no transcripts - filter all (in gene?)
  df.bQTL_gene_partitioning <- merge(df.bQTL_gene_partitioning, df.bQTL_RNAseq, by = "gene.ID")
  df.bQTL_gene_partitioning <- subset(df.bQTL_gene_partitioning, df.bQTL_gene_partitioning$mode %in% c("up", "down"))
  # df.bQTL_gene_partitioning <- subset(df.bQTL_gene_partitioning, df.bQTL_gene_partitioning$post_gene_5kb == "no")
  
  # ggplot2
  # ASB and RNASEQ ASB density plots - 2 in 1
  hist(df.bQTL_gene_partitioning$POSTfreq, 50)
  hist(df.bQTL_gene_partitioning$POSTfreq, 50)
  
}





# general differential expression between maize and ath
compare_zm_vs_ath_differential_expression <- function(){
  
  # BR regulated 
  path.rnaseq.down_regulated <- "data/expression/BR_repressed_adjpvalue_0.05&l2FC_0.5.txt"
  path.rnaseq.up_regulated <- "data/expression/BR_induced_adjpvalue_0.05&l2FC_0.5.txt"
  
  df.rnaseq.down_regulated <- read.table(path.rnaseq.down_regulated, header = TRUE, sep = "\t", fill = TRUE)
  df.rnaseq.up_regulated <- read.table(path.rnaseq.up_regulated, header = TRUE, sep = "\t", fill = TRUE)
  
  df.bQTL_RNAseq <- df.rnaseq.up_regulated[,c(1,3)]
  df.bQTL_RNAseq <- rbind(df.bQTL_RNAseq, df.rnaseq.down_regulated[,c(1,3)])
  
  
  ###
  length(intersect(test$gene.ID, df.rnaseq.up_regulated$X))
  length(intersect(test$gene.ID, df.rnaseq.down_regulated$X))
  
  
  # ATH regulation 
  df.expression_ath <- read.csv("datasets_paper/arabidopsis_overlap_genelists/Ath_BR_elife_03031_fig3_data1_v2.txt", stringsAsFactors = FALSE, sep = "\t")
  
  df.gene_orthologs
  
  df.gene_conversion.AGPv3_to_AGPv4
  
  
  df.gene_conversion = read.table("F:/junkDNA.ai/HaschSeq_Release/data/PlantGenomeOrthologs/PlantGenomeOrthologs", header = FALSE, stringsAsFactors = F, fill = T)
  
  grepl(pattern="SETIT", x=df.gene_conversion$V2[5])
  
  test = which("ARATH")
  
  file$V2 <- subset(file, file$V2 c("ARATH", "MAIZE"))
  
}


# target gene expression maize and ath comparison 
bQTL_on_ArabidopsisHomolog_geneexpression_evaluation <- function(){
  message("Conserved BRZ1 target genes in Maize and Arabidopsis - Table S4")
  # Figure 1 h, i
  
  
  
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



bQTL_on_brd1_RNASeq <- function(){
  
  
  ## loading new RNASeq (10.14.) 
  df.rnaseq.thomas <- read.csv("/Shared/Everyone/Michael_Thomas/brd1_BL_vs_brd1_Mock_deseq2_8samples.csv", stringsAsFactors = FALSE)
  names(df.rnaseq.thomas)[1] <- "id"
  
  # df.rnaseq <- read.table("/Shared/Everyone/Michael_Thomas/B73_Mo17_rnaseq/B73_Mo17_root_exp_TabS5.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
  
  # df.rnaseq.baldauf <- read.table("/Shared/Everyone/Michael_Thomas/B73_Mo17_rnaseq/Baldauf_2016_dataset.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
  
}




bQTL_on_RNASeq <- function(l.bQTL_gene_partitioning, df.gene_conversion.AGPv3_to_AGPv4){
  
  path.rnaseq.B73vsMo17 <- "data/expression/tpj13414-sup-0003-DataSetS1.txt"
  
  ## loading new RNASeq (10.14.) 
  message("loading rnaseq datasets...")
  strt<-Sys.time()
  df.rnaseq.B73vsMo17 <- read.table(path.rnaseq.B73vsMo17, header = TRUE, sep = "\t", fill = TRUE)
  
  for(i in 1:nrow(df.rnaseq.B73vsMo17)){
    idx.AGPv4 <- which(df.gene_conversion.AGPv3_to_AGPv4$gene.ID.AGPv3 == df.rnaseq.B73vsMo17$Gene_ID_AGPv3[i])
    if(length(idx.AGPv4) > 0){
      df.rnaseq.B73vsMo17$Gene_ID_AGPv4[i] <-  df.gene_conversion.AGPv3_to_AGPv4$gene.ID.AGPv4[idx.AGPv4[1]]
    }
  }
  print(Sys.time() - strt)
  
  
  df.bQTL_gene_partitioning <- l.bQTL_gene_partitioning[[1]]
  df.ASB_w_RNASeq <-  merge(df.bQTL_gene_partitioning, df.rnaseq.B73vsMo17, by.x = "gene.ID", by.y = "Gene_ID_AGPv4", all = FALSE)
  
  #### 
  
  df.rnaseq.B73vsMo17.bg <- unique(df.rnaseq.B73vsMo17[,c(2,20)])
  
  df.ASB_w_RNASeq <- unique(df.ASB_w_RNASeq[,c(1,46)])
  
  hist(log(df.ASB_w_RNASeq$B73_vs_Mo17), 50)
  hist(log(df.rnaseq.B73vsMo17.bg$B73_vs_Mo17), 50)
  
  df.up <- subset(df.bQTL_gene_partitioning, df.bQTL_gene_partitioning$mode == "up")
  df.down <- subset(df.bQTL_gene_partitioning, df.bQTL_gene_partitioning$mode == "down")
  
  v.gns.up <- unique(df.up$gene.ID)
  v.gns.down <- unique(df.down$gene.ID)
  
  v.gns.total <- v.gns.total[!v.gns.total %in% v.gns.up]
  v.gns.total <- v.gns.total[!v.gns.total %in% v.gns.down]
  
  df.ASB_w_RNASeq.up <- subset(df.ASB_w_RNASeq, df.ASB_w_RNASeq$gene.ID %in% v.gns.up)
  hist(df.ASB_w_RNASeq.up$B73_vs_Mo17)
  
  df.ASB_w_RNASeq.up$B73_vs_Mo17[df.ASB_w_RNASeq.up$B73_vs_Mo17 > 5] <- 5
  plot(df.ASB_w_RNASeq.up$POSTfreq, df.ASB_w_RNASeq.up$B73_vs_Mo17, col ="blue")
  points(df.ASB_w_RNASeq.down$POSTfreq, df.ASB_w_RNASeq.down$B73_vs_Mo17, col ="red")
  
  df.ASB_w_RNASeq.random <- df.ASB_w_RNASeq[sample(nrow(df.ASB_w_RNASeq), 5000),]
  
  points(df.ASB_w_RNASeq.random$POSTfreq, df.ASB_w_RNASeq.random$B73_vs_Mo17, col ="gray")
  
  df.ASB_w_RNASeq.down$B73_vs_Mo17[df.ASB_w_RNASeq.down$B73_vs_Mo17 > 5] <- 5
  df.ASB_w_RNASeq.down <- subset(df.ASB_w_RNASeq, df.ASB_w_RNASeq$gene.ID %in% v.gns.down)
  hist(df.ASB_w_RNASeq.down$B73_vs_Mo17)
  
  
  
  df.converted <- data.frame(AGPv4 = as.character(c(as.character(v.gns.up),as.character(v.gns.down))), AGPv3 = NA)
  df.converted <- data.frame(AGPv4 = as.character(v.gns.total), AGPv3 = NA)
  for(i in 1:nrow(df.converted)){
    idx.AGPv3 <- which(df.gene_conversion.AGPv3_to_AGPv4$gene.ID.AGPv4 == df.converted$AGPv4[i])
    if(length(idx.AGPv3) > 0){
      df.converted$AGPv3[i] <-  df.gene_conversion.AGPv3_to_AGPv4$gene.ID.AGPv3[idx.AGPv3[1]]
    }
  }
  
  
  write.csv(df.converted, paste("output/v.B73_targets_bound_non_regulated.",timeStamp, ".csv", sep = "")) # total
  
  
  write.csv(df.converted, paste("output/v.B73_targets_bound_regulated.",timeStamp, ".csv", sep = "")) # 
  
  
  df.bQTL_RNAseq.up <- subset(df.bQTL_RNAseq, df.bQTL_RNAseq$mode == "up")
  v.gns.up <- unique(df.bQTL_RNAseq.up$gene.ID)[!unique(df.bQTL_RNAseq.up$gene.ID) %in% v.gns.up]
  
  df.bQTL_RNAseq.down <- subset(df.bQTL_RNAseq, df.bQTL_RNAseq$mode == "down")
  v.gns.down <- unique(df.bQTL_RNAseq.down$gene.ID)[!unique(df.bQTL_RNAseq.down$gene.ID) %in% v.gns.down]
  
  write.csv(df.converted, paste("output/v.B73_targets_non_bound_non_regulated.",timeStamp, ".csv", sep = ""))
  
  # only annotate - exclusive 
  df.gene_function_partitioning <- merge(df.bQTL_gene_partitioning, df.gene_function, by = "gene.ID", all.x = TRUE)
  # df.geneID_conversion
  
  write.csv(df.gene_function_partitioning, paste("output/df.gene_function_partitioning_with_RNASeq.",timeStamp, ".csv", sep = ""))
  # write.csv(df.bQTL_gene_partitioning, paste("/Shared/Everyone/Michael_Thomas/results/df.bQTL_gene_partitioning_with_RNAseq.",timeStamp, ".csv", sep = ""))
  message("bQTL with rnaseq in gene partitioning " , nrow(df.bQTL_gene_partitioning))
  
  
  message("analyzing rnaseq gene expression profile comparison - prediction and validation")
  
  # post frequency zuordnen - werte groeser oder kleiner depending on the species
  #df.rnaseq.down_regulated <- read.csv(path.rnaseq.down_regulated)
  #df.rnaseq.up_regulated <- read.csv(path.rnaseq.up_regulated)
  
  message("rnaseq upregulated: ", nrow(df.rnaseq.up_regulated), ", downregulated: ", nrow(df.rnaseq.down_regulated), " genes")
  
  # indices in the rnaseq expression data - first orignal - second mutant
  v.indices.colExpressionLevel = c(10,11)
  
  v.modes.diffExpression = c("up", "down")
  
  l.rnaseq.diffExpProfiles = vector(mode = "list", length = 2)
  l.rnaseq.diffExpProfiles[[1]] <- df.rnaseq.up_regulated
  l.rnaseq.diffExpProfiles[[2]] <- df.rnaseq.down_regulated
  
  #df.rnaseq.up_regulated$B73.Shoot
  #df.rnaseq.up_regulated$Mo17.Shoot
  
  # if post frequency > 0.5 (Mo17) and mode up
  df.bQTL_gene_partitioning.set <- subset(df.bQTL_gene_partitioning, df.bQTL_gene_partitioning$mode != "")
  
  l.bQTL_gene_partitioning.ecotypes <- vector(mode = "list", length = 2)
  
  l.bQTL_gene_partitioning.ecotypes[[1]] <- subset(df.bQTL_gene_partitioning, df.bQTL_gene_partitioning$POSTfreq > 0.5)
  l.bQTL_gene_partitioning.ecotypes[[2]] <- subset(df.bQTL_gene_partitioning, df.bQTL_gene_partitioning$POSTfreq < 0.5)
  
  l.bQTL_gene_partitioning.ecotypes[[1]] <- subset(l.bQTL_gene_partitioning.ecotypes[[1]], l.bQTL_gene_partitioning.ecotypes[[1]]$mode != "")
  l.bQTL_gene_partitioning.ecotypes[[2]] <- subset(l.bQTL_gene_partitioning.ecotypes[[2]], l.bQTL_gene_partitioning.ecotypes[[2]]$mode != "")
  
  m.predictionValidation <- matrix(0, nrow = 2, ncol = 2, dimnames =  list(c("original", "mutant"),c("correct", "incorrect"))) 
  
  for(i in 1:2){
    
    for(j in 1:2){
      
      df.bQTL_gene_partitioning.ecotypes.j <- subset(l.bQTL_gene_partitioning.ecotypes[[i]], l.bQTL_gene_partitioning.ecotypes[[i]]$mode == v.modes.diffExpression[j])
      v.gns.ecotype.mode.diffExp <- unique(df.bQTL_gene_partitioning.ecotypes.j$gene.ID)
      
      # real differencial expression
      df.rnaseq.diffExpProfiles.j <- subset(l.rnaseq.diffExpProfiles[[j]], l.rnaseq.diffExpProfiles[[j]]$gene_name %in% v.gns.ecotype.mode.diffExp)
      
      if(v.modes.diffExpression[j] == "up"){
        
        for(k in 1:nrow(df.rnaseq.diffExpProfiles.j)){
          
          # predictions
          df.bQTL_gene_partitioning.ecotypes.jk <- subset(df.bQTL_gene_partitioning.ecotypes.j, df.bQTL_gene_partitioning.ecotypes.j$gene.ID == df.rnaseq.diffExpProfiles.j$gene_name[k])
          # real expression level based target value
          if(df.rnaseq.diffExpProfiles.j[k,v.indices.colExpressionLevel[1]] > df.rnaseq.diffExpProfiles.j[k,v.indices.colExpressionLevel[2]]){
            if(i == 1 && any(df.bQTL_gene_partitioning.ecotypes.jk$mode == "up")){
              m.predictionValidation[i,1] <- m.predictionValidation[i,1] + 1
            }else{
              m.predictionValidation[i,2] <- m.predictionValidation[i,2] + 1
            }
          }else{
            if(i == 2 && any(df.bQTL_gene_partitioning.ecotypes.jk$mode == "up")){
              m.predictionValidation[i,1] <- m.predictionValidation[i,1] + 1
            }else{
              m.predictionValidation[i,2] <- m.predictionValidation[i,2] + 1 
            }
            
          }
          
        }
        
      }else{
        
        for(k in 1:nrow(df.rnaseq.diffExpProfiles.j)){
          df.bQTL_gene_partitioning.ecotypes.jk <- subset(df.bQTL_gene_partitioning.ecotypes.j, df.bQTL_gene_partitioning.ecotypes.j$gene.ID == df.rnaseq.diffExpProfiles.j$gene_name[k])
          # real expression level based target value
          if(df.rnaseq.diffExpProfiles.j[k,v.indices.colExpressionLevel[1]] < df.rnaseq.diffExpProfiles.j[k,v.indices.colExpressionLevel[2]]){
            # predictions
            if(i == 1 & any(df.bQTL_gene_partitioning.ecotypes.jk$mode == "down")){
              m.predictionValidation[i,1] <- m.predictionValidation[i,1] + 1
            }else{
              m.predictionValidation[i,2] <- m.predictionValidation[i,2] + 1 
            }
          }else{
            if(i == 2 & any(df.bQTL_gene_partitioning.ecotypes.jk$mode == "down")){
              m.predictionValidation[i,1] <- m.predictionValidation[i,1] + 1
            }else{
              m.predictionValidation[i,2] <- m.predictionValidation[i,2] + 1 
            } 
          }
        }
      }
    }
  }
  
  
  
}


