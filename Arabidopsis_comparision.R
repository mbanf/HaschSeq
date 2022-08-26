

# Figure 1 F
AT_ZM_overlap = function(){
  
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










## feed into zm 
length(v.gns.overlap)



####

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



###


## GO enrichment tables for all areas
# df.GO.annot <- get_GOSlim_annotations(filename= "Athaliana_167_goslim.annot.txt")
# #saveRDS(df.GO.annot, "df.GO.annot.Ath.rds")
# df.GO.annot <- readRDS("df.GO.annot.Ath.rds")
# 
# #df.GO.annot <- readRDS("df.GO_Maize_Trimmed.rds")
# df.GO.annot$acc. <- gsub("\\..*","", df.GO.annot$acc.)
# df.GO.annot <- unique(df.GO.annot)
# 
# 
# lst.df.GO.annot.set <- vector(mode = "list", length = 3)
# names(lst.df.GO.annot.set) <- c("BP", "MF", "CC")
# for(i in 1:length(lst.df.GO.annot.set)){
#   lst.df.GO.annot.set[[i]] <- vector(mode = "list", length = 7)
#   names(lst.df.GO.annot.set[[i]]) <- names(lst.sets)
#   for(s in 1:length(lst.sets)){
#     df.GO.annot.set <- subset(df.GO.annot, df.GO.annot$acc. %in% lst.sets[[s]])
#     df.GO.annot.set <- unique(df.GO.annot.set)
#     df.GO.annot.set <- subset(df.GO.annot.set, df.GO.annot.set$Ontology == names(lst.df.GO.annot.set)[[i]])    
#     df.tb.annot <- as.data.frame(table(df.GO.annot$Term[df.GO.annot$Term != ""]), stringsAsFactors = FALSE)
#     df.tb.annot.set <- as.data.frame(table(df.GO.annot.set$Term[df.GO.annot.set$Term != ""]), stringsAsFactors = FALSE)
#     df.tb.annot.set <- merge(df.tb.annot.set, df.tb.annot, by = "Var1")
#     df.enrichment <- data.frame(GO = character(nrow(df.tb.annot.set)), p.val = numeric(nrow(df.tb.annot.set)),
#                                 percent = numeric(nrow(df.tb.annot.set)), genes = numeric(nrow(df.tb.annot.set)), stringsAsFactors = FALSE)
#     n.genes.sset <- length(unique(df.GO.annot.set$acc.))
#     n.genes <- length(unique(df.GO.annot$acc.))
#     for(j in 1:nrow(df.tb.annot.set)){  
#       n.inset <- df.tb.annot.set$Freq.x[j]
#       n.genomewide <- df.tb.annot.set$Freq.y[j]
#       #n.genomewide <- df.tb.annot$Freq[df.tb.annot$Var1 == df.tb.annot$Var1[j]]
#       ### global enrichment test - gene basis
#       hitInSample <- n.inset
#       sampleSize <- n.genes.sset
#       hitInPop <- n.genomewide #sum(tb.rate_limiting_domains$Freq)
#       failInPop <- n.genes - hitInPop #(nrow(df.global.domains) - hitInPop)
#       p.val <- print(phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE))
#       #fisher.test(matrix(c(hitInSample-1, hitInPop, failInPop, sampleSize), 2, 2), alternative='less');
#       foldChange <- (hitInSample / sampleSize) / (hitInPop / failInPop)
#       #   mat.count <- matrix(c(n.inset,n.genes.sset - n.inset,  n.genomewide ,n.genes - n.genomewide), ncol = 2, byrow = FALSE)
#       #   p.val <- fisher.test(mat.count)$p.value
#       #   
#       df.enrichment$GO[j] <- df.tb.annot.set$Var1[j]
#       df.enrichment$p.val[j] <- p.val
#       df.enrichment$percent[j] <- (n.inset / n.genes.sset)
#       df.enrichment$genes[j] = n.inset
#     }
#     df.enrichment <- df.enrichment[order(df.enrichment$p.val),]
#     #df.enrichment <- df.enrichment[order(df.enrichment$percent),]
#     lst.df.GO.annot.set[[i]][[s]] <- df.enrichment
#   }
# }
# 
# 
# v.sets <- c("BiologicalProcess.xlsx", "MolecularFunction.xlsx", "CellularCompartment.xlsx")
# # write genes 
# for(i in 1:length(v.sets)){
#   a123_466 <- lst.df.GO.annot.set[[i]][[1]]
#   a12_1038 <- lst.df.GO.annot.set[[i]][[2]]
#   a23_525 <- lst.df.GO.annot.set[[i]][[3]]
#   a13_420 <- lst.df.GO.annot.set[[i]][[4]]
#   a1_1487 <- lst.df.GO.annot.set[[i]][[5]]
#   a2_2272 <- lst.df.GO.annot.set[[i]][[6]]
#   a3_2784 <- lst.df.GO.annot.set[[i]][[7]]
#   save.xlsx(v.sets[i], a123_466, a12_1038, a23_525, a13_420, a1_1487, a2_2272, a3_2784)
#   
# }



