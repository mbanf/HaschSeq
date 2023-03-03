
message("loading gene annotation datasets...")
strt<-Sys.time()

path.geneID_conversion <- paste(folder_data,"GeneAnnotation/maize.v3TOv4.geneIDhistory.txt", sep ="/")
df.geneID_conversion <- read.table(path.geneID_conversion, header = FALSE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE, quote = "")[,1:2]
names(df.geneID_conversion) <- c("gene.ID.AGPv3", "gene.ID.AGPv4")

df.gene_conversion.AGPv3_to_AGPv4 <- df.geneID_conversion

path.gff <- paste(folder_data,"revision/B73_Mo17_CSHL.gff", sep ="/")

df.gff <- read.table(path.gff, header = FALSE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE, quote = "")
v.genePartitions <- c("gene", "five_prime_UTR", "CDS", "three_prime_UTR", "exon")
df.gff <- subset(df.gff, df.gff$V3 %in% v.genePartitions)

ids.list <- sapply(df.gff$V9, function(m) {strsplit(m, ";")})

df.gff["id"] <- sapply(df.gff$V9, function(m) {gsub("ID=", "", unlist(strsplit(m, ";"))[1])})
df.gff$id <- gsub("gene:","", df.gff$id)
df.gff$id <- gsub("Parent=","", df.gff$id)
df.gff$id <- gsub("transcript:","", df.gff$id) 
df.gff$id <- gsub("CDS:","", df.gff$id)
df.gff$id <- gsub("\\_.*","", df.gff$id)

# df.gff["transcript.id"] <- sapply(df.gff$V9, function(m) {gsub("ID=", "", unlist(strsplit(m, ";"))[1])})
# df.gff$transcript.id <- gsub("gene:","", df.gff$transcript.id)
# df.gff$transcript.id <- gsub("Parent=transcript:","", df.gff$transcript.id)
# df.gff$transcript.id <- gsub("CDS:","", df.gff$id)
# df.gff$transcript.id <- gsub("\\_.*","", df.gff$id)

df.gene_annotation <- df.gff
names(df.gene_annotation) <- c("chr", "source", "partition", "pos.start", "pos.stop", "aux1", "strand", "aux2", "gene_meta", "gene.ID")


print(Sys.time() - strt)


message("loading GWAS datasets...")
strt<-Sys.time()
path.phenotype_gwas.original <- paste(folder_data,"GWAS/AGPv4_uplifted_Wallace_etal_2014_PLoSGenet_GWAS_hits-150112.txt", sep ="/")
df.phenotype_gwas.original <- read.table(path.phenotype_gwas.original, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
df.phenotype_gwas.original <- subset(df.phenotype_gwas.original, df.phenotype_gwas.original$rmip >= 5)
df.phenotype_gwas <- df.phenotype_gwas.original#
v.traits  <- unique(df.phenotype_gwas$trait)
print(Sys.time() - strt)



message("loading Chip binding peaks") # ZmBZR1 Chip-seq binding peaks
l.path.ChipSeqB73 <- vector(mode = "list", length = 3)
l.path.ChipSeqB73[[1]] =paste(folder_data,"ChipSeq/B73_ChIP/AGTCAA/q2_5fold/q2_5fold.GEM_events.txt", sep ="/")
l.path.ChipSeqB73[[2]] =paste(folder_data,"ChipSeq/B73_ChIP/ATGTCA/q2_5fold/q2_5fold.GEM_events.txt", sep ="/")
l.path.ChipSeqB73[[3]] =paste(folder_data,"ChipSeq/B73_ChIP/CCGTCC/q2_5fold/q2_5fold.GEM_events.txt", sep ="/")



message("load gene orthologs")
# path.gene_orthologs <- paste(folder_data,"GeneAnnotation/Phytozome13_anno_1_1.txt", sep ="/")
# df.gene_orthologs = read.table(path.gene_orthologs, header=T, stringsAsFactors = FALSE, sep = "\t", fill = T) # TODO: when to use this vs the old ones? is this separated?



path.gene_orthologs <- paste(folder_data,"GeneAnnotation/Zmays_284_Ensembl-18_2010-01-MaizeSequence.annotation_info.txt", sep ="/")
df.gene_orthologs = read.table(path.gene_orthologs, header=T, stringsAsFactors = FALSE, sep = "\t", fill = T)



# FILES

path.rnaseq.down_regulated <- "data/expression/BR_repressed_adjpvalue_0.05&l2FC_0.5.txt"
path.rnaseq.up_regulated <- "data/expression/BR_induced_adjpvalue_0.05&l2FC_0.5.txt"
path.rnaseq.ATBES1_targets <- "data/expression/AtBES1_targets.txt"


# putative B73 enhancer (Fig. 4a) regions were extracted from Oka et al. 2017
df.enhancer_H3K9 <- read.table("data/Enhancer_HM/GSE94251_H3K9ac_ist.bedGraph", header = FALSE, sep ="\t", quote = "", stringsAsFactors = FALSE)

df.enhancer_candidates <- read.table("data/Enhancer_HM/HUSK_enhancer_candidates.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
df.enhancer_candidates <- rbind(df.enhancer_candidates, read.table("data/Enhancer_HM/V2_IST_enhancer_candidates.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE))

# 1495 intergenic enhancer regions identified in B73
df.enhancer_genes <- read.table("data/Enhancer_HM/Enhances_constutitive_expressed_specific_genes.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
df.enhancer_genes <- rbind(df.enhancer_genes, read.table("data/Enhancer_HM/Enhances_husk_specific_genes.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE))
df.enhancer_genes <- rbind(df.enhancer_genes, read.table("data/Enhancer_HM/Enhances_V2_IST_specific_genes.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE))

