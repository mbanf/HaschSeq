rm(list=ls()) # clear workspace 


folder_output = "output_test" 
folder_tmp = "tmp_revision"
folder_data = "data"

b.write_paper = FALSE

source("config.R")
source("utils.R")

## Extract condition
require(reshape2)
require(ggplot2)
library(dplyr)
require(seqinr)
require(Biostrings)
require(VariantAnnotation)
require(foreach)
require(doParallel)
require(ggplot2)
require(BSgenome)
require(VariantAnnotation)
require(plotly)
require(rtracklayer)
library(stringr)


write_bed <- function(df.QTLs){
  
  df.QTLs.bed <- df.QTLs[,c("B73-chr", "B73-pos", "B73-pos")]
  colnames(df.QTLs.bed) <- c('chrom', 'chromStart', 'chromEnd')
  
  for(i in 1:nrow(df.QTLs)){
    if(df.QTLs$POSTfreq[i] > 0.5){
      df.QTLs.bed$chrom[i] = df.QTLs$`B73-chr`[i]
      df.QTLs.bed$chromStart[i] = df.QTLs$`B73-pos`[i] - 1
      df.QTLs.bed$chromEnd[i] = df.QTLs$`B73-pos`[i]
      
    }else if(df.QTLs$POSTfreq[i] < 0.5){
      df.QTLs.bed$chrom[i] = df.QTLs$`Mo17-chr`[i]
      df.QTLs.bed$chromStart[i] = df.QTLs$`Mo17-pos`[i] - 1
      df.QTLs.bed$chromEnd[i] = df.QTLs$`Mo17-pos`[i]
    }
  }
  
  df.QTLs.bed
}

extract_sequences <- function(df.bQTLs,
                              genome_sequences,
                              s.env_bps = 100,
                              n.chromosomes = 10
                              ){
  
  df.bQTLs["id"] <- seq(1, nrow(df.bQTLs))
  si <- seqinfo(genome_sequences)
  message("extract sequences surrounding ASBs from diploid genome")
  
  v.sequences <- character(length = nrow(df.bQTLs))
  
  for(i in 1:n.chromosomes){ 
    print(i)
    
    chr_name <- paste("B73-chr", i, sep = "")
    chr_seq <- getSeq(genome_sequences, GRanges(chr_name, IRanges(1, seqlengths(si)[chr_name])))
    
    df.bQTLs.B73.i <- subset(df.bQTLs, df.bQTLs$`B73-chr` == chr_name)
    df.bQTLs.B73.i <- subset(df.bQTLs.B73.i, df.bQTLs.B73.i$POSTfreq > 0.5)
    
    for(j in 1:nrow(df.bQTLs.B73.i)){
      id <- df.bQTLs.B73.i$id[j]
      v.sequences[id] <- subseq(chr_seq, start=df.bQTLs.B73.i$`B73-pos`[j] - s.env_bps, end=df.bQTLs.B73.i$`B73-pos`[j] + s.env_bps)  
    }
    
    chr_name <- paste("Mo17-chr", i, sep = "")
    chr_seq <- getSeq(genome_sequences, GRanges(chr_name, IRanges(1, seqlengths(si)[chr_name])))
    
    df.bQTLs.Mo17.i <- subset(df.bQTLs, df.bQTLs$`Mo17-chr` == chr_name)
    df.bQTLs.Mo17.i <- subset(df.bQTLs.Mo17.i, df.bQTLs.Mo17.i$POSTfreq < 0.5)
    
    for(j in 1:nrow(df.bQTLs.Mo17.i)){
      id <- df.bQTLs.Mo17.i$id[j]
      v.sequences[id] <- subseq(chr_seq, start=df.bQTLs.Mo17.i$`Mo17-pos`[j] - s.env_bps, end=df.bQTLs.Mo17.i$`Mo17-pos`[j] + s.env_bps)  
    }
  }
  
  descs <-  paste("asb", df.bQTLs$id, sep = "")
  asb_environment <- DNAStringSet(v.sequences)
  names(asb_environment) <- descs
  
  asb_environment
}


load_gene_annotation <- function(path.gene_annotation){
  
  df.gff <- read.table(path.gene_annotation, header = FALSE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE, quote = "")
  v.genePartitions <- c("gene", "five_prime_UTR", "CDS", "three_prime_UTR", "exon")
  df.gff <- subset(df.gff, df.gff$V3 %in% v.genePartitions)
  
  ids.list <- sapply(df.gff$V9, function(m) {strsplit(m, ";")})
  
  df.gff["id"] <- sapply(df.gff$V9, function(m) {gsub("ID=", "", unlist(strsplit(m, ";"))[1])})
  df.gff$id <- gsub("gene:","", df.gff$id)
  df.gff$id <- gsub("Parent=","", df.gff$id)
  df.gff$id <- gsub("transcript:","", df.gff$id) 
  df.gff$id <- gsub("CDS:","", df.gff$id)
  df.gff$id <- gsub("\\_.*","", df.gff$id)

  names(df.gff) <- c("chr", "source", "partition", "pos.start", "pos.stop", "aux1", "strand", "aux2", "gene_meta", "gene.ID")
  
  df.gff
}


bQTL_scatterplot <- function(df.bQTLs, n.chromosomes = 10){

  ### scatter plot of the postfrequencies - should be updated to > , < 0.5
  v.offset <- numeric(n.chromosomes)
  v.start <- numeric(n.chromosomes)
  
  for(i in 1:n.chromosomes){
    chr <- paste("B73-chr", i, sep = "")
    df.bQTLs.i <- subset(df.bQTLs, df.bQTLs$`B73-chr` == chr)
    v.offset[i] <-  max(df.bQTLs.i$`B73-pos`)
    v.start[i] <- sum(v.offset[1:i])
  }
  
  v.start <- c(0, v.start)
  
  # directionality plot 
  df.scatterplot <- subset(df.bQTLs, df.bQTLs$`B73-chr` == "B73-chr1")
  
  for(i in 1:(n.chromosomes - 1)){
    chr <- paste("B73-chr", (i+1), sep = "")
    df.bQTLs.i <- subset(df.bQTLs, df.bQTLs$`B73-chr` == chr)
    df.bQTLs.i$`B73-pos` <- df.bQTLs.i$`B73-pos` + sum(v.offset[1:i])
    df.scatterplot <- rbind(df.scatterplot, df.bQTLs.i)
  }
  
  p6 <- ggplot(df.scatterplot, aes(x = `B73-pos`, y = POSTfreq, fill = POSTfreq)) +
    geom_point(shape = 21, size = 1,  alpha = 0.5, stroke = 0.0) + theme_bw() 
  
  p6 <- p6 + scale_fill_continuous(low = "blue", high = "green2")
  p6 <- p6 + geom_vline(xintercept = v.start, col='black', lwd=0.5, linetype="dashed")
  
  p6
}


# datasets
path.bQTLs <- "data/revision/B73.Mo17.WW.PF.uni_EG80.GT.RN.BZR1.GEM.csv"
path.bindingPeaks = "data/revision/BZR1_6.GEM_events.bed"
path.genomic_sequences = "data/revision/ref_B73Mo17.fasta"
path.gene_annotation <- "data/revision/B73_Mo17_CSHL.gff"
path.phenotypic_snps <- "data/revision/gwas/remapped_AGPv5_chr_uplifted_Wallace_etal_2014_PLoSGenet_GWAS_hits.bed"

path.raw_input_reads <- "data/revision/BZR7_12.RPGCq255.b1.ct.bigwig"
# path.raw_input_reads <- "data/revision/BZR7_12.10.bp.rand.EG.bedgraph"

# The data is generated from BZR1 ChIP in all six replicates in B73xMo17, mapped to a diploid genome, i.e. the two parental genomes pasted together. 
# Only unique mapping reads were kept, which means that only reads with sequence differences between the two inbred lines are mapped, 
# identical regions in the genome have no reads. However, since the ASB analysis only concerns SNPs that does not affect this part of the analysis. 
# The SNPs were generated by whole genome alignment, SNPs that map to several places in one of the parental genomes were excluded.
# That leaves us with 10.973.994 SNPs.
df.snps <- read.table(path.bQTLs,  sep ="\t", stringsAsFactors = FALSE)
names(df.snps) <- c("B73-chr", "B73-pos", "B73-base", "Mo17-base", "Mo17-position", "Genotype",
                      "CPM-B73", "CPM-Mo17", "Peak-B73", "Peak-Mo17", "post-frequency", "Reads-B73", "Reads-Mo17")

message(nrow(df.snps), " SNPs from BZR1 ChIP in all six replicates in B73xMo17")

df.snps["refCount"] <- df.snps["Reads-B73"]
df.snps["altCount"] <- df.snps["Reads-Mo17"]
df.snps["totalCount"] <- df.snps$refCount + df.snps$altCount
df.snps["refAllele"] <- df.snps["B73-base"]
df.snps["altAllele"] <- df.snps["Mo17-base"]


#### blacklisting ###
# TODO: filtering here or after min read allele 
message("Identify systematic bias in the input read data for blacklisting")


minReadsPerAllele = 1 # only work with covered snps
df.snps <- subset(df.snps, df.snps$refCount >= minReadsPerAllele & df.snps$altCount >= minReadsPerAllele)

df.snps["Mo17-chr"] <- as.character(unlist(lapply(strsplit(df.snps$`Mo17-position`, ":"), function(m) m[[1]])))
df.snps["Mo17-pos"] <- as.numeric(unlist(lapply(strsplit(df.snps$`Mo17-position`, ":"), function(m) m[[2]])))
df.snps <- subset(df.snps, df.snps$`B73-chr` %in% paste("B73-chr", 1:10, sep = "") | 
                    df.snps$`Mo17-chr` %in% paste("Mo17-chr", 1:10, sep = ""))
df.snps <- subset(df.snps, df.snps$`B73-chr` %in% paste("B73-chr", 1:10, sep = "")) # hard filter on B73 (no scaffolds)

df.snps["POSTfreq"] <- df.snps$refCount / df.snps$totalCount
df.snps["POSTallele"] <- ifelse(df.snps$altCount > df.snps$refCount, df.snps$altAllele, df.snps$refAllele)


source("utils/blacklisting.R")
df.snps.blacklisting <- blacklisting_background(df.snps, path.raw_input_reads, th.bp_ranges = 500)
saveRDS(df.snps.blacklisting, "df.snps.blacklisting.rds") # put in tmp folder for reusage 

# TODO: only use blacklisted 

# FIGURE TO PLOT SYSTEMATIC BIAS

# Correlation Postfrequency (SNP) vs Avg. Postfrequency (region +/- 500 bp )



#### Select background candidates #### 

minReadsPerAllele = 1
maxReadsPerAllele = 15
maxReadDepth = 20

df.bgSNPs.candidates <- subset(df.snps.blacklisting, df.snps.blacklisting$refCount >= minReadsPerAllele & df.snps.blacklisting$altCount >= minReadsPerAllele)
df.bgSNPs.candidates <- subset(df.bgSNPs.candidates, 
                               df.bgSNPs.candidates$refCount <= maxReadsPerAllele & df.bgSNPs.candidates$altCount <= maxReadsPerAllele)
df.bgSNPs.candidates <- subset(df.bgSNPs.candidates, df.bgSNPs.candidates$totalCount <= maxReadDepth)  
# 
# df.bgSNPs.candidates["Mo17-chr"] <- as.character(unlist(lapply(strsplit(df.bgSNPs.candidates$`Mo17-position`, ":"), function(m) m[[1]])))
# df.bgSNPs.candidates["Mo17-pos"] <- as.numeric(unlist(lapply(strsplit(df.bgSNPs.candidates$`Mo17-position`, ":"), function(m) m[[2]])))
# df.bgSNPs.candidates <- subset(df.bgSNPs.candidates, df.bgSNPs.candidates$`B73-chr` %in% paste("B73-chr", 1:10, sep = "") | 
#                                df.bgSNPs.candidates$`Mo17-chr` %in% paste("Mo17-chr", 1:10, sep = ""))
# df.bgSNPs.candidates <- subset(df.bgSNPs.candidates, df.bgSNPs.candidates$`B73-chr` %in% paste("B73-chr", 1:10, sep = "")) # hard filter on B73 (no scaffolds)
# 
# df.bgSNPs.candidates["POSTfreq"] <- df.bgSNPs.candidates$refCount / df.bgSNPs.candidates$totalCount
# df.bgSNPs.candidates["POSTallele"] <- ifelse(df.bgSNPs.candidates$altCount > df.bgSNPs.candidates$refCount, df.bgSNPs.candidates$altAllele, df.bgSNPs.candidates$refAllele)

df.bgSNPs.candidates <- subset(df.bgSNPs.candidates, df.bgSNPs.candidates$POSTfreq != 0.5) # This ignores bg SNPS with allelic bias == 0.5 !!!!

# message("background SNP candidates after min reads per allele filter: ", nrow(df.bgSNPs.candidates))
message(nrow(df.bgSNPs.candidates), " non-significant SNPs located outside high-confidence BZR1 binding peaks")
hist(df.bgSNPs.candidates$POSTfreq, breaks = 100, main = "Allelic bias of non-significant (background) SNPs located outside high-confidence BZR1 binding peaks", xlab = "allelic bias", cex.main = 0.7)

####

minReadsPerAllele = 1
minReadDepth = 25
df.snps <- subset(df.snps, df.snps$refCount >= minReadsPerAllele & df.snps$altCount >= minReadsPerAllele)
df.snps <- subset(df.snps, df.snps$refCount >= minReadDepth | df.snps$altCount >= minReadDepth)
message("input bQTL after min reads per allele filter: ", nrow(df.snps))
# df.snps <- subset(df.snps, df.snps$totalCount >= minReadDepth)  
df.snps["Mo17-chr"] <- as.character(unlist(lapply(strsplit(df.snps$`Mo17-position`, ":"), function(m) m[[1]])))
df.snps["Mo17-pos"] <- as.numeric(unlist(lapply(strsplit(df.snps$`Mo17-position`, ":"), function(m) m[[2]])))
df.snps <- subset(df.snps, df.snps$`B73-chr` %in% paste("B73-chr", 1:10, sep = "") | 
                    df.snps$`Mo17-chr` %in% paste("Mo17-chr", 1:10, sep = ""))
df.snps <- subset(df.snps, df.snps$`B73-chr` %in% paste("B73-chr", 1:10, sep = "")) # hard filter on B73 (no scaffolds)

message(nrow(df.snps), " SNPs from BZR1 ChIP in all six replicates in B73xMo17 (10 chromosomes only)")



prob.bias <- median(df.snps$refCount / df.snps$totalCount)
message("bias probability (median): ", prob.bias)
th.p.bQTL = 0.001
fdr <- p.adjust(sapply(1:nrow(df.snps), function(i) binom.test(as.integer(df.snps$altCount[i]), as.integer(df.snps$totalCount[i]), p = 1 - prob.bias)$p.value), "bonferroni")
df.snps["p-value (corrected)"] <- fdr
df.snps["POSTfreq"] <- df.snps$refCount / df.snps$totalCount
df.snps["POSTallele"] <- ifelse(df.snps$altCount > df.snps$refCount, df.snps$altAllele, df.snps$refAllele)


# Of the remaining 312656 SNPs, we determined significant variation of median allele
# frequency of 0.497 using a binomial test with a p-value cutoff of ≤ 0.001 adjusted for
# multiple testing using Bonferroni correction (n=46286 SNPs)
df.bQTLs <- subset(df.snps, df.snps$`p-value (corrected)` < th.p.bQTL)
df.bQTLs <- subset(df.bQTLs, df.bQTLs$POSTfreq < 1 & df.bQTLs$POSTfreq > 0)

message(nrow(df.bQTLs), " bQTLs (significance p < ", th.p.bQTL, " )")
hist(df.bQTLs$POSTfreq, breaks = 100, main = "Allelic bias of bQTL", xlab = "allelic bias")





message("Estimate significant bQTLs in binding peaks.")

df.peaks <- read.table(path.bindingPeaks, header = FALSE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE, quote = "")
names(df.peaks) <- c("seqnames", "start", "end", "id", "val", "nan")
df.peaks$seqnames <- str_replace(df.peaks$seqnames, "chr", "")
df.peaks["center"] <- as.numeric(unlist(lapply(strsplit(df.peaks$id, ":"), function(m) m[[2]])))


source("utils/binding_peaks_revision.R")
df.bQTLs <- merged_peaks(df.bQTLs, df.peaks)

message(nrow(df.bQTLs), " SNPs located within high-confidence BZR1 binding peaks (all 6 replicated merged)")
hist(df.bQTLs$POSTfreq, breaks = 100, main = "Allelic bias of bQTL located within high-confidence BZR1 binding peaks", xlab = "allelic bias", cex.main = 0.8)

write.csv(df.bQTLs, "32009_bQTLs_in_peaks.csv", quote = F, row.names = F)
write.table(write_bed(df.bQTLs), "32009_bQTLs_in_peaks.bed", row.names = F, col.names = F, quote = F, sep = "\t")


# TODO: put blacklisting here? ### 


source("utils/linkage_disequilibrium.R")
df.ASBs <- linkage_disequlilibrium(df.bQTLs,
                                    df.peaks,
                                    mode = "peak_distance",
                                    s.disequilibriumDistance = 75,
                                    n.chromosomes = 10)


message(nrow(df.ASBs), " ASBs (bQTLs in peaks after linkage disequilibrium)")
bQTL_scatterplot(df.ASBs)
hist(df.ASBs$POSTfreq, breaks = 100, main = paste("Allelic bias of", nrow(df.ASBs), "ASBs after linkage disequilibrium (min peak distance)"), 
     xlab = "allelic bias", cex.main = 0.7)

write.csv(df.ASBs, "7780_ASBs_75bp_min_peak_distance.csv", quote = F, row.names = F)
write.table(write_bed(df.ASBs), "7780_ASBs_75bp_min_peak_distance.bed", row.names = F, col.names = F, quote = F, sep = "\t")

asb_environment <- extract_sequences(df.ASBs,
                                    genome_sequences = readDNAStringSet(path.genomic_sequences),
                                    s.env_bps = 100,
                                    n.chromosomes = 10)

writeXStringSet(asb_environment, "200_bp_environment_sequences_from_7780_ASBs_75bp_min_peak_distance.fasta")


message("perform gene partitioning")
df.gene_annotation <- load_gene_annotation(path.gene_annotation)
  
# sanity check of gene annotation
#print(table(subset(df.gene_annotation, df.gene_annotation$chr %in% paste("B73-chr", seq(1,10), sep = ""))$partition))
#print(table(subset(df.gene_annotation, df.gene_annotation$chr %in% paste("Mo17-chr", seq(1,10), sep = ""))$partition))
# TODO: add UTR to the mo17 ?

source("utils/genomic_location_revision.R")
df.ASBs_genomic_location = add_genomic_location_bQTLs(df.ASBs,  
                                                      df.gene_annotation,
                                                      v.gene_partitions = c("gene", "five_prime_UTR",  "CDS", "three_prime_UTR", "exon"),
                                                      n.chromosomes = 10,
                                                      n.cpus = 5) 
                              
# saveRDS(df.ASBs_genomic_location, paste(folder_tmp, "7780_ASBs_w_genomic_location.rds", sep = "/"))                  
message(nrow(df.ASBs_genomic_location), " high confidence allele-specific BZR1 binding sites located near ", length(unique(df.ASBs_genomic_location$gene.ID)), " flanking genes")
write.csv(df.ASBs_genomic_location, "7780_ASBs_w_genomic_location.csv", quote = F, row.names = F)

df.distribution <- data.frame(ASBs = apply(df.ASBs_genomic_location[,v.partitions], 2, table)["yes",])
df.distribution["%ASBs"] <-  round(df.distribution$ASBs / sum(df.distribution$ASBs) * 100, 1)

message("Distribution of ASBs positions around / withing genes:")
print(df.distribution)



### BACKGROUND 


message("Control background SNP (bgSNP) sampling")
# Functional GWAS variants have been shown to be significantly enriched in gene proximal
# regions. Therefore, control background SNPs (bgSNPs) were proportionally sampled (excluding ASBs) 
# per chromosome and genomic location (i.e. 5 - 1kb upstream, 1kb upstream - TSS, 5’UTR, exon, intron, 3’UTR, TTS - 1 kb downstream, intergenic) to match
# the genomic distribution of the ASBs dataset. Additionally, we checked that ASBs and bgSNPs showed a similar minor allele-frequency . In total,
# 168950 bgSNPs were sampled, yielding approximately 50 times as many background SNPs per genome location, compared to the number of ASBs within each location
df.snps_genomic_location = add_genomic_location_bQTLs(df.bgSNPs.candidates,  
                                                      df.gene_annotation,
                                                      v.gene_partitions = c("gene", "five_prime_UTR",  "CDS", "three_prime_UTR", "exon"),
                                                      n.chromosomes = 10,
                                                      n.cpus = 2) 
# saveRDS(df.snps_genomic_location, paste(folder_tmp, "6205508_bgSNPs_candidates_w_genomic_location.rds", sep = "/"))


source("utils/qtl_background_sampling.R")

v.partitions <- c("promoter_5kb", "promoter_1kb", "gene", "five_prime_UTR", "exon", "intron", "three_prime_UTR", "post_gene_1kb", "non_genic") 

res <- create_background_QTLs(df.ASBs_genomic_location,
                              df.snps_genomic_location, 
                              multiplyer = 50, 
                              v.partitions = v.partitions, 
                              n.chromosomes = 10,
                              seed = 1234)
                                    
df.bgSNPs <- res$df.bgSNPs
message(nrow(df.bgSNPs), " sampled, non-significant (background) SNPs located outside high-confidence BZR1 binding peaks (all 6 replicated merged)")
hist(df.bgSNPs$POSTfreq, breaks = 100, main = paste("Allelic bias of", nrow(df.bgSNPs), "sampled, non-significant (background) SNPs located outside high-confidence BZR1 binding peaks"), 
     xlab = "allelic bias", cex.main = 0.7)


bQTL_scatterplot(df.bgSNPs)
# distribution - scatter
saveRDS(df.bgSNPs, "400435_bgSNPs_w_genomic_location.rds")          
write.csv(df.bgSNPs, "400435_background_SNPs_w_genomic_location.csv", quote = F, row.names = F)
write.table(write_bed(df.bgSNPs), "400435_background_SNPs.bed", row.names = F, col.names = F, quote = F, sep = "\t")



### GWAS ###

message("-------------- GWAS ------------- ")


# Fig. 4. ASBs of ZmBZR1 are linked to growth and disease related traits. 
# b) Association of ASBs with the curated 4041 significant GWAS hits for selected phenotypes of
# the NAM population. Abbreviations: Av: Average; intern: internode; len: length; w. pl.: whole plant;
# b.: below c) VCAP Variance component analysis results. Variance explained (h2) by the ASB-
#   SNP set (bars) and background SNP set (violin plots, derived from the permutation results).
# 29Orange color in the bars denotes a significantly higher variance explained (h2) by ASBs than
# expected by chance (one-sided permutation test < 0.1).

df.phenotypic_snps <- read.table(path.phenotypic_snps)
names(df.phenotypic_snps) <- c("chr", "pos.start", "pos.end", "trait")
v.traits  <- unique(df.phenotypic_snps$trait)
print(table(df.phenotypic_snps$trait))


source("utils/gwas_enrichment_revision.R")

perform_gwas <- function(df.ASBs,
                         df.bgSNPs, 
                         B73.only = F){
  
  if(B73.only){
    df.ASBs <- subset(df.ASBs, df.ASBs$POSTfreq > 0.5)
    df.bgSNPs <- subset(df.bgSNPs, df.bgSNPs$POSTfreq > 0.5)
  }
  
  res.ASBs <- gwas_enrichment(df.ASBs,
                             df.phenotypic_snps, 
                             v.partitions,
                             n.chromosomes = 10,
                             s.dist_to_phenotypic_snp=2000)
  
  res.bgSNPs <- gwas_enrichment(df.bgSNPs, 
                               df.phenotypic_snps, 
                               v.partitions,
                               n.chromosomes = 10,
                               s.dist_to_phenotypic_snp=2000)
  
  d.gsea <- context_gsea(tb.context.sample=res.ASBs$l.number_per_trait, 
                         tb.context.population=res.bgSNPs$l.number_per_trait, 
                         sampleSize = nrow(df.ASBs),
                         popSize = nrow(df.bgSNPs),
                         col_grad = c("red", "yellow"),
                         scale_factor = 10)
  
  d.gsea
}

d.gsea <- perform_gwas(df.ASBs, df.bgSNPs, B73.only = F)
write.csv(d.gsea, "2000bp_dist_phenotypic_gwas.csv",  row.names = FALSE, quote = FALSE)

d.gsea <- perform_gwas(df.ASBs, df.bgSNPs, B73.only = T)
write.csv(d.gsea, "2000bp_dist_phenotypic_gwas_B73_postfreq_only.csv",  row.names = FALSE, quote = FALSE)



message("Genomic feature profiling of ASBs (Methylation, Motifs and DNAseI, Enhancers) ")

df.ASBs <- readRDS(paste(folder_tmp, "7780_ASBs_w_genomic_location.rds", sep = "/"))
df.bgSNPs <- readRDS(paste(folder_tmp, "400435_bgSNPs_w_genomic_location.rds", sep = "/"))



message("------------ Methylation ---------")
# Methylation levels for CG, CHG, and CHH for B73 and Mo17 were extracted from Regulski et al. 2013
# Methylation frequency versus distance (up to +/- 2 kbp) around
# each ASB were averaged over 20bp bins, and visualized by regions bound by BZR1 with
# either high or low affinity levels depending on the inbred line. For B73, high and low
# affinity bound regions were defined by a post frequency of >= 0.85 or <= 0.15, respectively
# and oppositely for Mo17 by a post frequency <= 0.15 and >= 0.85), respectively.

source("utils/methylation_revision.R")
v.parents = c("B73", "Mo17", "B73", "Mo17", "B73", "Mo17")
v.variants = c("CG", "CG", "CHH", "CHH", "CHG", "CHG")
th.bp_offset = 20

# PERFORM METHYLATION OF ALL 
methylation_occupancy(df.ASBs, 
                     df.bgSNPs,
                     v.filenames = c("data/revision/methylation/B73_methylation_CG.bw", 
                                     "data/revision/methylation/Mo17_methylation_CG.bw",
                                     "data/revision/methylation/B73_methylation_CHH.bw",
                                     "data/revision/methylation/Mo17_methylation_CHH.bw",
                                     "data/revision/methylation/B73_methylation_CHG.bw",
                                     "data/revision/methylation/Mo17_methylation_CHG.bw"),
                     v.formats = c("bigWig", "bigWig", "bigWig", "bigWig", "bigWig", "bigWig"),
                     v.parents = v.parents,
                     v.variants = v.variants,
                     v.has_ranges = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
                     th.bp_offset = th.bp_offset,
                     n.chromosomes = 10,
                     n.cpus = 1,
                     folder = folder_tmp)


res <- load_methylation(df.ASBs, 
                        df.bgSNPs,
                        v.parents,
                        v.variants,
                        folder = folder_tmp)

df.ASBs <- res$df.ASBs
df.bgSNPs <- res$df.bgSNPs

dev.off()
k = 2
title = paste("allelic biases vs methylation levels (averaged per ASB using ",  th.bp_offset, "bp window)", sep = "")
par(mfrow=c(3,2))
scatter_plot(df.ASBs$POSTfreq, df.ASBs$`CG-B73`, xlab = "allelic bias", ylab = "methylation level (B73 CpG)")
scatter_plot(df.ASBs$POSTfreq, df.ASBs$`CG-Mo17`, xlab = "allelic bias", ylab = "methylation level (MO17 CpG)")
scatter_plot(df.ASBs$POSTfreq, df.ASBs$`CHG-B73`, xlab = "allelic bias", ylab = "methylation level (B73 CHG)")
scatter_plot(df.ASBs$POSTfreq, df.ASBs$`CHG-Mo17`, xlab = "allelic bias", ylab = "methylation level (MO17 CHG)")
scatter_plot(df.ASBs$POSTfreq, df.ASBs$`CHH-B73`, xlab = "allelic bias", ylab = "methylation level (B73 CHH)")
scatter_plot(df.ASBs$POSTfreq, df.ASBs$`CHH-Mo17`, xlab = "allelic bias", ylab = "methylation level (MO17 CHH)")
mtext(title, side = 3, line = -2, outer = TRUE, cex = 0.8)



# repeat with displaying only methylation explained ASBs
dev.off()
title = paste("allelic biases vs methylation levels (averaged per ASB using ",  th.bp_offset, "bp window)", sep = "")
par(mfrow=c(2,2))
idx = which(df.ASBs$`CG-Mo17` <= 0.3 & df.ASBs$`CG-B73` >= 0.7)
scatter_plot(df.ASBs$POSTfreq[idx], df.ASBs$`CG-B73`[idx], xlab = "allelic bias", ylab = "methylation level (B73 CpG)")

idx = which(df.ASBs$`CG-Mo17` >= 0.7 & df.ASBs$`CG-B73` <= 0.3)
scatter_plot(df.ASBs$POSTfreq[idx], df.ASBs$`CG-Mo17`[idx], xlab = "allelic bias", ylab = "methylation level (MO17 CpG)")

idx = which(df.ASBs$`CHG-Mo17` <= 0.3 & df.ASBs$`CHG-B73` >= 0.7)
scatter_plot(df.ASBs$POSTfreq[idx], df.ASBs$`CHG-B73`[idx], xlab = "allelic bias", ylab = "methylation level (B73 CHG)")

idx = which(df.ASBs$`CHG-Mo17` >= 0.7 & df.ASBs$`CHG-B73` <= 0.3)
scatter_plot(df.ASBs$POSTfreq[idx], df.ASBs$`CHG-Mo17`[idx], xlab = "allelic bias", ylab = "methylation level (MO17 CHG)")

mtext(title, side = 3, line = -2, outer = TRUE, cex = 0.8)


### 

dev.off()
k = 2
title = paste("allelic biases vs differential methylation levels (averaged per ASB using ",  th.bp_offset, "bp window)", sep = "")
par(mfrow=c(3,1))

idx1 = (df.ASBs$`CG-Mo17` <= 0.1 & df.ASBs$`CG-B73` >= 0.7)
idx2 = (df.ASBs$`CG-Mo17` >= 0.7 & df.ASBs$`CG-B73` <= 0.1)
idx = which(idx1 | idx2)
col_vals = rep("black", nrow(df.ASBs)) 
col_vals[idx] = "red"

scatter_plot(df.ASBs$POSTfreq, df.ASBs$`CG-Mo17` - df.ASBs$`CG-B73`, 
             xlab = "allelic bias", ylab = "methylation level difference (CpG)", 
             x_limit = c(0,1), y_limit = c(-1,1), 
             col_vals = col_vals)

### 
idx1 = (df.ASBs$`CHG-Mo17` <= 0.1 & df.ASBs$`CHG-B73` >= 0.7)
idx2 = (df.ASBs$`CHG-Mo17` >= 0.7 & df.ASBs$`CHG-B73` <= 0.1)
idx = which(idx1 | idx2)
col_vals = rep("black", nrow(df.ASBs)) 
col_vals[idx] = "red"

scatter_plot(df.ASBs$POSTfreq, 
             df.ASBs$`CHG-Mo17` - df.ASBs$`CHG-B73`,
             xlab = "allelic bias", ylab = "methylation level difference (CHG)", 
             x_limit = c(0,1), y_limit = c(-1,1), 
             col_vals = col_vals)



### 
idx1 = (df.ASBs$`CHH-Mo17` <= 0.1 & df.ASBs$`CHH-B73` >= 0.7)
idx2 = (df.ASBs$`CHH-Mo17` >= 0.7 & df.ASBs$`CHH-B73` <= 0.1)
idx = which(idx1 | idx2)
col_vals = rep("black", nrow(df.ASBs)) 
col_vals[idx] = "red"

scatter_plot(df.ASBs$POSTfreq, 
             df.ASBs$`CHH-Mo17` - df.ASBs$`CHH-B73`,
             xlab = "allelic bias", ylab = "methylation level difference (CHH)", 
             x_limit = c(0,1), y_limit = c(-1,1), 
             col_vals = col_vals)

mtext(title, side = 3, line = -2, outer = TRUE, cex = 0.8)



####

methylation_occupancy_distance_to_asb(df.ASBs, 
                                      df.bgSNPs,
                                      v.filenames = c("data/revision/methylation/B73_methylation_CG.bw", 
                                                      "data/revision/methylation/Mo17_methylation_CG.bw",
                                                      
                                                      "data/revision/methylation/B73_methylation_CHH.bw",
                                                      "data/revision/methylation/Mo17_methylation_CHH.bw",
                                                      
                                                      "data/revision/methylation/B73_methylation_CHG.bw",
                                                      "data/revision/methylation/Mo17_methylation_CHG.bw"),
                                      v.formats = c("bigWig", "bigWig", "bigWig", "bigWig", "bigWig", "bigWig"),
                                      v.parents = c("B73", "Mo17", "B73", "Mo17", "B73", "Mo17"),
                                      v.variants = c("CG", "CG", "CHH", "CHH", "CHG", "CHG"),
                                      v.has_ranges = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
                                      n.bin_width = 10,
                                      th.distance_to_ASB = 2000,
                                      th.high_affinity = 0.85,
                                      n.chromosomes = 10,
                                      group = "all",
                                      n.cpus = 1,
                                      folder = folder_tmp)




# TODO: plot distributions between ASBs and bgSNPs
# TODO: store tmp 

# TODO: what is needed to run this? 
# ASBs, bgSNPs (no genomic location needed!!)
# TODO: bind by additional column c("CpG", "CHG", "CHH")


# TODO: load and add as columns or whatever

# saveRDS(l.SNPs_w_chromatin, paste(folder_tmp, "l.SNPs_w_chromatin_ASBs_and_bgSNPs.rds", sep = "/"))



# asb_distance_vs_methylation(df.ASBs)


# df.ASBs # readRDS(paste("tmp/l.bQTL_gene_partitioning_withGeneDistances_backgroundSampled_", timeStamp, ".rds", sep = ""))[[1]]


# To identify ASBs which overlapped with significant variation in either CG, CHG, or CHH
# methylation between B73 and Mo17, we first, per ASB, assigned averaged methylation
# levels of Mo17 and B73 methylation levels (separately for the CG, CHG, or CHH
#  methylation datasets) within a given window of +/- 40 bp around the ASB position.
# Differentially methylated alleles were defined as described previously 13 . Accordingly, we
# identified ASBs as overlapping with differentially methylated regions if in the B73 or Mo17
# methylation datasets, the methylation level of one allele would be >= 70% while the level
# of the corresponding allele was <= 10%.


# TODO: create subdirectories for plots !!!
methylation_occupancy(l.bQTL_gene_partitioning, 
                      th.distance_to_ASB = 2000,
                      n.chromosomes = 10,
                      n.cpus = 4,
                      n.binWidth = 20,
                      b.val = TRUE,
                      v.species = c("Mo17", "B73"),
                      v.groups = c("all", "genic", "non_genic"),
                      v.species = c("Mo17", "Mo17" , "B73", "B73"),
                      v.stringency <- c("< 0.5 | > 0.5", "< 0.15 | > 0.85", "< 0.5 | > 0.5", "< 0.15 | > 0.85"),
                      v.filenames = c(paste(folder_data,"methylation/B73_CHG.bw", sep = "/"),
                                      paste(folder_data,"methylation/B73_CHH.bw", sep = "/"),
                                      paste(folder_data,"methylation/B73_CpG.bw", sep = "/"),
                                      paste(folder_data,"methylation/MO17_CHG.bw", sep = "/"),
                                      paste(folder_data,"methylation/MO17_CHH.bw", sep = "/"),
                                      paste(folder_data,"methylation/MO17_CpG.bw", sep = "/")),
                      v.datasets = c("B73_CHG", "B73_CHH", "B73_CpG", "MO17_CHG", "MO17_CHH", "MO17_CpG"),
                      v.formats = c("bigWig", "bigWig", "bigWig", "bigWig", "bigWig", "bigWig"),
                      folder_tmp = folder_tmp,
                      folder_output = folder_output,
                      do.plot = F)


# rename methylation, open chromatin?
l.SNPs_w_chromatin <- RDS(paste(folder_tmp, "l.SNPs_w_chromatin_ASBs_and_bgSNPs.rds", sep = "/"))

# Fig. 4 l) Correlation of CG methylation with allele-specific BZR1 binding. Average
# CG methylation of the B73 - Mo17 allele of the 40bp surrounding each ASB are plotted against
# the allelic bias (expressed in percentage of B73 read counts). Differential CG methylation is
# indicated by red dots.
differential_methylation_vs_allelic_bias(l.SNPs=l.bQTL_gene_partitioning,
                                         th.bp_offset = 20,
                                         th.distance_to_ASB = 2000,
                                         n.chromosomes = 10,
                                         n.cpus = 1,
                                         b.val = TRUE,
                                         v.species = c("Mo17", "B73"),
                                         degrees = 15,
                                         th.distance_to_ASB = 5000,
                                         th.padding = 50,
                                         width = 6,
                                         height = 4,
                                         group = "genic_and_nongenic",
                                         v.filenames = c("data/dnase/GSE94291_DNase_ist.bedGraph", 
                                                         "data/methylation/B73_CpG.bw", 
                                                         "data/methylation/MO17_CpG.bw",
                                                         "data/methylation/B73_CHG.bw",
                                                         "data/methylation/MO17_CHG.bw",
                                                         "data/methylation/B73_CHH.bw",
                                                         "data/methylation/MO17_CHH.bw"),
                                         v.datasets = c("DNase", "B73_CpG",  "MO17_CpG", "B73_CHG", "MO17_CHG", "B73_CHH", "MO17_CHH"),
                                         v.formats = c("bedGraph", "bigWig", "bigWig", "bigWig", "bigWig", "bigWig", "bigWig"),
                                         v.has_ranges = c("yes",  "no", "no",  "no", "no",  "no", "no")
)


# TODO: get absolute and differential figure ... 


l.SNPs_w_chromatin <- readRDS(paste(folder_tmp, "l.SNPs_w_chromatin_ASBs_and_bgSNPs.rds", sep = "/"))


df.ASBs



l.SNPs=l.bQTL_gene_partitioning
th.bp_offset = 20
th.distance_to_ASB = 2000
n.chromosomes = 10
n.cpus = 1
b.val = TRUE
v.species = c("Mo17", "B73")
degrees = 15
th.distance_to_ASB = 5000
th.padding = 50
width = 6
height = 4
group = "genic_and_nongenic"
v.filenames = c("data/dnase/GSE94291_DNase_ist.bedGraph", 
                "data/methylation/B73_CpG.bw", 
                "data/methylation/MO17_CpG.bw",
                "data/methylation/B73_CHG.bw",
                "data/methylation/MO17_CHG.bw",
                "data/methylation/B73_CHH.bw",
                "data/methylation/MO17_CHH.bw")
v.datasets = c("DNase", "B73_CpG",  "MO17_CpG", "B73_CHG", "MO17_CHG", "B73_CHH", "MO17_CHH")
v.formats = c("bedGraph", "bigWig", "bigWig", "bigWig", "bigWig", "bigWig", "bigWig")
v.has_ranges = c("yes",  "no", "no",  "no", "no",  "no", "no")


# 
# Fig. 2. Variation in DNA sequence and methylation underlie differential BZR1 binding. a)
# Schematics of HASCh-seq approach and possible causes for allele-specific binding events.
# Chromatin-IP is performed in F1 hybrid plants. From top to bottom, three possible scenarios for
# TF binding to the parental genomes (green and blue) are depicted, where binding strength is
# represented by the line width of the black arrows: when no alterations in motif and chromatin
# structure occurs, binding is expected to be equal. Lower binding is expected If the binding motif
# is altered by a SNP or epigenetic features like DNA methylation for one of the alleles. b) An
# example of allele-specific binding of ZmBZR1 to target gene Zm00001d031434. BZR1 bound
# target reads that uniquely map to either B73 (green) or Mo17 (blue) are shown. c) Allelic and
# spatial distribution of allele-specific BZR1 binding sites with significant bias towards B73 (green)
# or Mo17 (blue) along the maize B73 reference genome. Allelic bias is expressed as percentage
# of B73 read counts. Chromosome borders and length are depicted by dashed lines and arrows at
# the bottom, respectively. Centromere locations are indicated by orange rectangles. A red box on
# 26chromosome 1 marks the ASB upstream of Zm00001d031434 displayed in Fig. 2b. d) Genomic
# distribution of ASBs. 3297 ASBs were classified according to their location relative to genes. In
# case of two genes in the proximity of an ASB, the priority given was exon>intron>UTR>1 kb
# upstream>1 kb downstream>1-5 kb upstream. e) Frequency of BRRE (CGTG[C/T]G), G-Box
# (CACGTG), and the control motif GTACGG (SBP-box 12 ) around ASBs of the alleles with higher
# BZR1 binding. f) Fraction of ASBs overlapping with motifs, for which the allele with canonical
# BRRE, G-Box or a control motif GCCGCC (GCC-box 12 ), showed higher ZmBZR1 affinity. Both
# BR-related BRRE and G-box motifs, but not the control motif, diverge significantly (p<0.001,
# Fisher’s exact test) from the expected 50% random distribution. g) Distribution of lead SNPs of
# ASBs within motifs. Among ASBs which overlapped with BRRE motif, SNPs were enriched in the
# core CG bases. h) ASBs affecting BZR1 motifs and/or overlapping with differentially methylated
# (CpG, CHG or CHH) sites between B73 and Mo17. i)-k) Average CpG (h), CHG (i) and CHH(j)
# methylation frequency in B73 (green) and Mo17 (blue) over ASB loci with a least 85% binding
# bias towards B73 or Mo17. High affinity ( ___ ) and low affinity ( …. ) alleles are considered separately
# for each genotype. l) Correlation of CG methylation with allele-specific BZR1 binding. Average
# CG methylation of the B73 - Mo17 allele of the 40bp surrounding each ASB are plotted against
# the allelic bias (expressed in percentage of B73 read counts). Differential CG methylation is
# indicated by red dots.







message("------------ Motifs ---------")
# To identify ASBs which may be explained by motif variation (Fig. 2h), we extracted the
# +/- 5 bp of the high affinity BZR1 bound allele surrounding ASBs. Using R we scanned
# those 11 bp fragments for canonical BRRE (CGTG[T/C]G, C[G/A]CACG), allowing a
# single base pair mismatch outside the core motif (CGTG), or G-box (CACGTG) motifs
# and determined ASBs where the SNP changed a BRRE or G-box motif into an altered
# (non BRRE or G-box) motif. 

source("utils/motifs_revision.R")

if(!b.load_motif_analysis){
  source("motifs/predefined_motif_analysis.R")
  
  
  
  predefined_motif_analysis(l.bQTL_gene_partitioning, df.peaks, motifs, v.motif_offset)
  
  l.motif_analysis <- readRDS(paste(folder_tmp, "l.motif_analysis.rds", sep = "/"))
  l.nucleotideInPeaks <- readRDS(paste(folder_tmp, "l.nucleotideInPeaks.rds", sep = "/"))
  
  
  source("motifs/predefined_motifs_in_peaks.R")
  predefined_motifs_in_peaks(l.motif_analysis, l.nucleotideInPeaks)
  
  l.motif_analysis.postprocessed = readRDS(paste(folder_tmp, "l.motif_analysis.postprocessed.rds", sep = "/"))
  
  
  source("motifs/directionality_and_snp_distribution.R")
  directionality_and_snp_distribution(l.motif_analysis.postprocessed)
  
  df.motif_nucleotidePositionPartitions <- readRDS(paste(folder_tmp, "df.motif_nucleotidePositionPartitions.rds", sep = "/"))
  
  
  source("motifs/motif_snp_position_significance.R")
  motif_snp_position_significance(l.bQTL_gene_partitioning, motifs, v.motif_offset, df.peaks, genome, l.genome.mutant, n.chromosomes)
  
  
  source("motifs/evaluating_BZR1_core_motif.R")
  evaluating_BZR1_core_motif(df.ASBs, core_motif = "CGTG", s.motif_offset = 5, n.chromosomes = 10)
  
}else{
  l.motif_analysis <- readRDS(paste(folder_tmp, "l.motif_analysis.rds", sep = "/"))
  l.nucleotideInPeaks <- readRDS(paste(folder_tmp, "l.nucleotideInPeaks.rds", sep = "/"))
  l.motif_analysis.postprocessed = readRDS(paste(folder_tmp, "l.motif_analysis.postprocessed.rds", sep = "/"))
  df.motif_nucleotidePositionPartitions <- readRDS(paste(folder_tmp, "df.motif_nucleotidePositionPartitions.rds", sep = "/"))
  df.motifPositionAnalysis <- readRDS(paste(folder_tmp, "df.motifPositionAnalysis.rds", sep = "/"))
  df.ASBs.core_motifs <- readRDS(paste(folder_tmp, "df.ASBs.core_motifs.rds", sep = "/"))
}

if(b.write_paper){
  write.table(l.motif_analysis[[1]], paste(folder_output, "/motifs/S0_1.txt", sep = ""),  row.names = FALSE, quote = FALSE, sep ="\t")
  
  print(table(l.motif_analysis.postprocessed[[1]]$direction))
  write.table(l.motif_analysis.postprocessed[[1]], paste(folder_output, "/motifs/S0_2.txt", sep = ""),  row.names = FALSE, quote = FALSE, sep ="\t")
  write.table(df.motif_directionality, paste(folder_output, "/motifs/S0_3.txt", sep = ""),  row.names = FALSE, quote = FALSE, sep ="\t")
  write.table(df.motif_nucleotidePositionPartitions, paste(folder_output, "/motifs/S0_4.txt", sep = ""),  row.names = FALSE, quote = FALSE, sep ="\t")
  write.table(df.motifPositionAnalysis, paste(folder_output, "/motifs/S0_5.txt", sep = ""),  row.names = FALSE, quote = FALSE, sep ="\t")
  write.table(df.ASBs.core_motifs, paste(folder_output, "/SX3.txt", sep = ""), quote = FALSE, row.names = FALSE, sep ="\t")
  
}


