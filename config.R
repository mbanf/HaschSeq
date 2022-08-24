b.firstRun = FALSE 
n.chromosomes <- 10

b.peak_analysis <- TRUE
b.disequilibriumLinkageHandling <- TRUE
b.load_filtered_binding_peaks = TRUE
b.doInitialHomozygousFilter <- FALSE
b.save_intermediate_results <- TRUE
b.load_motif_analysis <- TRUE
b.load_chip_analysis <- TRUE

minReadDepth = 50
th.p.bQTL <- 0.001
s.disequilibriumDistance <- 75
b.criteria_linkageDisequilibrium = "pval"

n.cpus <- 1

# peak width 
s.half_window_size <- 100
s.half_window_size.input <- 201

#### 

v.partitions <- c("promoter_5kb", "promoter_1kb", "gene", "five_prime_UTR", "exon", "intron", "three_prime_UTR", "post_gene_1kb", "non_genic") 
s.multiplyer <- 50

# population sampling
n.bgSnp.multiplier <- 100
b.load_snp_to_gene_partitioning_and_background_sampling = TRUE
b.duplicateRemoval <- FALSE
b.loadGenePartitioning <- FALSE


#"BRRE1", "BRRE1_RC", "BRRE2","BRRE2_RC", "GBOX" => nur diese sind targets (fuer ASB explanation), others are controls
v.motifs = c("CGTGCG", "CGCACG", 
           "GTACGG", "CCGTAC",
           "GCCGCC","CGGCGG",
           "CGTGTG","CACACG", "CACGTG")
names(v.motifs) = c("BRRE1", "BRRE1_RC", "SBP","SBP","OTHER", "OTHER","BRRE2","BRRE2_RC", "GBOX")
v.motif_offset = c(5,5,5,5,5,5,5,5,5)
v.motifs = toupper(v.motifs)

# GWAS
s.dist_ASB_to_GWAS = 2000

# 
# # nucleotide distribution
# prior.params = c(0.237493, 0.262507, 0.262507, 0.237493)
# names(prior.params) <- c("A","C", "G", "T")
# 
# 
# motifs <- toupper(c("cgtgcg", "cacgtg", "cgcacg", "cgtgca", "tgcacg", "cagcag", "ctgctg", "cgtgtg",
#                     "aaaaaa", "tttttt", "gccgcc", "ggcggc", "gtgcgg", "ccgcac")) 
#  
# 
# motifs.BRRE_GBOX <- c("cgtgcg", "cgtgtg", "cacgtg", "cgcacg","cacacg")
# # motifs.BRRE_GBOX <- c("cgtgcg","cgtgtg", "cacgtg", "cgcacg","cgcaca")
# 
# # Parameters
# s.half_window_size <- 100
# s.half_window_size.input <- 201
# 
# s.multiplyer <- 50
# 
# timeStamp <- 0117
# 
# 
# 
# n.cpus <- 2
# verbose <- FALSE
# 
# s.AGPv3_or_AGPv4 = "AGPv4"
# 
# # define background size
# n.bgSnp.multiplier <- 100
# 
# b.doInitialHomozygousFilter <- FALSE
# b.disequilibriumLinkageHandling <- TRUE
# b.denovo_motif_analysis <- TRUE
# b.predefined_motif_analysis <- TRUE
# 
# minReadDepth = 50
# th.p.bQTL <- 0.001
# s.disequilibriumDistance <- 75
# 
# 
# b.load.background_with_gene_partitioning = TRUE
# b.criteria_linkageDisequilibrium = "pval"
# 
# b.load_filtered_binding_peaks = TRUE
# b.load_ChipSeqB73_filtered_peaks = TRUE
# b.load_snp_to_gene_partitioning_and_background_sampling = TRUE
# 
# 
# s.dist_ASB_to_GWAS = 2000
# 
# 
# 
# #### dataset definition 
# 
# 
# v.filenames <- c("datasets_paper/Enhancer_HM/GSE94251_H3K9ac_ist.bedGraph", 
#                  "datasets_paper/Enhancer_HM/GSE94291_DNase_ist.bedGraph", 
#                  "datasets_paper/Methylation/B73_CHG.bw", 
#                  "datasets_paper/Methylation/B73_CHH.bw", 
#                  "datasets_paper/Methylation/B73_CpG.bw", 
#                  "datasets_paper/Methylation/MO17_CHG.bw", 
#                  "datasets_paper/Methylation/MO17_CHH.bw", 
#                  "datasets_paper/Methylation/MO17_CpG.bw")
# 
# v.datasets <- c("H3K9", "DNase",  "B73_CHG", "B73_CHH", "B73_CpG", "MO17_CHG", "MO17_CHH", "MO17_CpG")
# v.formats <- c("bedGraph", "bedGraph", "bigWig", "bigWig", "bigWig", "bigWig", "bigWig", "bigWig")
# v.has_ranges = c("yes", "yes", "no", "no", "no", "no", "no", "no")
# 
# 
# l.SNPs <- readRDS(paste("tmp/l.bQTL_gene_partitioning_withGeneDistances_backgroundSampled_", timeStamp, ".rds", sep = ""))
# 
