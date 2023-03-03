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
