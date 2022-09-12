### specialized re-annotation for GWAS / VCAP analysis (compensate for missing UTR in Mo17 annotation) ### 

### for VCAP we need as well UTR annotations for both parents for ASBs and background SNPs ### 
### Hence, we need to resample bgSNPs, and use B73 annotations (chr and pos) only ### 


source("utils/genomic_location_revision.R")

df.ASBs <- readRDS(paste(folder_tmp, "6143_ASBs_w_genomic_location.rds", sep = "/"))

df.ASBs_genomic_location = add_genomic_location_bQTLs(df.ASBs,  
                                                      df.gene_annotation,
                                                      n.cpus = 5,
                                                      b.low_affinity = F,
                                                      b.B73_only = T)


df.snps_genomic_location = add_genomic_location_bQTLs(df.bgSNPs.candidates,  
                                                      df.gene_annotation,
                                                      n.cpus = 5,
                                                      b.low_affinity = F,
                                                      b.B73_only = T)

# saveRDS(df.ASBs_genomic_location, paste(folder_tmp, "6143_ASBs_w_genomic_location.B73.rds", sep = "/"))   
# saveRDS(df.snps_genomic_location, paste(folder_tmp, "6181870_bgSNPs_candidates_w_genomic_location.B73.rds", sep = "/"))

df.ASBs_genomic_location <- readRDS(paste(folder_tmp, "6143_ASBs_w_genomic_location.B73.rds", sep = "/"))
df.snps_genomic_location <- readRDS(paste(folder_tmp, "6181870_bgSNPs_candidates_w_genomic_location.B73.rds", sep = "/"))


source("utils/qtl_background_sampling.R")

v.partitions <- c("promoter_5kb", "promoter_1kb", "gene", "five_prime_UTR", "exon", "intron", "three_prime_UTR", "post_gene_1kb", "non_genic") 

res <- create_background_QTLs(df.ASBs_genomic_location,
                              df.snps_genomic_location, 
                              n.bg.multiplier = 40,  # reduced multiplyer
                              v.partitions = v.partitions, 
                              seed = 1234)

df.bgSNPs <- res$df.bgSNPs
message(nrow(df.bgSNPs), " sampled, non-significant (background) SNPs located outside high-confidence BZR1 binding peaks (all 6 replicated merged)")
hist(df.bgSNPs$POSTfreq, breaks = 100, main = paste("Allelic bias of", nrow(df.bgSNPs), "sampled, non-significant (background) SNPs located outside high-confidence BZR1 binding peaks"), 
     xlab = "allelic bias", cex.main = 0.7)

genomic_distribution(df.ASBs_genomic_location, df.bgSNPs, v.partitions)

bQTL_scatterplot(df.bgSNPs)

# write.csv(df.ASBs_genomic_location, "6143_ASBs_w_genomic_location.B73.csv", quote = F, row.names = F)
# write.csv(df.bgSNPs, "254017_background_SNPs_w_genomic_location.B73.csv", quote = F, row.names = F)
