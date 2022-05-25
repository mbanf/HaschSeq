message("loading binding QTL datasets...")

strt<-Sys.time()
postTotal = read.csv(path.bQTLs, header=T, stringsAsFactors = FALSE)
df.bQTLsInput = read.csv(path.bQTLsInput, header=T, stringsAsFactors = FALSE)
print(Sys.time() - strt)
#postTotal = read.table(pd, sep="\t", header=T, stringsAsFactors = FALSE)

df.bQTLsInput <- subset(df.bQTLsInput, df.bQTLsInput$refCount > 0 & df.bQTLsInput$altCount > 0)
df.bQTLsInput["id"] <- paste(df.bQTLsInput$contig, "_", df.bQTLsInput$position, sep = "")
postTotal["id"] <- paste(postTotal$contig, "_", postTotal$position, sep = "")

# store all snps 
df.bQTLsOutputformat <- df.bQTLsInput[,c("id","contig", "position", "refAllele", "altAllele")]

# input numbers
n.bQTLs_readDepth50 <- nrow(postTotal)
n.bQTLs_readDepth5 <- nrow(df.bQTLsInput)

print(n.bQTLs_readDepth50)
print(n.bQTLs_readDepth5)

# re-insert ASBs with at least minReadsPerAllele reads in ref and alt (based on the input with presence in at least 5 in 6 repeats)

# blacklist areas 
minReadsPerAllele = minReadDepth * 0.1 
postTotalSafe <- subset(postTotal, postTotal$refCount >= minReadsPerAllele & postTotal$altCount >= minReadsPerAllele)
postTotalCheck <- subset(postTotal, !postTotal$id %in% postTotalSafe$id)
postTotalCheck <- subset(postTotalCheck, postTotalCheck$id %in% df.bQTLsInput$id)
postTotal <- rbind(postTotalSafe, postTotalCheck)

# if(b.save_intermediate_results){
#   write.csv(postTotal, "paper_tmp_files/bQTL_after_input.csv")
# }
# 

# filtering
n.bQTLs_readDepth50.both_alleles_10Percent_bias <- nrow(postTotal)
n.bQTLs_readDepth5.both_alleles_nonZero_cound <- nrow(df.bQTLsInput)

print(n.bQTLs_readDepth50.both_alleles_10Percent_bias)
print(n.bQTLs_readDepth5.both_alleles_nonZero_cound)

message("input bQTL after min reads per allele filter: ", nrow(postTotal))

v.postFreq <- postTotal$refCount / postTotal$totalCount
v.postFreq <- v.postFreq[which(v.postFreq != 0)]
v.postFreq <- v.postFreq[which(v.postFreq != 1)]
prob.bias <- median(v.postFreq)

message("bias probability: ", prob.bias)

#p = 0.5 # 0.551948 # expectancy
p = 1 - prob.bias # 0.5120976 # post frequency - ref / all
test <- sapply(1:nrow(postTotal), function(i) binom.test(as.integer(postTotal$altCount[i]), as.integer(postTotal$totalCount[i]), p = p)$p.value) #  analytic.pv(preProps[i],preVars[i],postProps[i],depths[i]))
fdr <- p.adjust(test, "bonferroni")
postTotal["p-value (corrected)"] <- fdr

postTotal <- subset(postTotal, postTotal$totalCount >= minReadDepth)

postTotal["POSTfreq"] <- postTotal$refCount / postTotal$totalCount
postTotal["POSTallele"] <- ifelse(postTotal$altCount > postTotal$refCount, postTotal$altAllele, postTotal$refAllele)

# postTotal <- subset(postTotal, postTotal$POSTfreq < 1 & postTotal$POSTfreq > 0)
# median(postTotal$POSTfreq)

# message("make figure per paper")
df.bQTLsInput["POSTfreq"] <- df.bQTLsInput$refCount / df.bQTLsInput$totalCount
# df.bQTLsInput <- subset(df.bQTLsInput, df.bQTLsInput$totalCount > 5)

# perform chromosome based - blacklisting

#message("input bQTL: ", nrow(postTotal))
#bQTL_scatterplot(postTotal=postTotal)

## manual artifact removal
message("remove manually labelled artifact regions - further reducing significant ")
# make figure 
# bQTL_scatterplot_chr(postTotal=postTotal, chr = 3)
chr = 3
cut.left <-  186240476
cut.right <- 210080072
print(abs(cut.left - cut.right) / 1e6)

idx.cut <- which(df.bQTLsInput$contig == chr & df.bQTLsInput$position < cut.right & df.bQTLsInput$position > cut.left)
df.bQTLsInput <- df.bQTLsInput[!seq(1:nrow(df.bQTLsInput)) %in%  idx.cut,]

chr = 10
cut.left <- 2965337
cut.right <- 4124958
print(abs(cut.left - cut.right) / 1e6)

idx.cut <- which(df.bQTLsInput$contig == chr & df.bQTLsInput$position < cut.right & df.bQTLsInput$position > cut.left)
df.bQTLsInput <- df.bQTLsInput[!seq(1:nrow(df.bQTLsInput)) %in%  idx.cut,]

# double check
# bQTL_scatterplot(postTotal=df.bQTLsInput, chr = 10)

postTotal <- subset(postTotal, postTotal$id %in% df.bQTLsInput$id)
postTotal <- postTotal[,!names(postTotal) %in% "id"]

message("input bQTL after min reads per allele filter: ", nrow(postTotal))

# if(b.save_intermediate_results){
#   write.csv(postTotal, "paper_tmp_files/bQTL_after_input_blacklisting.csv")
# }