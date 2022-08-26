# install.packages( paste(folder_data,"BSgenome.Zmays.NCBI.AGPv4/", sep ="/"), repos = NULL, type="source")

#require("BSgenome.Zmays.NCBI.AGPv4")
#genome <- BSgenome.Zmays.NCBI.AGPv4


# only once - assemble the mutant genome - define output directory 
#v.mo17_files <- c("chr1.vcf","chr2.vcf","chr3.vcf","chr4.vcf","chr5.vcf","chr6.vcf","chr7.vcf","chr8.vcf","chr9.vcf","chr10.vcf")
#path.snptables <- paste(folder_data, "mo17_snptables/", sep ="/")


#path.bQTLs <- paste(folder_data,"bQTL/C1_C22_YFP_mindepth50_q13_BaseQ13", sep ="/") #"C1_C22_YFP_q20.30.POSTth.txt"
#path.bQTLsInput <- paste(folder_data,"bQTL/C1_C22_Input_mindepth5_q13_BaseQ13", sep ="/") #"C1_C22_YFP_q20.30.POSTth.txt" - control - alle dna 50 - 50 % => hole genau die mit tf

#path.bindingPeaks <- paste(folder_data,"GEM/C1_C22_YFP/q5_10fold/q5_10fold.GEM_events.txt", sep ="/")

# TODO: where are these? - general questions
# l.path.bindingPeaksAll <- vector(mode = "list", length = 6)
# l.path.bindingPeaksAll[[1]] ="/Volumes/Shared/Everyone/Michael_Thomas/AGPv4/GEM/Output/BZR1_1_S1_L002_R/q2_5fold/q2_5fold.GEM_events.bed"
# l.path.bindingPeaksAll[[2]] ="/Volumes/Shared/Everyone/Michael_Thomas/AGPv4/GEM/Output/BZR1_2_S2_L002_R/q2_5fold/q2_5fold.GEM_events.bed"
# l.path.bindingPeaksAll[[3]] ="/Volumes/Shared/Everyone/Michael_Thomas/AGPv4/GEM/Output/BZR1_3_S3_L002_R/q2_5fold/q2_5fold.GEM_events.bed"
# l.path.bindingPeaksAll[[4]] ="/Volumes/Shared/Everyone/Michael_Thomas/AGPv4/GEM/Output/BZR1_4_S4_L002_R/q2_5fold/q2_5fold.GEM_events.bed"
# l.path.bindingPeaksAll[[5]] ="/Volumes/Shared/Everyone/Michael_Thomas/AGPv4/GEM/Output/BZR1_5_S5_L002_R/q2_5fold/q2_5fold.GEM_events.bed"
# l.path.bindingPeaksAll[[6]] ="/Volumes/Shared/Everyone/Michael_Thomas/AGPv4/GEM/Output/BZR1_6_S6_L002_R/q2_5fold/q2_5fold.GEM_events.bed"


###


message("loading gene annotation datasets...")
strt<-Sys.time()

#path.gene_function<- paste(folder_data,"GeneAnnotation/B73v4_gene_function.txt", sep ="/")
path.geneID_conversion <- paste(folder_data,"GeneAnnotation/maize.v3TOv4.geneIDhistory.txt", sep ="/")
#df.gene_function <- read.table(path.gene_function, header = FALSE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE, quote = "")
#names(df.gene_function) <- c("gene.ID", "gene.function")

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



"Parent="

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


# gene model weg, 
# unique filter ... 



# 
# 
# 
# # install.packages("BSgenome.Zmays.NCBI.AGPv4/", repos = NULL, type="source")
# 
# # genome as delivered - 
# if(s.AGPv3_or_AGPv4 == "AGPv3"){
#   require("BSgenome.Zmays.NCBI.AGPv3")
#   genome <- BSgenome.Zmays.NCBI.AGPv3
# }else{
#   require("BSgenome.Zmays.NCBI.AGPv4")
#   genome <- BSgenome.Zmays.NCBI.AGPv4
# }
# 
# 
# # group into specific categories 
# # check if folders are empty 
# 
# # ath chip seq - dark
# df.At_BZR1_ChipSeqTargets <- read.table("datasets_paper/ArabidopsisScriptsAndDatasets/At_BZR1_ChipSeqTargets.txt", header = TRUE, sep ="\t", quote = "", stringsAsFactors = FALSE)
# names(df.At_BZR1_ChipSeqTargets) <- "locus"
# 
# # ath chIP chip - light
# df.At_BZR1_ChIPchipTargets <- read.table("datasets_paper/ArabidopsisScriptsAndDatasets/At_BZR1_ChIPchip.txt", header = TRUE, sep ="\t", quote = "", stringsAsFactors = FALSE)
# 
# path.bindingPeaks <- "datasets_paper/GEM/C1_C22_YFP/q5_10fold/q5_10fold.GEM_events.txt"
# 
# l.path.bindingPeaksAll <- vector(mode = "list", length = 6)
# l.path.bindingPeaksAll[[1]] ="/Volumes/Shared/Everyone/Michael_Thomas/AGPv4/GEM/Output/BZR1_1_S1_L002_R/q2_5fold/q2_5fold.GEM_events.bed"
# l.path.bindingPeaksAll[[2]] ="/Volumes/Shared/Everyone/Michael_Thomas/AGPv4/GEM/Output/BZR1_2_S2_L002_R/q2_5fold/q2_5fold.GEM_events.bed"
# l.path.bindingPeaksAll[[3]] ="/Volumes/Shared/Everyone/Michael_Thomas/AGPv4/GEM/Output/BZR1_3_S3_L002_R/q2_5fold/q2_5fold.GEM_events.bed"
# l.path.bindingPeaksAll[[4]] ="/Volumes/Shared/Everyone/Michael_Thomas/AGPv4/GEM/Output/BZR1_4_S4_L002_R/q2_5fold/q2_5fold.GEM_events.bed"
# l.path.bindingPeaksAll[[5]] ="/Volumes/Shared/Everyone/Michael_Thomas/AGPv4/GEM/Output/BZR1_5_S5_L002_R/q2_5fold/q2_5fold.GEM_events.bed"
# l.path.bindingPeaksAll[[6]] ="/Volumes/Shared/Everyone/Michael_Thomas/AGPv4/GEM/Output/BZR1_6_S6_L002_R/q2_5fold/q2_5fold.GEM_events.bed"
# 
# l.path.ChipSeqB73 <- vector(mode = "list", length = 3)
# l.path.ChipSeqB73[[1]] ="datasets_paper/ChipSeq/B73_ChIP/AGTCAA/q2_5fold/q2_5fold.GEM_events.txt"
# l.path.ChipSeqB73[[2]] ="datasets_paper/ChipSeq/B73_ChIP/ATGTCA/q2_5fold/q2_5fold.GEM_events.txt"
# l.path.ChipSeqB73[[3]] ="datasets_paper/ChipSeq/B73_ChIP/CCGTCC/q2_5fold/q2_5fold.GEM_events.txt"
# 
# path.bQTLs <- "datasets_paper/bQTL/C1_C22_YFP_mindepth50_q13_BaseQ13" #"C1_C22_YFP_q20.30.POSTth.txt"
# path.bQTLsInput <- "datasets_paper/bQTL/C1_C22_Input_mindepth5_q13_BaseQ13" #"C1_C22_YFP_q20.30.POSTth.txt" - control - alle dna 50 - 50 % => hole genau die mit tf
# path.bQTLs.parentalLines <- c("datasets_paper/parental_origin/B73xMo17_mindepth25_q13_BaseQ13", "/datasets_paper/parental_origin/Mo17xB73_mindepth25_q13_BaseQ13")
# path.eQTLs <- "datasets_paper/eQTL/AGPv4_uplifted_Liu_et_al_2016_eQTL_AGPv2_mmc2.bed"
# path.phenotype_gwas.original <- "datasets_paper/GWAS/AGPv4_uplifted_Wallace_etal_2014_PLoSGenet_GWAS_hits-150112.txt"
# path.MNase <- "datasets_paper/MNase/AGPv4_uplifted.AP.bfthresh1.1.MNaseHS.Ranges.bed"
# 
# path.gene_function<- "datasets_paper/GeneAnnotation/B73v4_gene_function.txt"
# path.geneID_conversion <- "datasets_paper/GeneAnnotation/maize.v3TOv4.geneIDhistory.txt"
# path.gff <- "datasets_paper/GeneAnnotation/Zea_mays.AGPv4.36.chr.gff3"
# 
# path.gene_orthologs <- "datasets_paper/GeneAnnotation/Zmays_284_Ensembl-18_2010-01-MaizeSequence.annotation_info.txt"
# 
# 
# # only once - assemble the mutant genome - define output directory 
# v.mo17_files <- c("chr1.vcf","chr2.vcf","chr3.vcf","chr4.vcf","chr5.vcf","chr6.vcf","chr7.vcf","chr8.vcf","chr9.vcf","chr10.vcf")
# 
# # predefined motifs to be used
# motifs <- toupper(c("cgtgcg", "cacgtg", "cgcacg", "cgtgca", "tgcacg", "cagcag", "ctgctg", "cgtgtg",
#                     "aaaaaa", "tttttt", "gccgcc", "ggcggc", "gtgcgg", "ccgcac")) 
# 
# meme_non_motifs <- toupper(c("cgtgcg", "cacgtg", "cgcacg")) 
# 
# vec.chroms.snp <- c("chr1", "chr2", "chr3", "chr4", "chr5","chr6","chr7","chr8", "chr9", "chr10")
# 
# path.snptables <- "datasets_paper/mo17_snptables/"
# 
# df.enhancer_candidates <- read.table("datasets_paper/other/HUSK_enhancer_candidates.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
# df.enhancer_candidates <- rbind(df.enhancer_candidates, read.table("datasets_paper/other/V2_IST_enhancer_candidates.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE))
# 
# df.enhancer_genes <- read.table("datasets_paper/other/Enhances_constutitive_expressed_specific_genes.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
# df.enhancer_genes <- rbind(df.enhancer_genes, read.table("datasets_paper/other/Enhances_husk_specific_genes.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE))
# df.enhancer_genes <- rbind(df.enhancer_genes, read.table("datasets_paper/other/Enhances_V2_IST_specific_genes.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE))
# 
# 
# df.gene_orthologs = read.table(path.gene_orthologs, header=T, stringsAsFactors = FALSE, sep = "\t", fill = T)
# 
# 
# 
# ###? 

# # TODO: make unique - how the tmp stuff has been produced
# message("preparing mutant genome")
# 
# df.snp_positions <- c()
# l.genome.mutant <- vector(mode = "list", length = n.chromosomes)
# # welche haben laenge
# for(i in 1:n.chromosomes){
# 
#   print(paste("processing chromosome ", i))
# 
#   df.snp.pos <- readRDS(paste("tmp/df.snp.pos_",i, ".rds")) # TODO: fix version in onedrive etc... 
# 
#   # n.snps <- n.snps + length(df.snp.pos)
# 
#   # remove heterozygote snps
#   # idx.snps <- which(lapply(df.snp.pos@elementMetadata$ALT, length) == 1)
# 
#   # A - generic SNPs
#   #vec.snp.pos <- as.numeric(df.snp.pos@ranges@start[idx.snps])
#   vec.snp.pos <- as.numeric(df.snp.pos@ranges@start)
#   # vec.snp.bases <- as.character(unlist(df.snp.pos@elementMetadata$ALT)) #
# 
#   # vec.snp.bases <- as.character(unlist(lapply(df.snp.pos@elementMetadata$ALT, function(m) m[1])))
#   vec.ref.bases <- as.character(df.snp.pos@elementMetadata$REF)
# 
#   vec.snp.bases <- CharacterList(df.snp.pos@elementMetadata$ALT)
#   vec.snp.bases = unstrsplit(vec.snp.bases, sep = ",")
#   vec.snp.bases <- gsub("\\,.*", "", vec.snp.bases)
# 
#   #vec.snp.bases <- as.character(unlist(df.snp.pos@elementMetadata$ALT[idx.snps]))
# 
#   vec.snp.bases <- ifelse(vec.snp.bases == "<DEL>", "N",vec.snp.bases) # both represented separately
#   vec.snp.bases <- ifelse(vec.snp.bases == "<INS>", "N",vec.snp.bases)
# 
#   # remove empty strings #
#   idx.snp.exceptions <- which(vec.snp.bases == "")
# 
#   if(length(idx.snp.exceptions) > 0){
#     vec.snp.pos <- vec.snp.pos[-idx.snp.exceptions]
#     vec.snp.bases <- vec.snp.bases[-idx.snp.exceptions]
#     vec.ref.bases <- vec.ref.bases[-idx.snp.exceptions]
#   }
# 
#   df.snp_positions.i <- data.frame(ID = paste(i ,"_", vec.snp.pos, sep ="") , chromosome = rep(paste(i ,"_maternal", sep ="") , length(vec.snp.pos)) ,
#                                    position = vec.snp.pos, strand =  rep(1, length(vec.snp.pos)),
#                                    Ref_SNP = paste(vec.ref.bases ,"/", vec.snp.bases, sep =""))
# 
#   df.snp_positions <- rbind(df.snp_positions, df.snp_positions.i)
# 
#   genome.reference <- DNAString(genome[[i]])
#   genome.mutant    <- replaceLetterAt(genome.reference, vec.snp.pos, vec.snp.bases) # vergleich vor austausch
# 
#   l.genome.mutant[[i]] <- genome.mutant
# 
# 
#   #  v.sequences <- unlist(lapply(l.sequences, function(m) {m[[1]]}))
#   if(FALSE){
#     Sequences = DNAStringSet(genome.mutant)
#     #Sequences <- read.DNAStringSet(FastaFile, "fasta")
#     names(Sequences) <- as.character(paste("chr",i,sep = ""))
#     writeXStringSet(Sequences, paste("output/mo17Genome/chr",i,".fasta", sep = ""), format="fasta")
#   }
# 
# }
# 
# if(FALSE){
#   write.table(df.snp_positions,  paste("output/df.snp_positions.csv", sep = ""), col.names = FALSE, row.names = FALSE, sep ="\t")
# }
# 
# 
# 
# # binding peaks
# message("loading binding peak datasets...")
# 
# if(!b.load_filtered_binding_peaks){
#   
#   strt<-Sys.time()
#   df.peaks <- read.table(path.bindingPeaks, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
#   df.peaks["seqnames"] <- as.numeric(unlist(lapply(strsplit(df.peaks$Position, ":"), function(m) m[[1]])))
#   df.peaks["posPeak"] <- as.numeric(unlist(lapply(strsplit(df.peaks$Position, ":"), function(m) m[[2]])))
#   df.peaks["start"] <- df.peaks$posPeak - s.half_window_size
#   df.peaks["end"] <- df.peaks$posPeak + s.half_window_size
#   
#   l.df.peaks.all <- vector(mode = "list", length = length(l.path.bindingPeaksAll))
#   for(i in 1:length(l.path.bindingPeaksAll)){
#     l.df.peaks.all[[i]] <- read.table(l.path.bindingPeaksAll[[i]], header = F, sep = "\t", stringsAsFactors = FALSE, fill = TRUE, quote = "")  
#     l.df.peaks.all[[i]] <- l.df.peaks.all[[i]][2:nrow(l.df.peaks.all[[i]]),]
#   }
#   
#   # 6 von 6 
#   l.df.peaks.all[[1]]["left"] <- l.df.peaks.all[[1]]$V2 - s.half_window_size.input
#   l.df.peaks.all[[1]]["right"] <- l.df.peaks.all[[1]]$V3 + s.half_window_size.input
#   
#   if(FALSE){
#     strt<-Sys.time() 
#     #cl<-makeCluster(min(n.chromosomes, n.cpus))
#     #registerDoParallel(cl)
#     # l.peaks.filtered <-  foreach(i = 1:n.chromosomes, .packages=c("seqinr", "VariantAnnotation", "Biostrings")) %dopar% { 
#     df.peaks.filtered <- c()
#     for(i in 1:n.chromosomes){
#       l.peaks.chromosome <- vector(mode = "list", length = 6)
#       l.idx.blacklist <- vector(mode = "list", length = 6)
#       for(l in 1:length(l.path.bindingPeaksAll)){
#         l.peaks.chromosome[[l]] <- subset(l.df.peaks.all[[l]], l.df.peaks.all[[l]]$V1 == paste("chr", i, sep = ""))
#         l.idx.blacklist[[l]] <- numeric(nrow(l.peaks.chromosome[[l]]))
#       }
#       v.idx.keep <- numeric(nrow(l.peaks.chromosome[[1]]))
#       pb <- txtProgressBar(min = 0, max = nrow(l.peaks.chromosome[[1]]), style = 3)
#       for(j in 1:nrow(l.peaks.chromosome[[1]])){
#         setTxtProgressBar(pb, j)
#         n.found <- 0
#         for(l in 2:length(l.peaks.chromosome)){
#           b.found = FALSE
#           for(k in 1:nrow(l.peaks.chromosome[[l]])){
#             
#             if(!k %in% l.idx.blacklist[[l]]){        
#               
#               if( l.peaks.chromosome[[1]]$V1[j] == l.peaks.chromosome[[l]]$V1[k]){
#                  # give me a middle 
#                  if(l.peaks.chromosome[[l]]$V3[k] - 100 > l.peaks.chromosome[[1]]$left[j] & l.peaks.chromosome[[l]]$V3[k] - 100 < l.peaks.chromosome[[1]]$right[j]){
#                     n.found = n.found + 1
#                     b.found = TRUE
#                     l.idx.blacklist[[l]][k] <- k
#                     break
#                  }
#               }
#             }
#           }
#           if(b.found == FALSE){
#             break
#           }
#         }
#         if(n.found == 5){
#           v.idx.keep[j] <- j  
#         }
#       }
#       close(pb)
#     
#       v.idx.keep <- v.idx.keep[v.idx.keep != 0]
#       l.peaks.chromosome[[1]] <- l.peaks.chromosome[[1]][v.idx.keep,]
#       
#       df.peaks.filtered <- rbind(df.peaks.filtered, l.peaks.chromosome[[1]])
#       # l.peaks.chromosome[[1]]
#     }
#     #stopCluster(cl)
#     print(Sys.time()-strt)
#     
#   }else{
#   
#     l.peaks.filtered <- readRDS("tmp/l.peaks.filtered.rds")  
#     
#     df.peaks.filtered <- c()
#     for(i in 1:n.chromosomes){
#       df.peaks.filtered <- rbind(df.peaks.filtered, l.peaks.filtered[[i]])
#     }
#     # 
#     df.peaks.final <- c()
#     for(i in 1:n.chromosomes){
#       
#       df.peaks.chromosome <- subset(df.peaks, df.peaks$seqnames == i)
#       df.peaks.filtered.chromosome <- subset(df.peaks.filtered, df.peaks.filtered$V1 == paste("chr", i, sep = ""))
#       
#       v.idx.keep <- numeric(nrow(df.peaks.chromosome))
#       pb <- txtProgressBar(min = 0, max = nrow(df.peaks.chromosome), style = 3)
#       for(j in 1:nrow(df.peaks.chromosome)){
#         setTxtProgressBar(pb, j)
#         for(k in 1:nrow(df.peaks.filtered.chromosome)){
#           dist <- abs(df.peaks.chromosome$posPeak[j]  -  df.peaks.filtered.chromosome$V3[k] - 100)  
#           if(dist < 500){
#             v.idx.keep[j] <- j
#             break
#           }
#         }
#       }
#       close(pb)
#       
#       v.idx.keep <- v.idx.keep[v.idx.keep>0]
#       df.peaks.chromosome <- df.peaks.chromosome[v.idx.keep, ]
#       df.peaks.final <- rbind(df.peaks.final, df.peaks.chromosome)
#     }
# 
#     if(FALSE){
#       saveRDS(df.peaks.final, "tmp/df.peaks.final.rds")  
#     }
#     
#   }
# }else{
#   df.peaks.final <- readRDS("tmp/df.peaks.final.rds")
# }
# 
# df.peaks <- df.peaks.final
# df.peaks["start"] <- df.peaks$posPeak - s.half_window_size
# df.peaks["end"] <- df.peaks$posPeak + s.half_window_size
# 
# message("number of binding peaks: ", nrow(df.peaks))
# 
# 
# message("loading binding QTL datasets...")
# strt<-Sys.time()
# postTotal = read.csv(path.bQTLs, header=T, stringsAsFactors = FALSE)
# df.bQTLsInput = read.csv(path.bQTLsInput, header=T, stringsAsFactors = FALSE)
# print(Sys.time() - strt)
# #postTotal = read.table(pd, sep="\t", header=T, stringsAsFactors = FALSE)
# 
# df.bQTLsInput <- subset(df.bQTLsInput, df.bQTLsInput$refCount > 0 & df.bQTLsInput$altCount > 0)
# df.bQTLsInput["id"] <- paste(df.bQTLsInput$contig, "_", df.bQTLsInput$position, sep = "")
# postTotal["id"] <- paste(postTotal$contig, "_", postTotal$position, sep = "")
# 
# # store all snps 
# df.bQTLsOutputformat <- df.bQTLsInput[,c("id","contig", "position", "refAllele", "altAllele")]
# 
# # input numbers
# n.bQTLs_readDepth50 <- nrow(postTotal)
# n.bQTLs_readDepth5 <- nrow(df.bQTLsInput)
# 
# print(n.bQTLs_readDepth50)
# print(n.bQTLs_readDepth5)
# 
# # re-insert ASBs with at least minReadsPerAllele reads in ref and alt (based on the input with presence in at least 5 in 6 repeats)
# 
# # blacklist areas 
# minReadsPerAllele = minReadDepth * 0.1 
# postTotalSafe <- subset(postTotal, postTotal$refCount >= minReadsPerAllele & postTotal$altCount >= minReadsPerAllele)
# postTotalCheck <- subset(postTotal, !postTotal$id %in% postTotalSafe$id)
# postTotalCheck <- subset(postTotalCheck, postTotalCheck$id %in% df.bQTLsInput$id)
# postTotal <- rbind(postTotalSafe, postTotalCheck)
# 
# if(b.save_intermediate_results){
#   write.csv(postTotal, "paper_tmp_files/bQTL_after_input.csv")
# }
# 
# 
# # filtering
# n.bQTLs_readDepth50.both_alleles_10Percent_bias <- nrow(postTotal)
# n.bQTLs_readDepth5.both_alleles_nonZero_cound <- nrow(df.bQTLsInput)
# 
# print(n.bQTLs_readDepth50.both_alleles_10Percent_bias)
# print(n.bQTLs_readDepth5.both_alleles_nonZero_cound)
# 
# message("input bQTL after min reads per allele filter: ", nrow(postTotal))
# 
# v.postFreq <- postTotal$refCount / postTotal$totalCount
# v.postFreq <- v.postFreq[which(v.postFreq != 0)]
# v.postFreq <- v.postFreq[which(v.postFreq != 1)]
# prob.bias <- median(v.postFreq)
# 
# message("bias probability: ", prob.bias)
# 
# #p = 0.5 # 0.551948 # expectancy
# p = 1 - prob.bias # 0.5120976 # post frequency - ref / all
# test <- sapply(1:nrow(postTotal), function(i) binom.test(as.integer(postTotal$altCount[i]), as.integer(postTotal$totalCount[i]), p = p)$p.value) #  analytic.pv(preProps[i],preVars[i],postProps[i],depths[i]))
# fdr <- p.adjust(test, "bonferroni")
# postTotal["p-value (corrected)"] <- fdr
# 
# postTotal <- subset(postTotal, postTotal$totalCount >= minReadDepth)
# 
# postTotal["POSTfreq"] <- postTotal$refCount / postTotal$totalCount
# postTotal["POSTallele"] <- ifelse(postTotal$altCount > postTotal$refCount, postTotal$altAllele, postTotal$refAllele)
# 
# # postTotal <- subset(postTotal, postTotal$POSTfreq < 1 & postTotal$POSTfreq > 0)
# median(postTotal$POSTfreq)
# 
# message("make figure per paper")
# df.bQTLsInput["POSTfreq"] <- df.bQTLsInput$refCount / df.bQTLsInput$totalCount
# # df.bQTLsInput <- subset(df.bQTLsInput, df.bQTLsInput$totalCount > 5)
# 
# # perform chromosome based - blacklisting
# 
# message("input bQTL: ", nrow(postTotal))
# bQTL_scatterplot(postTotal=postTotal)
# 
# ## manual artifact removal
# message("remove manually labelled artifact regions - further reducing significant ")
# # make figure 
# # bQTL_scatterplot_chr(postTotal=postTotal, chr = 3)
# chr = 3
# cut.left <-  186240476
# cut.right <- 210080072
# print(abs(cut.left - cut.right) / 1e6)
# 
# idx.cut <- which(df.bQTLsInput$contig == chr & df.bQTLsInput$position < cut.right & df.bQTLsInput$position > cut.left)
# df.bQTLsInput <- df.bQTLsInput[!seq(1:nrow(df.bQTLsInput)) %in%  idx.cut,]
# 
# chr = 10
# cut.left <- 2965337
# cut.right <- 4124958
# print(abs(cut.left - cut.right) / 1e6)
# 
# idx.cut <- which(df.bQTLsInput$contig == chr & df.bQTLsInput$position < cut.right & df.bQTLsInput$position > cut.left)
# df.bQTLsInput <- df.bQTLsInput[!seq(1:nrow(df.bQTLsInput)) %in%  idx.cut,]
# 
# # double check
# # bQTL_scatterplot(postTotal=df.bQTLsInput, chr = 10)
# 
# postTotal <- subset(postTotal, postTotal$id %in% df.bQTLsInput$id)
# postTotal <- postTotal[,!names(postTotal) %in% "id"]
# 
# message("input bQTL after min reads per allele filter: ", nrow(postTotal))
# 
# if(b.save_intermediate_results){
#   write.csv(postTotal, "paper_tmp_files/bQTL_after_input_blacklisting.csv")
# }
# 
# 
# ## automatically run until here
# 
# 
# #df.bQTLsInput
# # background blacklisting and top significant filter
# 
# # significant too 
# 
# # install.packages("plotly")
# 
# 
# 
# # 15 x 5 pdf 
# # save current workspace
# # save.image("tmp/workspace_til_bQTL.RData")
# # load("tmp/workspace_til_bQTL.RData")
# 
# #df.transcripts <- read.fasta("Zmays_284_5b+.transcript_primaryTranscriptOnly.fa")
# #v.gns <- names(df.transcripts)
# 
# # laoding Z.mays annotation information, including gene identifiers and positions 
# # df.annotation <- read.table(path.annotation, header = FALSE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
# 
# #l.gffSets <- vector(mode = "list", length = 2)
# message("loading gene annotation datasets...")
# strt<-Sys.time()
# 
# df.gff <- read.table(path.gff, header = FALSE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE, quote = "")
# # 
# # df.gff["id"] <- sapply(df.gff$V9, function(m) {gsub("Name=", "", unlist(strsplit(m, ";"))[2])})
# # df.gff$id <- gsub("Parent=","", df.gff$id)
# # df.gff$id <- gsub(".v6a","", df.gff$id)
# 
# v.genePartitions <- c("gene", "five_prime_UTR", "CDS", "three_prime_UTR", "exon")
# df.gff <- subset(df.gff, df.gff$V3 %in% v.genePartitions)
# 
# 
# ids.list <- sapply(df.gff$V9, function(m) {strsplit(m, ";")})
# 
# 
# df.gff["id"] <- sapply(df.gff$V9, function(m) {gsub("ID=", "", unlist(strsplit(m, ";"))[1])})
# #df.gff <- subset(df.gff, df.gff$V3 == "gene")
# df.gff$id <- gsub("gene:","", df.gff$id)
# df.gff$id <- gsub("Parent=transcript:","", df.gff$id)
# df.gff$id <- gsub("CDS:","", df.gff$id)
# df.gff$id <- gsub("\\_.*","", df.gff$id)
# 
# df.gff["transcript.id"] <- sapply(df.gff$V9, function(m) {gsub("ID=", "", unlist(strsplit(m, ";"))[1])})
# #df.gff <- subset(df.gff, df.gff$V3 == "gene")
# df.gff$transcript.id <- gsub("gene:","", df.gff$transcript.id)
# df.gff$transcript.id <- gsub("Parent=transcript:","", df.gff$transcript.id)
# df.gff$transcript.id <- gsub("CDS:","", df.gff$id)
# df.gff$transcript.id <- gsub("\\_.*","", df.gff$id)
# 
# 
# df.gene_annotation <- df.gff
# names(df.gene_annotation) <- c("chr", "source", "partition", "pos.start", "pos.stop", "aux1", "strand", "aux2", "gene_meta", "gene.ID")
# 
# 
# 
# df.gene_function <- read.table(path.gene_function, header = FALSE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE, quote = "")
# names(df.gene_function) <- c("gene.ID", "gene.function")
# 
# df.geneID_conversion <- read.table(path.geneID_conversion, header = FALSE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE, quote = "")[,1:2]
# names(df.geneID_conversion) <- c("gene.ID.AGPv3", "gene.ID.AGPv4")
# 
# df.gene_conversion.AGPv3_to_AGPv4 <- df.geneID_conversion
# 
# print(Sys.time() - strt)
# 
# 
# 
# 
# 
# #l.gffSets[[1]] <- df.gff
# 
# # df.gff <- read.table("Zmays_284_5b+.gene.gff3", header = FALSE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
# # df.gff["id"] <- sapply(df.gff$V9, function(m) {gsub("ID=", "", unlist(strsplit(m, ";"))[1])})
# # df.gff <- subset(df.gff, df.gff$V3 == "CDS")
# # df.gff$id <- gsub(".v6a","", df.gff$id)
# # l.gffSets[[2]] <- df.gff
# # names(l.gffSets) <- c("in promoter", "not in gene")
# 
# # test <- subset(df.gff, df.gff$id %in% v.gns)
# message("loading GWAS datasets...")
# strt<-Sys.time()
# # df.report_gwas <- read.table(path.report_gwas, header = FALSE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
# # df.report_gwas <- subset(df.report_gwas, df.report_gwas$V18 == "Primary Assembly")
# # df.report_gwas <- df.report_gwas[,c(4,5,8,13)]
# #names(df.report_gwas) <- c("chr.old", "chr", "pos.old", "pos")
# # 
# df.phenotype_gwas.original <- read.table(path.phenotype_gwas.original, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
# # names(df.phenotype_gwas.original) <- c("trait", "V4", "V8", "allele", "rmip", "source")
# # df.phenotype_gwas.original <- subset(df.phenotype_gwas.original, df.phenotype_gwas.original$source == "Hapmap2")
# df.phenotype_gwas.original <- subset(df.phenotype_gwas.original, df.phenotype_gwas.original$rmip >= 5)
# df.phenotype_gwas <- df.phenotype_gwas.original# 
# v.traits  <- unique(df.phenotype_gwas$trait)
# 
# # df.phenotype_gwas <- read.table(path.phenotype_gwas, header = FALSE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
# # df.phenotype_gwas <- merge(df.phenotype_gwas, df.phenotype_gwas.original, by = c("V4","V8"))
# # 
# # df.phenotype_gwas <- subset(df.phenotype_gwas, df.phenotype_gwas$source == "Hapmap2")
# # df.phenotype_gwas <- subset(df.phenotype_gwas, df.phenotype_gwas$rmip >= 5)
# print(Sys.time() - strt)
# 
# 
# message("loading ChipSeq B73  datasets...")
# 
# l.df.ChipSeqB73 <- vector(mode = "list", length = length(l.path.ChipSeqB73))
# for(i in 1:length(l.path.ChipSeqB73)){
#   l.df.ChipSeqB73[[i]] <- read.table(l.path.ChipSeqB73[[i]], header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE, quote = "")  
#   l.df.ChipSeqB73[[i]]["chr"] <- as.numeric(gsub("\\:.*", "", l.df.ChipSeqB73[[i]]$Position))
#   l.df.ChipSeqB73[[i]]["pos"] <- as.numeric(gsub(".*:", "", l.df.ChipSeqB73[[i]]$Position))
#   l.df.ChipSeqB73[[i]]["left"] <- l.df.ChipSeqB73[[i]]$pos - 201
#   l.df.ChipSeqB73[[i]]["right"] <- l.df.ChipSeqB73[[i]]$pos + 201
# }
# 
# 
# if(!b.load_ChipSeqB73_filtered_peaks){
#   strt<-Sys.time() 
#   #cl<-makeCluster(min(n.chromosomes, n.cpus))
#   #registerDoParallel(cl)
#   # l.peaks.filtered <-  foreach(i = 1:n.chromosomes, .packages=c("seqinr", "VariantAnnotation", "Biostrings")) %dopar% { 
#   df.ChipSeq.filtered <- c()
#   for(i in 1:n.chromosomes){
#     l.peaks.chromosome <- vector(mode = "list", length = 3)
#     l.idx.blacklist <- vector(mode = "list", length = 3)
#     for(l in 1:length(l.df.ChipSeqB73)){
#       l.peaks.chromosome[[l]] <- subset(l.df.ChipSeqB73[[l]], l.df.ChipSeqB73[[l]]$chr == i)
#       l.idx.blacklist[[l]] <- numeric(nrow(l.peaks.chromosome[[l]]))
#     }
#     v.idx.keep <- numeric(nrow(l.peaks.chromosome[[1]]))
#     pb <- txtProgressBar(min = 0, max = nrow(l.peaks.chromosome[[1]]), style = 3)
#     for(j in 1:nrow(l.peaks.chromosome[[1]])){
#       setTxtProgressBar(pb, j)
#       n.found <- 0
#       for(l in 2:length(l.peaks.chromosome)){
#         b.found = FALSE
#         for(k in 1:nrow(l.peaks.chromosome[[l]])){
#           if(!k %in% l.idx.blacklist[[l]]){        
#             if(l.peaks.chromosome[[1]]$chr[j] == l.peaks.chromosome[[l]]$chr[k]){
#               if(l.peaks.chromosome[[l]]$pos[k] > l.peaks.chromosome[[1]]$left[j] & l.peaks.chromosome[[l]]$pos[k] < l.peaks.chromosome[[1]]$right[j]){
#                 n.found = n.found + 1
#                 b.found = TRUE
#                 l.idx.blacklist[[l]][k] <- k
#                 break
#               }
#             }
#           }
#         }
#         if(b.found == FALSE){
#           break
#         }
#       }
#       if(n.found == 2){
#         v.idx.keep[j] <- j  
#       }
#     }
#     close(pb)
#     v.idx.keep <- v.idx.keep[v.idx.keep != 0]
#     l.peaks.chromosome[[1]] <- l.peaks.chromosome[[1]][v.idx.keep,]
#     df.ChipSeq.filtered <- rbind(df.ChipSeq.filtered, l.peaks.chromosome[[1]])
#     # l.peaks.chromosome[[1]]
#   }
#   #stopCluster(cl)
#   print(Sys.time()-strt)
#   
#   saveRDS(df.ChipSeq.filtered, "tmp/df.ChipSeq.filtered.rds")
# }else{
#   df.ChipSeq.filtered <- readRDS("tmp/df.ChipSeq.filtered.rds")
#   df.ChipSeqB73 <- df.ChipSeq.filtered
#   write.table(df.ChipSeqB73, paste(folder_output, "/SX2.txt", sep = ""),  row.names = FALSE, quote = FALSE, sep ="\t")
# }
# 
# 
# 
# 
# 
# 
# 
# # 
# 
# 
# # eQTL
# message("loading eQTL datasets...")
# strt<-Sys.time()
# # df.eQTL_hotspots <- read.table(path.eQTL_hotspots, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
# # 
# # df.eqtl.region <- read.table(path.eQTL_regions, header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
# # df.eqtl.region$eQTL_region_start_AGPv2 <- df.eqtl.region$eQTL_region_start_AGPv2 - 1
# # df.eqtl.region$eQTL_region_start_AGPv3 <- df.eqtl.region$eQTL_region_start_AGPv3 - 1
# # df.eqtl.region <- subset(df.eqtl.region, df.eqtl.region$recip == "First Pass")
# # 
# # df.eqtl.leadSnp <- read.table(path.eQTL_leadSNP, header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
# # df.eqtl.leadSnp <- subset(df.eqtl.leadSnp, df.eqtl.leadSnp$recip == "First Pass")
# # #df.eqtl.leadSnp$Lead_SNP_pos_AGPv2 <- df.eqtl.leadSnp$Lead_SNP_pos_AGPv2 - 1
# # #df.eqtl.leadSnp$Lead_SNP_pos_AGPv3 <- df.eqtl.leadSnp$Lead_SNP_pos_AGPv3 - 1
# # df.eqtl.agpv3 <- read.table(path.eQTL_AGPV3, header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
# # names(df.eqtl.agpv3)[12:13] <- c("Lead_SNP_chr_AGPv3", "Lead_SNP_pos_AGPv3")
# # df.eqtl.agpv3["eQTL_region_chr_AGPv3"] <- NA
# 
# df.eqtls <- read.table(path.eQTLs, header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
# 
# # df.eqtls.liu.agp3 <- read.table(path.eQTLs_liu_AGP3, header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
# # df.eqtls.liu.agp2 <- read.table(path.eQTLs_liu_AGP2, header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
# # names(df.eqtls.liu.agp2)[1:4] <- c("feat_name", "source_id", "source_start", "source_stop")
# # df.eqtls.liu.agp3 <- merge(df.eqtls.liu.agp3, df.eqtls.liu.agp2, by = c("feat_name", "source_id", "source_start", "source_stop"), all = TRUE)
# # 
# # df.eqtls.liu.agp3 <- unique(df.eqtls.liu.agp3[,c(1,7,13,14)])
# # df.eqtls.liu.agp3 <- subset(df.eqtls.liu.agp3, !is.na(df.eqtls.liu.agp3$mapped_id))
# print(Sys.time() - strt)
# message("loading hypersensitivity datasets...")
# df.MNase <- read.table(path.MNase, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
# names(df.MNase) <- c("chr", "start", "end")
# print(Sys.time() - strt)
# 
# # 
# # df.enhancer_H3K9 <- read.table("datasets_paper/Enhancer_HM/GSE94251_H3K9ac_ist.bedGraph", header = FALSE, sep ="\t", quote = "", stringsAsFactors = FALSE)
# # df.enhancer_DNase <- read.table("datasets_paper/Enhancer_HM/GSE94291_DNase_ist.bedGraph", header = FALSE, sep ="\t", quote = "", stringsAsFactors = FALSE)
# # 
# # # MNase
# # df.dataset2 <- import.bw("datasets_paper/OneDrive-2018-02-21/MartBS123_CHG.bw")
# # df <- snpSummary(df.dataset)
# # 
# # 
# # df <- as.data.frame(df.dataset)
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # df.snp.pos <- rowRanges(df.dataset)
# # #save.image("tmp/workspace_til_hypersensitivity.RData")
# # 
# # # STANDARDIZED DATASETS - finish # 
# # message("standardizing all dataset formats - needs pre-configuration!...")
# # # standardized eqtl datasets # 
# # # df.eqtls <- df.eqtls.liu.agp3
# # # names(df.eqtls) <- c("targetGene", "chr", "pos.start", "pos.end")
# # 
# # # gene annotation
# # 
# # # rnaseq gene expression
# # names(df.rnaseq.down_regulated)[1] <- "gene.ID"
# # names(df.rnaseq.up_regulated)[1] <- "gene.ID"
# # 
# # #df.gene_expression <- df.rnaseq
# #names(df.gene_expression) <- c("gene.ID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
# 
# #df.gene_expression["diffExp"] <- ifelse(df.gene_expression$pvalue < 0.05, TRUE, FALSE)
# # df.gene_expression["mode"] <- ifelse(df.gene_expression$log2FoldChange >= 0.5 & df.gene_expression$pvalue < 0.05, "up", ifelse(df.gene_expression$log2FoldChange <= - 0.5 & df.gene_expression$pvalue < 0.05, "down", ""))
# 
# # GWAS datasets 
# # df.phenotype_gwas <- df.phenotype_gwas[,c("V5", "V13", "V14", "V15", "trait")]
# # names(df.phenotype_gwas) <- c("chr", "pos.start", "pos.end", "strand", "trait")
# 
# # paste(Sys.time(), sep = "")
# # format(Sys.time(), "%a %b %d %X %Y")
# 
# #Sys.Date()
# 
# # for(i in 1:nrow(df.eqtl.agpv3)){
# #   df.tmp <- subset(df.eqtl.region, df.eqtl.region$eQTL_region_chr_AGPv2 == df.eqtl.agpv3$eQTL_region_chr_AGPv2[i])
# #   df.tmp <- subset(df.tmp, df.tmp$eQTL_region_start_AGPv2 == df.eqtl.agpv3$eQTL_region_start_AGPv2[i] & df.tmp$eQTL_region_end_AGPv2 == df.eqtl.agpv3$eQTL_region_end_AGPv2[i])
# #   df.eqtl.agpv3$eQTL_region_start_AGPv3[i] <- df.tmp$eQTL_region_start_AGPv3[1]
# #   df.eqtl.agpv3$eQTL_region_end_AGPv3[i] <- df.tmp$eQTL_region_end_AGPv3[1]
# #   df.eqtl.agpv3$eQTL_region_chr_AGPv3[i] <- df.tmp$eQTL_region_chr_AGPv3[1]
# #   df.tmp <- subset(df.eqtl.leadSnp, df.eqtl.leadSnp$Lead_SNP_chr_AGPv2 == df.eqtl.agpv3$Lead_SNP_chr_AGPv2[i])
# #   df.tmp <- subset(df.tmp, df.tmp$Lead_SNP_pos_AGPv2 == df.eqtl.agpv3$Lead_SNP_pos_AGPv2[i])
# #   df.eqtl.agpv3$Lead_SNP_chr_AGPv3[i] <- df.tmp$Lead_SNP_chr_AGPv3[1]
# #   df.eqtl.agpv3$Lead_SNP_pos_AGPv3[i] <- df.tmp$Lead_SNP_pos_AGPv3[1]
# # }
# # 
# # df.eqtl.agpv3 <- subset(df.eqtl.agpv3, !is.na(df.eqtl.agpv3$eQTL_region_chr_AGPv3))
# # df.eqtl.agpv3 <- subset(df.eqtl.agpv3, !is.na(df.eqtl.agpv3$eQTL_region_start_AGPv3))
# # df.eqtl.agpv3 <- subset(df.eqtl.agpv3, !is.na(df.eqtl.agpv3$eQTL_region_end_AGPv3))
# # df.eqtl.agpv3 <- subset(df.eqtl.agpv3, !is.na(df.eqtl.agpv3$Lead_SNP_chr_AGPv3))
# # df.eqtl.agpv3 <- subset(df.eqtl.agpv3, !is.na(df.eqtl.agpv3$Lead_SNP_pos_AGPv3))
# # 
# # save.image("tmp/workspace_all_datasets.RData")
# 
# #df.tmp <- merge(df.eqtl.agpv3, df.eqtl.region, by = c("eQTL_region_chr_AGPv2", "eQTL_region_start_AGPv2", "eQTL_region_end_AGPv2"))
# 
# #df.phenotype_gwas <- read.table("/Shared/Everyone/Michael_Thomas/GWAS/Wallace_etal_2014_PLoSGenet_GWAS_hits-150112.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
# #df.phenotype_gwas <- subset(df.phenotype_gwas, df.phenotype_gwas$source == "Hapmap2")
# #df.phenotype_gwas <- subset(df.phenotype_gwas, df.phenotype_gwas$rmip >= 5)
# 
# # report_gwas <- "/Shared/Everyone/Michael_Thomas/GWAS/report_GWAS_hits_filtered.txt.xls"
# # 
# # df.report_gwas <- read.table(report_gwas, header = FALSE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
# # df.report_gwas <- subset(df.report_gwas, df.report_gwas$V18 == "Primary Assembly")
# # df.report_gwas <- df.report_gwas[,c(4,5,8,13)]
# # 
# # df.report_gwas <- subset(df.report_gwas, df.report_gwas$V8 %in% df.phenotype_gwas$pos)
# 
# # hypersensitivty
# 
# # # laoding Z.mays annotation information, including gene identifiers and positions 
# # df.annotation <- read.table("Zmays_284_5b+.annotation_info.txt", header = FALSE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
# # 
# # df.gff <- read.table("Zmays_284_5b+.gene.gff3", header = FALSE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
# # df.gff["id"] <- sapply(df.gff$V9, function(m) {gsub("ID=", "", unlist(strsplit(m, ";"))[1])})
# # #df.gff <- subset(df.gff, df.gff$V3 == "gene")
# # df.gff$id <- gsub(".v6a","", df.gff$id)
# # 
# # # test all chromosomes - normal SNPs (no deletions or insertions considered)
# # if(FALSE){
# #   for(i in 1:10){
# #     
# #     print(paste("processing chromosome ", i))
# #     
# #     # assemble the entire genome (all 12 chromosomes)
# #     vcf <- readVcf(paste("mo17_snptables/",v.mo17_files[i],sep=""),"BSgenome.Zmays.NCBI.AGPv3")
# #     
# #     if(FALSE){
# #       df <- snpSummary(vcf)
# #       idx.homozygous <- which(df$g11 == 1) # remove all non homozygous
# #       
# #       df.snp.pos <- rowRanges(vcf)
# #       df.snp.pos <- df.snp.pos[idx.homozygous,]
# #       
# #       saveRDS(df.snp.pos, paste("df.snp.pos_",i, ".rds"))
# #     }
# #     
# #     # new 11 20  
# #     vcf2 <- vcf[geno(vcf)$GT == "1/1" & ref(vcf) %in% c("A", "G", "T", "C")] # & alt(vcf) %in% c("A", "G", "T", "C")]
# #     vcf2 <- vcf2[unlist(lapply(alt(vcf2), function(m) {length(m) == 1}))]  #%in% c("A", "G", "T", "C")]
# #     vcf2 <- vcf2[unlist(lapply(alt(vcf2), function(m) {m %in% c("A", "G", "T", "C")}))]
# #     
# #     writeVcf(vcf2, paste("/Shared/Everyone/Michael_Thomas/recode_vcf/snps_homozygous_chr",i, ".vcf"))
# #     
# #   }
# # }
# # 
# # 
