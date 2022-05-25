message("Build mutated genome from input SNPs")

for(i in 1:n.chromosomes){
  
  print(paste("processing chromosome ", i))
  
  filename <- paste(path.snptables,v.mo17_files[i],sep="")
  
  # assemble the entire genome (all 12 chromosomes)
  vcf <- readVcf(filename,"BSgenome.Zmays.NCBI.AGPv4")
  
  df <- snpSummary(vcf)
  df.snp.pos <- rowRanges(vcf)
  
  if(b.doInitialHomozygousFilter){
    idx.homozygous <- which(df$g11 == 1) # remove all non homozygous
    df.snp.pos <- df.snp.pos[idx.homozygous,]
  }
  
  # subset to first alternative 
  saveRDS(df.snp.pos, paste("tmp/df.snp.pos_",i, ".rds"))
  
}

df.snp_positions <- c()
l.genome.mutant <- vector(mode = "list", length = n.chromosomes)
# welche haben laenge 
for(i in 1:n.chromosomes){
  
  print(paste("processing chromosome ", i))
  
  df.snp.pos <- readRDS(paste("tmp/df.snp.pos_",i, ".rds"))
  
  # n.snps <- n.snps + length(df.snp.pos)
  
  # remove heterozygote snps
  # idx.snps <- which(lapply(df.snp.pos@elementMetadata$ALT, length) == 1)
  
  # A - generic SNPs
  #vec.snp.pos <- as.numeric(df.snp.pos@ranges@start[idx.snps])
  vec.snp.pos <- as.numeric(df.snp.pos@ranges@start)
  # vec.snp.bases <- as.character(unlist(df.snp.pos@elementMetadata$ALT)) # 
  
  # vec.snp.bases <- as.character(unlist(lapply(df.snp.pos@elementMetadata$ALT, function(m) m[1])))
  vec.ref.bases <- as.character(df.snp.pos@elementMetadata$REF)
  
  vec.snp.bases <- CharacterList(df.snp.pos@elementMetadata$ALT)
  vec.snp.bases = unstrsplit(vec.snp.bases, sep = ",")
  vec.snp.bases <- gsub("\\,.*", "", vec.snp.bases)
  
  #vec.snp.bases <- as.character(unlist(df.snp.pos@elementMetadata$ALT[idx.snps]))
  
  vec.snp.bases <- ifelse(vec.snp.bases == "<DEL>", "N",vec.snp.bases) # both represented separately
  vec.snp.bases <- ifelse(vec.snp.bases == "<INS>", "N",vec.snp.bases)
  
  # remove empty strings # 
  idx.snp.exceptions <- which(vec.snp.bases == "")
  
  if(length(idx.snp.exceptions) > 0){
    vec.snp.pos <- vec.snp.pos[-idx.snp.exceptions]
    vec.snp.bases <- vec.snp.bases[-idx.snp.exceptions]
    vec.ref.bases <- vec.ref.bases[-idx.snp.exceptions]
  }
  
  df.snp_positions.i <- data.frame(ID = paste(i ,"_", vec.snp.pos, sep ="") , chromosome = rep(paste(i ,"_maternal", sep ="") , length(vec.snp.pos)) , 
                                   position = vec.snp.pos, strand =  rep(1, length(vec.snp.pos)),
                                   Ref_SNP = paste(vec.ref.bases ,"/", vec.snp.bases, sep =""))
  
  df.snp_positions <- rbind(df.snp_positions, df.snp_positions.i)
  
  genome.reference <- DNAString(genome[[i]])
  genome.mutant    <- replaceLetterAt(genome.reference, vec.snp.pos, vec.snp.bases) # vergleich vor austausch
  
  l.genome.mutant[[i]] <- genome.mutant
  
  #  v.sequences <- unlist(lapply(l.sequences, function(m) {m[[1]]}))
  # if(FALSE){
  #   Sequences = DNAStringSet(genome.mutant)
  #   #Sequences <- read.DNAStringSet(FastaFile, "fasta")
  #   names(Sequences) <- as.character(paste("chr",i,sep = ""))
  #   writeXStringSet(Sequences, paste("output/mo17Genome/chr",i,".fasta", sep = ""), format="fasta")
  # }
  
}

write.table(df.snp_positions,  paste("data/df.snp_positions.csv", sep = ""), col.names = FALSE, row.names = FALSE, sep ="\t")
saveRDS(l.genome.mutant, "tmp/genome_mutant.rds")
