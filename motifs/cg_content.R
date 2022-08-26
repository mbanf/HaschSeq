
# obsolete?
message("Computing nucleotide content surrounding binding peaks")

# mutual motif analysis in binding peaks
# bQTL 
# non bQTL 
# random selections (200 bp windows)


# analysis surrounding nucleotides of binding peak motifs 
windowSise <- 100
nucleotides <- c("A","G","C","T")
ecotypes <-  c("both", "mutant", "reference")

l.content <- vector(mode = "list", length = length(motifs))
names(l.content) <- motifs

for(k in 1:length(motifs)){
  l.content[[k]] <- vector(mode = "list", length = 3)
  names(l.content[[k]]) <- c("both", "mutant", "reference")
  
  for(j in 1:3){
    l.content[[k]][[j]] <- vector(mode = "list", length = 4)
    names(l.content[[k]][[j]]) <- nucleotides
    
    for(n in 1:4){
      l.content[[k]][[j]][[n]] <- as.data.frame(t(numeric(windowSise*2)))
    }
  }
}

# analyze motif (in binding peaks) surrounding nucleotide contents
for(i in 1:n.chromosomes){
  
  print(paste("processing chromosome ", i))
  
  df.results.binding <- read.table(paste("df.target_binding_results_genomewide_",i, ".xls", sep = ""), sep = "\t", header = TRUE, stringsAsFactors = FALSE) 
  df.snp.pos <- readRDS(paste("df.snp.pos_",i, ".rds"))
  
  # remove heterozygote snps
  idx.snps <- which(lapply(df.snp.pos@elementMetadata$ALT, length) == 1)
  
  # A - generic SNPs #
  vec.snp.pos <- as.numeric(df.snp.pos@ranges@start[idx.snps])
  vec.snp.bases <- as.character(unlist(df.snp.pos@elementMetadata$ALT[idx.snps]))
  
  vec.snp.bases <- ifelse(vec.snp.bases == "<DEL>", "N",vec.snp.bases) # both represented separately
  vec.snp.bases <- ifelse(vec.snp.bases == "<INS>", "N",vec.snp.bases)
  
  # remove empty strings # 
  idx.snp.exceptions <- which(vec.snp.bases == "")
  
  if(length(idx.snp.exceptions) > 0){
    vec.snp.pos <- vec.snp.pos[-idx.snp.exceptions]
    vec.snp.bases <- vec.snp.bases[-idx.snp.exceptions]
  }
  
  genome.reference <- DNAString(genome[[i]])
  genome.mutant    <- replaceLetterAt(genome.reference, vec.snp.pos, vec.snp.bases) # vergleich vor austausch
  
  ecotypes <- unique(df.results.binding$ecotype)
  
  for(m in 1:length(motifs)){
    
    print(paste(m, "of", length(motifs)))  
    
    df.results.binding.k <- subset(df.results.binding, df.results.binding$motif == motifs[m])
    
    for(s in 1:length(ecotypes)){
      
      df.results.binding.ks <- subset(df.results.binding.k, df.results.binding.k$ecotype == ecotypes[s])
      
      if(s %in% c(1,3)){
        genome.tmp <- genome.reference  
      }else{
        genome.tmp <- genome.mutant
      }
      
      for(j in 1:nrow(df.results.binding.ks)){
        
        # print(paste(j, "of", nrow(df.results.binding.ks)))  
        
        motif.start <- df.results.binding.ks$motif_start_position[j]
        
        for(n in 1:4){
          v.content <- numeric(windowSise*2)
          
          v.seq.reference.left <- subseq(genome.tmp, start=motif.start - windowSise, end=motif.start - 1)  
          v.seq.reference.right <- subseq(genome.tmp, start=motif.start + 6, end=motif.start + 5 + windowSise)  
          
          v.content[1:windowSise] <- ifelse(unlist(strsplit(as.character(v.seq.reference.left), "")) == nucleotides[n], 1, 0)
          v.content[(windowSise+1):(2*windowSise)] <- ifelse(unlist(strsplit(as.character(v.seq.reference.right), "")) == nucleotides[n], 1, 0)
          
          l.content[[m]][[s]][[n]] <- rbind(l.content[[m]][[s]][[n]], as.data.frame(t(v.content)))
          
        }
        
      }
      
    }
    
  }
  
}

saveRDS(l.content, "l.content.rds")


l.content <- readRDS("l.content.rds")

# print motif surroundings 
for(k in 1:length(motifs)){
  
  print(motifs[k])
  
  for(j in 1:3){
    print(ecotypes[j])
    
    cols <- c("black", "red", "blue", "green")
    plot(colSums(l.content[[k]][[j]][[1]])/nrow(l.content[[k]][[j]][[1]]), type = "l", ylim = c(0,0.5), xlim = c(85,115), main = paste(motifs[k], "(", ecotypes[j], ")", sep =""))
    
    for(n in 2:4){
      lines(colSums(l.content[[k]][[j]][[n]])/nrow(l.content[[k]][[j]][[n]]), col = cols[n])
    }
    
    legend('topright', nucleotides, lty=1, col=cols, bty='n', cex=.75)
    
  }
  
}





#### analyze random nucleotide contents
l.randomBackground <- vector(mode = "list", length =4)
for(n in 1:4)
  l.randomBackground[[n]] <- matrix(0, 10000, 200)

for(i in 1:n.chromosomes){
  
  print(paste("processing chromosome ", i))
  
  df.results.binding <- read.table(paste("df.target_binding_results_genomewide_",i, ".xls", sep = ""), sep = "\t", header = TRUE, stringsAsFactors = FALSE) 
  df.snp.pos <- readRDS(paste("df.snp.pos_",i, ".rds"))
  
  # remove heterozygote snps
  idx.snps <- which(lapply(df.snp.pos@elementMetadata$ALT, length) == 1)
  
  # A - generic SNPs #
  vec.snp.pos <- as.numeric(df.snp.pos@ranges@start[idx.snps])
  vec.snp.bases <- as.character(unlist(df.snp.pos@elementMetadata$ALT[idx.snps]))
  
  vec.snp.bases <- ifelse(vec.snp.bases == "<DEL>", "N",vec.snp.bases) # both represented separately
  vec.snp.bases <- ifelse(vec.snp.bases == "<INS>", "N",vec.snp.bases)
  
  # remove empty strings # 
  idx.snp.exceptions <- which(vec.snp.bases == "")
  
  if(length(idx.snp.exceptions) > 0){
    vec.snp.pos <- vec.snp.pos[-idx.snp.exceptions]
    vec.snp.bases <- vec.snp.bases[-idx.snp.exceptions]
  }
  
  genome.reference <- DNAString(genome[[i]])
  
  # genome.mutant    <- replaceLetterAt(genome.reference, vec.snp.pos, vec.snp.bases) # vergleich vor austausch
  
  i.start <- (i - 1) * 1000
  
  for(n in 1:4){
    
    i.samples <- sample((length(genome.reference) - windowSise*2), 1000)
    v.content <- numeric(windowSise*2)
    
    for(j in 1:1000){
      
      v.seq.reference.random <- subseq(genome.reference, start=i.samples[j], end=i.samples[j] + windowSise*2 - 1)  
      l.randomBackground[[n]][i.start + j,] <- ifelse(unlist(strsplit(as.character(v.seq.reference.random), "")) == nucleotides[n], 1, 0)
      
    }
    
  }
  
}


# plot random background 
cols <- c("black", "red", "blue", "green")
plot(colSums(l.randomBackground[[1]])/nrow(l.randomBackground[[1]]), type = "l", ylim = c(0,0.5), xlim = c(1,200))
for(n in 2:4){
  lines(colSums(l.randomBackground[[n]])/nrow(l.randomBackground[[n]]), col = cols[n])
}

####






###


l.content <- vector(mode = "list", length = length(motifs))
names(l.content) <- motifs

for(k in 1:length(motifs)){
  l.content[[k]] <- vector(mode = "list", length = 3)
  names(l.content[[k]]) <- c("mutual", "reference", "mutant")
  
  for(j in 1:3){
    l.content[[k]][[j]] <- vector(mode = "list", length = 4)
    names(l.content[[k]][[j]]) <- nucleotides
  }
}

df.results.binding <- data.frame(motif <- character(), chromosome = character(), start = numeric(), end = numeric(), motif_start_position =numeric(), ecotype = numeric(), regulation = character())
names(df.results.binding) <- c("motif",  "chromosome", "start", "end", "motif_start_position", "ecotype", "regulation")

print(paste("processing chromosome ", i))

tf_target_bind.sset <- subset(df.peaks, df.peaks$seqnames == i)
print(paste("processing chromosome ", i))

for(j in 1:length(motifs)){
  
  print(paste("processing motif ", j))
  
  vec.intersect <- read.table(paste("cisvarAnalysis/positions/df.",vec.chroms.snp[i],"_",motifs[j],"_intersect.xls", sep = ""), header = FALSE, sep = "\t")
  vec.col.unique <- read.table(paste("cisvarAnalysis/positions/df.",vec.chroms.snp[i],"_",motifs[j],"_reference.xls", sep = ""), header = FALSE, sep = "\t")
  vec.cvi.unique <- read.table(paste("cisvarAnalysis/positions/df.",vec.chroms.snp[i],"_",motifs[j],"_mutant.xls", sep = ""), header = FALSE, sep = "\t")
  
  #vec.intersect <- read.table(paste("positions/df.",vec.chroms.snp[i],"_",motifs[j],"_intersect.xls", sep = ""), header = FALSE, sep = "\t")
  #vec.col.unique <- read.table(paste("positions/df.",vec.chroms.snp[i],"_",motifs[j],"_reference.xls", sep = ""), header = FALSE, sep = "\t")
  #vec.cvi.unique <- read.table(paste("positions/df.",vec.chroms.snp[i],"_",motifs[j],"_mutant.xls", sep = ""), header = FALSE, sep = "\t")
  
  for(m in 1:nrow(tf_target_bind.sset)){
    print(paste(m, "of", nrow(tf_target_bind.sset)))     
    
    #             for(l in 1:nrow(vec.intersect)){
    #               motif.start <- vec.intersect$V1[l]
    #               if(tf_target_bind.sset$start[m] <= motif.start && motif.start <= tf_target_bind.sset$end[m]){
    #                 newrow <- data.frame(motif <- motifs[j], chromosome = tf_target_bind.sset$seqnames[m], 
    #                                      start = tf_target_bind.sset$start[m], end = tf_target_bind.sset$end[m],
    #                                      motif_start_position = motif.start, ecotype = "both")
    #                 names(newrow) <- c("motif", "chromosome", "start", "end", "motif_start_position", "ecotype")
    #                 df.results.binding <- rbind(df.results.binding, newrow)
    #               }
    #             }
    
    
    for(l in 1:nrow(vec.col.unique)){
      motif.start <- vec.col.unique$V1[l]
      
      if(tf_target_bind.sset$start[m] <= motif.start && motif.start <= tf_target_bind.sset$end[m]){ # motif in peak
        
        newrow <- data.frame(motif <- motifs[j], chromosome = tf_target_bind.sset$seqnames[m], 
                             start = tf_target_bind.sset$start[m], end = tf_target_bind.sset$end[m],
                             motif_start_position = motif.start, ecotype = "reference")
        names(newrow) <- c("motif", "chromosome", "start", "end", "motif_start_position", "ecotype")
        
        df.results.binding <- rbind(df.results.binding, newrow)
        
        for(n in 1:4){
          v.content <- numeric(200)
          for(k in 1:100){
            v.seq.reference.left <- subseq(genome.reference, start=motif.start - k, end=motif.start - k)  
            v.seq.reference.right <- subseq(genome.reference, start=motif.start + 5 + k, end=motif.start + 5 + k)  
            v.content[k] <- ifelse(unlist(strsplit(as.character(v.seq.reference.left), "")) == nucleotides[n], 1, 0)
            v.content[k+100] <- ifelse(unlist(strsplit(as.character(v.seq.reference.right), "")) == nucleotides[n], 1, 0)
          }
          l.content[[j]][[2]][[n]] <- rbind(l.content[[j]][[2]][[n]], as.data.frame(t(v.content)))
        }
        
      }
    }
    
    
    #                 v.surrounding.start <- c(1,10,20,30,40,50,60,70,80,90)
    #                 v.surrounding.end <- c(10,20,30,40,50,60,70,80,90,100)
    #                 
    #                 for(k in 1:length(v.surrounding.start)){
    #                   
    #                   v.seq.reference.left <- subseq(genome.reference, start=motif.start - v.surrounding.end[k], end=motif.start - v.surrounding.start[k])  
    #                   v.seq.reference.right <- subseq(genome.reference, start=motif.start + 5 + v.surrounding.start[k], end=motif.start + 5 + v.surrounding.end[k])  
    #                   
    #                   tb.percentage <- table(c(unlist(strsplit(as.character(v.seq.reference.left), "")), 
    #                                            
    #                                            
    #                                            
    #                                            unlist(strsplit(as.character(v.seq.reference.right), "")))) / 20
    #                 }
    #                   
    #                 
    #                 
    #                 
    
    
    
    
    for(l in 1:nrow(vec.cvi.unique)){
      motif.start <- vec.cvi.unique$V1[l]
      if(tf_target_bind.sset$start[m] <= motif.start && motif.start <= tf_target_bind.sset$end[m]){
        newrow <- data.frame(motif <- motifs[j], chromosome = tf_target_bind.sset$seqnames[m], 
                             start = tf_target_bind.sset$start[m], end = tf_target_bind.sset$end[m],
                             motif_start_position = motif.start, ecotype = "mutant")
        names(newrow) <- c("motif", "chromosome", "start", "end", "motif_start_position", "ecotype")
        df.results.binding <- rbind(df.results.binding, newrow)
        
        ###
        
        for(n in 1:4){
          v.content <- numeric(200)
          for(k in 1:100){
            v.seq.reference.left <- subseq(genome.reference, start=motif.start - k, end=motif.start - k)  
            v.seq.reference.right <- subseq(genome.reference, start=motif.start + 5 + k, end=motif.start + 5 + k)  
            v.content[k] <- ifelse(unlist(strsplit(as.character(v.seq.reference.left), "")) == nucleotides[n], 1, 0)
            v.content[k+100] <- ifelse(unlist(strsplit(as.character(v.seq.reference.right), "")) == nucleotides[n], 1, 0)
          }
          l.content[[j]][[3]][[n]] <- rbind(l.content[[j]][[3]][[n]], as.data.frame(t(v.content)))
        }
        
      }
    }
    
    
    
  }  
}

write.table(df.results.binding, paste("df.target_binding_results_genomewide_",i, ".xls", sep = ""), sep = "\t", row.names = FALSE) 

}





for(m in 1:length(motifs)){
  
  
  # motifs - binding peak 
  
  # analyze the nucleotide content around motif 
  
  
  for(j in 1:length(vec.reference.unique)){
    
    
    subseq(genome.reference, start=56730620, end =  56730633)
    subseq(genome.mutant, start=56730625, end =  56730633)
    
    v.seq.mutant <- subseq(genome.mutant, start=vec.mutant.unique[j], end=vec.mutant.unique[j] + 5)  
    v.seq.reference <- subseq(genome.reference, start=vec.mutant.unique[j], end=vec.mutant.unique[j] + 5)  
    
  }
  
}



#                 v.surrounding.start <- c(1,10,20,30,40,50,60,70,80,90)
#                 v.surrounding.end <- c(10,20,30,40,50,60,70,80,90,100)
#                 
#                 for(k in 1:length(v.surrounding.start)){
#                   
#                   v.seq.reference.left <- subseq(genome.reference, start=motif.start - v.surrounding.end[k], end=motif.start - v.surrounding.start[k])  
#                   v.seq.reference.right <- subseq(genome.reference, start=motif.start + 5 + v.surrounding.start[k], end=motif.start + 5 + v.surrounding.end[k])  
#                   
#                   tb.percentage <- table(c(unlist(strsplit(as.character(v.seq.reference.left), "")), 
#                                            
#                                            
#                                            
#                                            unlist(strsplit(as.character(v.seq.reference.right), "")))) / 20
#                 }
#                   
#                 
#                 
#                 