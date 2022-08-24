# To determine potential enrichments, +/- 10 kbps surrounding enhancer regions were intersected with ASBs and bgSNPs.

snp_to_enhancer_distance <- function(df.bQTLs, 
                                     df.enhancer,
                                     n.binWidth = 10,
                                     n.maxDist = 10000){
  
  v.chromosomes <- unique(df.enhancer$chr)
  
  df.bQTLs_w_enhancer_annotation <- c() 
  for(c in v.chromosomes){
    
    df.enhancer.chr <- subset(df.enhancer, df.enhancer$chr == c)
    df.bQTLs.chr <- subset(df.bQTLs, df.bQTLs$`B73-chr` == c) 
    
    df.bQTLs.chr <- subset(df.bQTLs.chr, df.bQTLs.chr$non_genic == "yes")
    df.bQTLs.chr["in_enhancer_region"] <- F
    df.bQTLs.chr["dist_to_closest_enhancer_region"] <- NA # distance 
    
    for(j in 1:nrow(df.bQTLs.chr)){
      b.in_enhancer_region <- any(df.enhancer.chr$start < df.bQTLs.chr$`B73-pos`[j] & df.bQTLs.chr$`B73-pos`[j] < df.enhancer.chr$end)
      if(b.in_enhancer_region){
        df.bQTLs.chr$in_enhancer_region[j] = T
        df.bQTLs.chr$dist_to_closest_enhancer_region[j] <- 0
        
      }else{
        # identify the closest enhancer region to (before and behind)
        v.dist.abs <- c(abs(df.bQTLs.chr$`B73-pos`[j] - df.enhancer.chr$start), abs(df.bQTLs.chr$`B73-pos`[j] - df.enhancer.chr$end)) 
        v.dist.raw <- c((df.bQTLs.chr$`B73-pos`[j] -  df.enhancer.chr$start), (df.bQTLs.chr$`B73-pos`[j] - df.enhancer.chr$end))
        idx <- which(v.dist.abs == min(v.dist.abs))[1]
        df.bQTLs.chr$dist_to_closest_enhancer_region[j] <- v.dist.raw[idx]
      }
    }
    df.bQTLs_w_enhancer_annotation <- rbind(df.bQTLs_w_enhancer_annotation, df.bQTLs.chr)
  }

  # is not in enhancer region 
  v.dist_to_ASB <- df.bQTLs_w_enhancer_annotation$dist_to_closest_enhancer_region
  v.dist_to_ASB <- v.dist_to_ASB[v.dist_to_ASB > -n.maxDist & v.dist_to_ASB < n.maxDist] 
  v.breaks_sets <- seq(-n.maxDist, n.maxDist, n.binWidth)
  v.number_per_break <- numeric(length(v.breaks_sets))
  names(v.number_per_break) <- v.breaks_sets
  
  for(j in 1:(length(v.breaks_sets) - 1)){
    idx <- which(v.breaks_sets[j] <= v.dist_to_ASB & v.dist_to_ASB <= v.breaks_sets[j + 1])
    if(length(idx) > 0)
      v.number_per_break[j] <- length(idx)
  }
  v.number_per_break <- v.number_per_break[-length(v.number_per_break)]

  return(list(df.bQTLs_w_enhancer_annotation=df.bQTLs_w_enhancer_annotation, 
              v.dist_to_ASB=v.dist_to_ASB,
              v.number_per_break=v.number_per_break))
  
}


snp_to_enhancer_distance <- function(df.ASBs,
                                     df.bgSNPs,
                                     path.enhancer,
                                     df.enhancer_genes,
                                     n.binWidth = 10,
                                     n.maxDist = 100000 ){
  
  message("Identify ASB and bgSnps distance to enhancer regulation ")
  
  df.enhancer <- read.table(path.enhancer, header = FALSE, sep ="\t", quote = "", stringsAsFactors = FALSE)
  names(df.enhancer) <- c("chr", "start", "end")
  df.enhancer$chr <- paste("B73-chr", df.enhancer$chr, sep = "")
  
  print(table(df.enhancer$chr))

  res.ASBs.enhancer <- snp_to_enhancer_distance(df.ASBs, df.enhancer, n.binWidth, n.maxDist)
  res.bgSNPs.enhancer <- snp_to_enhancer_distance(df.bgSNPs, df.enhancer, n.binWidth, n.maxDist)
  
  df.distPlot <- rbind(data.frame(bin = as.numeric(names(res.ASBs.enhancer$v.number_per_break)), val = res.ASBs.enhancer$v.number_per_break, set = "ASBs"), 
                       data.frame(bin = as.numeric(names(res.bgSNPs.enhancer$v.number_per_break)), val = res.bgSNPs.enhancer$v.number_per_break, set = "bgSNPs"))
  df.distDensityPlot <- rbind(data.frame(val = res.ASBs.enhancer$v.dist_to_ASB, set = "ASBs"), 
                              data.frame(val = res.bgSNPs.enhancer$v.dist_to_ASB, set = "bgSNPs"))
  
  # Fig. 4. a) ASBs of ZmBZR1 are overrepresented at enhancer sites as compared to randomly selected background (bg) SNPs
  p.dist_density <- ggplot(df.distDensityPlot, aes(x=val, fill = set)) + geom_density(aes(group=set, colour=set), alpha = 0.2) + theme_bw() + scale_colour_manual(values = c("red", "black"))
  p.dist <- ggplot(df.distPlot, aes(x=bin, y=val,  fill = set, colour = set)) + stat_smooth(aes(x = bin, y = val), method = "lm",  formula = y ~ poly(x, 21), se = FALSE) + theme_bw() + scale_colour_manual(values = c("red", "black"))
  
  return(list(p.dist_density=p.dist_density, p.dist=p.dist)) 
  
}

# 
# snp_to_enhancer_distance <- function(l.bQTL_gene_partitioning,
#                                      df.enhancer_genes,
#                                      v.sets = c("ASB", "bgSNP"),
#                                      n.binWidth = 10,
#                                      n.maxDist = 100000 ){
#   
#   message("Identify ASB and bgSnps distance to enhancer regulation ")
#   
#   l.numbers <- vector(mode = "list", length = 2)
#   df.distPlot <- c()
#   df.distDensityPlot <- c()
#   
#   for(s in 1:2){
#     
#     postTotal.significant <- l.bQTL_gene_partitioning[[s]]
#     postTotal.significant_with_enhancer_annotation <- c()
#     
#     v.enhancer_regulated_genes <- c()
#     
#     for(i in 1:n.chromosomes){
#       
#       df.enhancer_genes.i <- subset(df.enhancer_genes, df.enhancer_genes$chr == i)
#       postTotal.significant.i <- subset(postTotal.significant, postTotal.significant$contig == i) 
#       postTotal.significant.i <- subset(postTotal.significant.i, postTotal.significant.i$non_genic == "yes")
#       postTotal.significant.i["is_in_enhancerRegion"] <- "no"
#       postTotal.significant.i["dist_to_closest_enhancer_region"] <- NA # distance 
#       
#       for(j in 1:nrow(postTotal.significant.i)){
#         b.ASB_in_enhancer_region <- any(df.enhancer_genes.i$start < postTotal.significant.i$position[j] & postTotal.significant.i$position[j] < df.enhancer_genes.i$end)
#         if(b.ASB_in_enhancer_region){
#           postTotal.significant.i$is_in_enhancerRegion[j] = "yes"
#           postTotal.significant.i$dist_to_closest_enhancer_region[j] <- 0
#           
#         }else{
#           # identify the closest enhancer region to (before and behind)
#           v.dist.abs <- c(abs(postTotal.significant.i$position[j] - df.enhancer_genes.i$start), abs(postTotal.significant.i$position[j] - df.enhancer_genes.i$end)) 
#           v.dist.raw <- c((postTotal.significant.i$position[j] -  df.enhancer_genes.i$start), (postTotal.significant.i$position[j] - df.enhancer_genes.i$end))
#           idx <- which(v.dist.abs == min(v.dist.abs))[1]
#           dist <- v.dist.raw[idx]
#           postTotal.significant.i$dist_to_closest_enhancer_region[j] <- dist
#         }
#       }
#       postTotal.significant_with_enhancer_annotation <- rbind(postTotal.significant_with_enhancer_annotation, postTotal.significant.i)
#     }
#     
#     l.numbers[[s]] <- table(postTotal.significant_with_enhancer_annotation$is_in_enhancerRegion)
#     sum(table(postTotal.significant_with_enhancer_annotation$is_in_enhancerRegion))
#     
#     # is not in enhancer region 
#     v.dist_to_ASB <- postTotal.significant_with_enhancer_annotation$dist_to_closest_enhancer_region
#     v.dist_to_ASB <- v.dist_to_ASB[v.dist_to_ASB > -n.maxDist & v.dist_to_ASB < n.maxDist] 
#     v.breaks_sets <- seq(-n.maxDist, n.maxDist, n.binWidth)
#     v.number_per_break <- numeric(length(v.breaks_sets))
#     names(v.number_per_break) <- v.breaks_sets
#     
#     for(j in 1:(length(v.breaks_sets) - 1)){
#       idx <- which(v.breaks_sets[j] <= v.dist_to_ASB & v.dist_to_ASB <= v.breaks_sets[j + 1])
#       if(length(idx) > 0)
#         v.number_per_break[j] <- length(idx)
#     }
#     
#     v.number_per_break <- v.number_per_break[-length(v.number_per_break)]
#     df.distPlot <- rbind(df.distPlot, data.frame(bin = as.numeric(names(v.number_per_break)), val = v.number_per_break, set = v.sets[s]))
#     df.distDensityPlot <- rbind(df.distDensityPlot, data.frame(val = v.dist_to_ASB, set = v.sets[s]))
#   }
#   
#   # comparison with bgSNPs
#   hitInSample <- l.numbers[[1]][2]
#   sampleSize <- l.numbers[[1]][1] + l.numbers[[1]][2]
#   hitInPop <- l.numbers[[2]][2]
#   popSize <- l.numbers[[2]][1] + l.numbers[[2]][2]
#   failInPop <- popSize - hitInPop
#   
#   # fc <- (l.numbers[[1]][2] / l.numbers[[1]][1]) / ( l.numbers[[2]][2] / l.numbers[[2]][1])
#   m <- matrix(c(l.numbers[[1]][2], l.numbers[[2]][2], l.numbers[[1]][1], l.numbers[[2]][1]), 2, 2)
#   dimnames(m) <- list(
#     group = c("ASBs", "bgSNPs"),
#     enhancer_in_max_distance = c("yes", "no")
#   )
#   res <- fisher.test(m, alternative='greater') # TODO: make correct test 
#   
#   message("Comparing ASBs and bgSNPs with respect to distance to nearest enhancer (max. distance is ", n.maxDist, " bps)")
#   print(m)
#   
#   message("Fold change of overrepresentation in ASBs vs bgSNPS is ", res$estimate , " with p-value of ", res$p.value, " (Fisher's exact test)")
#   
#   # Fig. 4. a) ASBs of ZmBZR1 are overrepresented at enhancer sites as compared to randomly selected background (bg) SNPs
#   p.dist_density <- ggplot(df.distDensityPlot, aes(x=val, fill = set)) + geom_density(aes(group=set, colour=set), alpha = 0.2) + theme_bw() + scale_colour_manual(values = c("red", "black"))
#   p.dist <- ggplot(df.distPlot, aes(x=bin, y=val,  fill = set, colour = set)) + stat_smooth(aes(x = bin, y = val), method = "lm",  formula = y ~ poly(x, 21), se = FALSE) + theme_bw() + scale_colour_manual(values = c("red", "black"))
#  
#   return(list(p.dist_density=p.dist_density, p.dist=p.dist)) 
#   
# }

# df.enhancer <- df.enhancer_H3K9
# names(df.enhancer) <- c("chr", "start",  "end", "val")
# 
# 
# # funktion to measure 
# message("non genic ASBs to nearest genes")
# df.distPlot <- c()
# #select gene - number 1
# df.gene_annotation.gene_only <- subset(df.gene_annotation, df.gene_annotation$partition == "gene")
# 
# # select specific ASB partition
# df.ASBs_nonGenic <- subset(l.bQTL_gene_partitioning[[1]], l.bQTL_gene_partitioning[[1]]$non_genic == "yes")
# # df.bgSNPs_nonGenic <- subset(l.bQTL_gene_partitioning[[2]], l.bQTL_gene_partitioning[[2]]$non_genic == "yes")
# 
# v.dist <- numeric()
# v.strands <- c("+", "-")
# for(i in 1:n.chromosomes){
#   df.data.i <- subset(df.ASBs_nonGenic, df.ASBs_nonGenic$contig == i)
#   df.gene_annotation.i <- subset(df.gene_annotation.gene_only, df.gene_annotation.gene_only$chr == i) 
#   for(t in 1:length(v.strands)){ # strand sensitive  
#     df.gene_annotation.i.j <- subset(df.gene_annotation.i, df.gene_annotation.i$strand == v.strands[t])
#     for(j in 1:nrow(df.data.i)){
#       if(v.strands[t] == "+"){      # before or after depends on the strand
#         diff_sign <- ((df.data.i$position[j] - df.gene_annotation.i.j$pos.start))
#         diff <- (abs(df.data.i$position[j] - df.gene_annotation.i.j$pos.start))
#       }else{
#         diff_sign <- ((df.gene_annotation.i.j$pos.stop))
#         diff <- (abs(df.data.i$position[j] - df.gene_annotation.i.j$pos.stop))
#       }
#       idx <- which(diff == min(diff))
#       v.dist <- c(v.dist, diff_sign[idx])  
#     }
#   }
# }
# 
# v.dist <- v.dist[v.dist > -10000 & v.dist < 10000]
# ggplot(data.frame(vals = v.dist), aes(vals)) + geom_density() + theme_bw() 
# ggplot(data.frame(vals = v.dist), aes(vals)) + geom_histogram(binwidth=100,fill="white",colour="black") + geom_density(aes(y=..count..*100)) + theme_bw() 



#####




# message("combined enhancer with gene expression from qTeller - up and down-regulation")
# 
# df.expression_enhancer_genes <- read.csv("data/Enhancer_HM/qTeller_enhancer.csv", header = TRUE, fill = TRUE, sep = ",")
# df.expression_enhancer_genes["r_apex_B73_vs_M17"] <- df.expression_enhancer_genes$B73_Apex / df.expression_enhancer_genes$Mo17_Apex
# df.expression_enhancer_genes["r_shoot_B73_vs_M17"] <- df.expression_enhancer_genes$B73_shoot / df.expression_enhancer_genes$Mo17_shoot
# 
# df.expression_enhancer_genes.B73 <- subset(df.expression_enhancer_genes, df.expression_enhancer_genes$r_apex_B73_vs_M17 > 2 | df.expression_enhancer_genes$r_shoot_B73_vs_M17 > 2)
# df.expression_enhancer_genes.Mo17 <-  subset(df.expression_enhancer_genes, df.expression_enhancer_genes$r_apex_B73_vs_M17 < 0.5 | df.expression_enhancer_genes$r_shoot_B73_vs_M17 < 0.5)
# nrow(df.expression_enhancer_genes)
# 
# (155 / 1425) / (3949 / 1346186) # TODI: resulting statistics 37 fold increase
# 
# 
# df.bQTL_RNAseq <- df.rnaseq.up_regulated[,c(1,3)]
# df.bQTL_RNAseq <- rbind(df.bQTL_RNAseq, df.rnaseq.down_regulated[,c(1,3)])
# 
# df.bQTL_RNAseq["mode"] <- "down"
# df.bQTL_RNAseq$mode[1:nrow(df.rnaseq.up_regulated)] <- "up"
# names(df.bQTL_RNAseq)[1] <- "gene.ID"
# 
# 
# subset(postTotal.significant, postTotal.significant$gene.ID %in% (intersect(df.expression_enhancer_genes.B73$gene_name, df.rnaseq.down_regulated$X)))
# 
# length(intersect(df.expression_enhancer_genes.Mo17$gene_name, df.rnaseq.up_regulated$X))
# #
# hitInSample <- 155
# sampleSize <- 1425
# hitInPop <- 3949
# popSize <- 1346186
# #
# failInPop <- popSize - hitInPop #(nrow(df.global.domains) - hitInPop)
# p.val <- (phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE))
# #
# 
# 
# fisher.test(matrix(c(hitInSample-1, hitInPop, failInPop, sampleSize), 2, 2), alternative='less')
# #     foldChange <- (hitInSample / sampleSize) / (hitInPop / popSize)
# 
# postTotal.significant_with_enhancer_annotation <- subset(postTotal.significant_with_enhancer_annotation, postTotal.significant_with_enhancer_annotation$is_in_enhancerRegion == "yes")
# 
# #postTotal.significant_with_enhancer_annotation <- subset(postTotal.significant_with_enhancer_annotation, postTotal.significant_with_enhancer_annotation$non_genic == "yes")
# #table(postTotal.significant_with_enhancer_annotation$is_in_enhancerRegion)
# 
# #####
# 
# postTotal.significant <- l.bQTL_gene_partitioning[[1]]
# postTotal.significant_with_enhancer_annotation <- c()
# 
# for(i in 1:n.chromosomes){
#   # df.enhancer_candidates.i <- subset(df.enhancer_candidates, df.enhancer_candidates$chr == i)
#   df.enhancer_candidates.i <- subset(df.enhancer_genes, df.enhancer_genes$chr == i)
#   postTotal.significant.i <- subset(postTotal.significant, postTotal.significant$contig == i) 
#   postTotal.significant.i["is_in_enhancerRegion"] <- "no"
#   postTotal.significant.i <- subset(postTotal.significant.i, postTotal.significant.i$non_genic == "yes")
#   for(j in 1:nrow(postTotal.significant.i)){
#     b.ASB_in_enhancer_region <- any(df.enhancer_candidates.i$start < postTotal.significant.i$position[j] & postTotal.significant.i$position[j] < df.enhancer_candidates.i$end)
#     if(b.ASB_in_enhancer_region){
#       postTotal.significant.i$is_in_enhancerRegion[j] = "yes"
#     }
#   }
#   postTotal.significant_with_enhancer_annotation <- rbind(postTotal.significant_with_enhancer_annotation, postTotal.significant.i)
# }
# 
# table(postTotal.significant_with_enhancer_annotation$is_in_enhancerRegion)
# 
# # work only on asb partitioning data (remove duplicates) 
# #postTotal.significant <- l.bQTL_gene_partitioning[[1]]
# postTotal.significant_with_enhancer_annotation <- c()
# for(i in 1:n.chromosomes){
#   df.enhancer_candidates.i <- subset(df.enhancer_candidates, df.enhancer_candidates$chr == i)
#   df.peaks.i <- subset(df.peaks, df.peaks$seqnames == i)
#   df.peaks.i["is_in_enhancerRegion"] <- "no"
#   for(j in 1:nrow(df.peaks.i)){
#     b.ASB_in_enhancer_region <- any(df.enhancer_candidates.i$start < df.peaks.i$posPeak[j] & df.peaks.i$posPeak[j] < df.enhancer_candidates.i$end)
#     if(b.ASB_in_enhancer_region){
#       df.peaks.i$is_in_enhancerRegion[j] = "yes"
#     }
#   }
#   postTotal.significant_with_enhancer_annotation <- rbind(postTotal.significant_with_enhancer_annotation, df.peaks.i)
# }
# 
# 
# table(postTotal.significant_with_enhancer_annotation$is_in_enhancerRegion)
# 
# 
# (149 / 1569) / (1555 / 28347)
# (155 / 1270) / (1555 / 28347)
# 
# 
# df.enhancer_genes_in_ASB <- c()
# for(i in 1:n.chromosomes){
#   df.enhancer_genes.i <- subset(df.enhancer_genes, df.enhancer_genes$chr == i)
#   df.enhancer_genes.i["contains_ASB"] <- "no"
#   df.enhancer_genes.i["ASB_pos"] <- NA
#   df.enhancer_genes.i["ASB_postfreq"] <- NA
#   postTotal.significant.i <- subset(postTotal.significant, postTotal.significant$contig == i) 
#   for(j in 1:nrow(df.enhancer_genes.i)){
#     b.ASB_in_enhancer_region <- any(df.enhancer_genes.i$start[j] < postTotal.significant.i$position & postTotal.significant.i$position < df.enhancer_genes.i$end[j])
#     if(b.ASB_in_enhancer_region){
#       idx.j <- which(postTotal.significant.i$position > df.enhancer_genes.i$start[j] & postTotal.significant.i$position < df.enhancer_genes.i$end[j])
#       df.enhancer_genes.i$contains_ASB[j] = "yes"
#       df.enhancer_genes.i$ASB_pos[j] <- paste(postTotal.significant.i$position[idx.j], collapse = ", ")
#       df.enhancer_genes.i$ASB_postfreq[j] <- paste(postTotal.significant.i$POSTfreq[idx.j], collapse = ", ")
#     }
#   }
#   df.enhancer_genes_in_ASB <- rbind(df.enhancer_genes_in_ASB, df.enhancer_genes.i)
# }
# df.enhancer_genes_in_ASB <- subset(df.enhancer_genes_in_ASB, df.enhancer_genes_in_ASB$contains_ASB == "yes")
# write.table(df.enhancer_genes_in_ASB, "output/df.enhancer_genes_in_ASB.txt", sep = "\t")
# 
# 
# #
# df.enhancer_genes_in_ASB <- read.table("output/df.enhancer_genes_in_ASB.txt", sep = "\t", stringsAsFactors = FALSE)
# 
# 
# df.enhancer_genes_in_ASB["Upstream_geneID_AGP3"] <- NA
# df.enhancer_genes_in_ASB["Downstream_geneID_AGP3"] <- NA
# 
# for(j in 1:nrow(df.enhancer_genes_in_ASB)){
#   
#   if(!is.na(df.enhancer_genes_in_ASB$Upstream_gene_geneID[j])){
#     idx <- which(df.geneID_conversion$gene.ID.AGPv4 == df.enhancer_genes_in_ASB$Upstream_gene_geneID[j])
#     if(length(idx) > 0){
#       df.enhancer_genes_in_ASB$Upstream_geneID_AGP3[j] <- df.geneID_conversion$gene.ID.AGPv3[idx[1]]  
#     }
#   }
#   
#   if(!is.na(df.enhancer_genes_in_ASB$Downstream_gene_Downstream_gene_geneID[j])){
#     idx <- which(df.geneID_conversion$gene.ID.AGPv4  == df.enhancer_genes_in_ASB$Downstream_gene_Downstream_gene_geneID[j])
#     if(length(idx) > 0){
#       df.enhancer_genes_in_ASB$Downstream_geneID_AGP3[j] <- df.geneID_conversion$gene.ID.AGPv3[idx[1]]  
#     }
#   }
#   
# }
# 
# write.table(df.enhancer_genes_in_ASB, "output/df.enhancer_genes_in_ASB_w_AGP3.txt", sep = "\t")
