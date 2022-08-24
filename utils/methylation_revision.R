scatter_plot <- function(x, y, xlab = "allelic bias", ylab = "dnase level", x_limit = NULL, y_limit = NULL, col_vals = NULL){
  
  df <- data.frame(x,y)
  
  ## Use densCols() output to get density at each point
  v <- densCols(x,y, colramp=colorRampPalette(c("black", "white")))
  df$dens <- col2rgb(v)[1,] + 1L
  
  ## Map densities to colors
  cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", 
                              "#FCFF00", "#FF9400", "#FF3100"))(256)
  
  if(is.null(col_vals)){
    col_vals <- cols[df$dens]
  }
  df$col <- col_vals

  ## Plot it, reordering rows so that densest points are plotted on top
  plot(y~x, data=df[order(df$dens),], pch=20, col=col, cex=0.7, xlab = xlab, ylab = ylab, xlim=x_limit, ylim=y_limit)
}


get_asb_methylation_distance <- function(df.bQTLs.parent, 
                                         df.methylation,
                                         chrs,
                                         pos,
                                         chr.parent,
                                         n.bin_width = 10,
                                         th.distance_to_ASB = 2000,
                                         group = "all",
                                         n.chromosomes = 10,
                                         n.snps.split = 50){
  
  if(group == "genic"){
    df.bQTLs.parent <- subset(df.bQTLs.parent, df.bQTLs.parent$non_genic == "no")
  }else if(group == "non_genic"){
    df.bQTLs.parent <- subset(df.bQTLs.parent, df.bQTLs.parent$non_genic == "yes") 
  }
  
  v.bins <- seq(-th.distance_to_ASB, th.distance_to_ASB, n.bin_width)
  v.hist <- numeric(length(v.bins))
  names(v.hist) <- v.bins
  v.counter <- numeric(length(v.bins))
  
  for(i in 1:n.chromosomes){
    
    df.bQTLs.i <- subset(df.bQTLs.parent, df.bQTLs.parent[,chr.parent] == chrs[i])
    v.pos.bQTLs.i <- df.bQTLs.i[,pos]
    
    df.methylation.i <- subset(df.methylation, df.methylation$chr == chrs[i])
    df.methylation.i <- df.methylation.i[order(df.methylation.i$start),]
    df.methylation.i <- subset(df.methylation.i, df.methylation.i$start >= min(df.bQTLs.i[,pos] - th.distance_to_ASB) 
                               & df.methylation.i$start <= max(df.bQTLs.i[,pos]) + th.distance_to_ASB)
    
    v.pos.methylation <- df.methylation.i$start
    v.vals.methylation <- df.methylation.i$val 
    
    p.start <- v.pos.bQTLs.i[1] - th.distance_to_ASB
    p.end <- v.pos.bQTLs.i[min(length(v.pos.bQTLs.i), n.snps.split)] + th.distance_to_ASB
    
    idx.split <- which(p.start <= v.pos.methylation & v.pos.methylation <= p.end)
    v.pos.methylation.split <- v.pos.methylation[idx.split]
    v.vals.methylation.split <- v.vals.methylation[idx.split]
    
    for(j in 1:length(v.pos.bQTLs.i)){
      
      if(j %% n.snps.split == 0){
        p.start <- v.pos.bQTLs.i[j] - th.distance_to_ASB
        p.end <- v.pos.bQTLs.i[min(length(v.pos.bQTLs.i), (j + n.snps.split))] + th.distance_to_ASB
        idx.split <- which(p.start <= v.pos.methylation & v.pos.methylation <= p.end)
        v.pos.methylation.split <- v.pos.methylation[idx.split]
        v.vals.methylation.split <- v.vals.methylation[idx.split]
      }  
      
      idx <- which(v.pos.methylation.split >= v.pos.bQTLs.i[j] - th.distance_to_ASB & v.pos.methylation.split <= v.pos.bQTLs.i[j] + th.distance_to_ASB)
      if(length(idx) == 0)
        next 
      
      v.dist_to_asb <- v.pos.methylation.split[idx] - v.pos.bQTLs.i[j]
      v.vals <- v.vals.methylation.split[idx]
      
      for(k in 1:(length(v.bins) - 1)){
        idx <- which(v.bins[k] <= v.dist_to_asb & v.dist_to_asb <= v.bins[k + 1])
        if(length(idx) == 0)
          next
        v.hist[k] <- v.hist[k] + mean(v.vals[idx])
        v.counter[k] <- v.counter[k] + 1
      }
      
    }
    
  }
  
  v.hist <- v.hist / v.counter
  
  v.hist[is.infinite(v.hist)] <- 0
  
  v.hist
  
}


asb_methylation_frequency_at_distance <- function(df.bQTLs, 
                                                  df.methylation,
                                                  n.bin_width = 10,
                                                  th.distance_to_ASB = 2000,
                                                  th.post_frequency = 0.85,
                                                  parent = "B73",
                                                  group = "all",
                                                  n.chromosomes = 10,
                                                  n.snps.split = 50){
  if(parent == "B73"){ 
    chrs = paste("B73-chr",1:10, sep = "")
    pos = "B73-pos"
    chr.parent = "B73-chr"
    
  }else{
    chrs = paste("Mo17-chr",1:10, sep = "")
    pos = "Mo17-pos"
    chr.parent = "Mo17-chr"
  }
  
  
  df.bQTLs.parent <- subset(df.bQTLs, df.bQTLs$POSTfreq > th.post_frequency)
  
  v.hist.B73 <- get_asb_methylation_distance(df.bQTLs.parent, 
                                             df.methylation,
                                             chrs,
                                             pos,
                                             chr.parent,
                                             n.bin_width,
                                             th.distance_to_ASB,
                                             group,
                                             n.chromosomes,
                                             n.snps.split)
  
  
  df.bQTLs.parent <- subset(df.bQTLs, df.bQTLs$POSTfreq < 1 - th.post_frequency)
  
  v.hist.Mo17 <- get_asb_methylation_distance(df.bQTLs.parent, 
                                              df.methylation,
                                              chrs,
                                              pos,
                                              chr.parent,
                                              n.bin_width,
                                              th.distance_to_ASB,
                                              group,
                                              n.chromosomes,
                                              n.snps.split)

  if(parent == "B73"){
  
    return(list(v.hist.high_affinity=v.hist.B73,
                v.hist.low_affinity=v.hist.Mo17))
    
  }else{
    
    return(list(v.hist.high_affinity=v.hist.Mo17,
                v.hist.low_affinity=v.hist.B73))
    
  }

}




#' A result statistics
#'
#' This function analysis the results
#' @param df.cluster_annotations inferred cluster condition annotations 
#' @param m.functionality gene cluster acitivity heatmap
#' @param th.bp_offset (default = 20)
#' @keywords 
#' @export
#' @examples
#' cat_function()
methylation_occupancy_distance_to_asb  = function(df.ASBs, 
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
                                                 folder = "tmp"){
  
  message("(differential) methylation values versus post frequency")
  
  for(k in 1:length(v.filenames)){
    gc()
    message("Processing, ", v.filenames[k])
    chrs <- paste(v.parents[k], "-chr", 1:n.chromosomes, sep = "")
    
    if(v.formats[k] == "bedGraph"){
      df.methylation <- read.table(v.filenames[k], header = FALSE, sep ="\t", quote = "", stringsAsFactors = FALSE)
    }else if(v.formats[k] == "bigWig"){
      df.methylation <- import.bw(v.filenames[k])
      df.methylation <- as.data.frame(df.methylation)
      df.methylation <- df.methylation[,c(1,2,6)]
    }
    names(df.methylation) <- c("chr",  "start", "val")
    df.methylation <- subset(df.methylation, df.methylation$chr %in% chrs)
    df.methylation = subset(df.methylation, df.methylation$val > 0) 
    gc()
    
    message("Checking ", nrow(df.methylation), " entries of methylation for ", nrow(df.ASBs), " ASBs")

    res <- asb_methylation_frequency_at_distance(df.ASBs, 
                                                df.methylation,
                                                n.bin_width = n.bin_width,
                                                th.distance_to_ASB = th.distance_to_ASB,
                                                th.post_frequency = th.high_affinity,
                                                parent = v.parents[k],
                                                group = group,
                                                n.chromosomes = n.chromosomes,
                                                n.snps.split = 50)

    gc()
    saveRDS(res$v.hist.high_affinity, paste(folder, paste("v.hist.methylation.ASBs.high_affinity.",v.parents[k],".", v.variants[k],".",th.high_affinity,".rds", sep = ""), sep ="/"))  
    saveRDS(res$v.hist.low_affinity, paste(folder, paste("v.hist.methylation.ASBs.low_affinity.",v.parents[k],".", v.variants[k],".",th.high_affinity,".rds", sep = ""), sep ="/"))  
    
    rm(df.methylation)
    gc()
  }
  
}


plot_methylation_occupancy_distance_to_asb  = function(v.variant_sets = c("CG", "CHG", "CHH"),
                                                  th.high_affinity = 0.85,
                                                  folder = "tmp"){
  
  v.affinity = c("low affinity", "high affinity", "low affinity", "high affinity")
  v.parent_sets = c("B73", "B73", "Mo17", "Mo17")
  
  for(i in 1:length(v.variant_sets)){
    
    v.hist.methylation <- c(readRDS(paste(folder, paste("v.hist.methylation.ASBs.low_affinity.B73.", v.variant_sets[i],".",th.high_affinity, ".rds", sep = ""), sep ="/")),
                            readRDS(paste(folder, paste("v.hist.methylation.ASBs.high_affinity.B73.", v.variant_sets[i],".",th.high_affinity,".rds", sep = ""), sep ="/")),
                            readRDS(paste(folder, paste("v.hist.methylation.ASBs.low_affinity.Mo17.", v.variant_sets[i],".",th.high_affinity, ".rds", sep = ""), sep ="/")),
                            readRDS(paste(folder, paste("v.hist.methylation.ASBs.high_affinity.Mo17.", v.variant_sets[i],".",th.high_affinity,".rds", sep = ""), sep ="/")))
    
    n = length(v.hist.methylation)/4
    
    df.hist <- data.frame(val = as.numeric(v.hist.methylation), 
                          bin = as.numeric(names(v.hist.methylation)), 
                          color = rep(v.parent_sets, each = n),
                          line_type = rep(v.affinity, each = n))
    
    # ggplot(df.hist, aes(x=bin, y=val, colour = color, linetype = line_type)) + stat_smooth(aes(x = bin, y = val), method = "lm",  formula = y ~ poly(x, 21), se = FALSE, n = 1000) + scale_colour_manual(values = c("green", "blue")) + theme_bw() + xlim(-950, 950)
    p <- ggplot(df.hist, aes(x=bin, y=val, colour = color, linetype = line_type)) + geom_line() + scale_colour_manual(values = c("green", "blue")) + theme_bw() + xlim(-1990, 1990)
    # p <- p + ggtitle(v.variant_sets[i])
    plot(p)
    #ggplot(df.hist, aes(x=bin, y=val, colour = color, linetype = line_type)) + scale_colour_manual(values = c("green", "blue")) + theme_bw() + xlim(-2000, 2000)
    # pdf(paste(folder_output, paste("Novel_Methylation/", names(l.res_plots)[i],".pdf", sep = ""), sep ="/"), width = 6, height = 4)
    
    file = paste(folder, paste("df.hist.methylation_distance_to_ASBs.", v.variant_sets[i], ".csv", sep = ""), sep ="/")
    write.csv(df.hist, file, quote = F, row.names = F)

  }
  

  # mtext(title, side = 3, line = -2, outer = TRUE, cex = 0.8)
  # 
  # 
  # df.distSet <- rbind(df.distSet, data.frame(dist = as.numeric(v.dist_to_ASB), 
  #                                            val = as.numeric(v.vals.dist_to_ASB), 
  #                                            color = v.species[s], 
  #                                            line_type = v.stringency[s]))
  # 
  # 
  # for(k in 1:length(l.res_plots)){
  #   write.table(l.res_plots[[k]], paste(folder_output, paste("df.methylation_frequency_vs_distance_to_asb_binned_",v.datasets[k], "_", group,"_", n.binWidth, "_", timeStamp, ".txt", sep = ""), sep ="/"), sep = "\t", row.names = FALSE)
  # }
  # 
  # # 6 x 4
  # if(do.plot){
  #   # Fig. 4 i)-k) Average CpG (h), CHG (i) and CHH(j) methylation frequency in B73 (green) and Mo17 (blue) over ASB loci with a least 85% binding
  #   # bias towards B73 or Mo17. High affinity ( ___ ) and low affinity ( …. ) alleles are considered separately for each genotype.
  #   for(i in 1:length(l.res_plots)){
  #     p <- ggplot(l.res_plots[[i]], aes(x=bin, y=val, colour = color, linetype = line_type)) + stat_smooth(aes(x = bin, y = val), method = "lm",  formula = y ~ poly(x, 21), se = FALSE, n = 1000) + scale_colour_manual(values = c("green", "blue")) + theme_bw() + xlim(-950, 950)
  #     pdf(paste(folder_output, paste("Novel_Methylation/", names(l.res_plots)[i],".pdf", sep = ""), sep ="/"), width = 6, height = 4)
  #     print(p)
  #     dev.off()
  #   }
  # }
  # 
  # for(i in 1:length(l.res_plots)){
  #   p <- ggplot(l.res_plots[[i]], aes(x=bin, y=val, colour = color, linetype = line_type)) + stat_smooth(aes(x = bin, y = val), method = "lm",  formula = y ~ poly(x, 21), se = FALSE, n = 1000) + scale_colour_manual(values = c("green", "blue")) + theme_bw() + xlim(-950, 950)
  #   pdf(paste(folder_output, paste("Novel_Methylation/", names(l.res_plots)[i],".pdf", sep = ""), sep ="/"), width = 6, height = 4)
  #   print(p)
  #   dev.off()
  # }
  # 
}



check_methylation <- function(df.bQTLs, 
                              df.methylation,
                              th.bp_offset = 20,
                              parent = "B73",
                              has_ranges = FALSE,
                              n.cpus = 2, 
                              n.chromosomes = 10,
                              n.snps.split = 50){
  
  df.bQTLs["id"] <- seq(1, nrow(df.bQTLs))
  
  if(parent == "B73"){
    chrs = paste("B73-chr",1:10, sep = "")
    pos = "B73-pos"
    chr.parent = "B73-chr"
  }else{
    chrs = paste("Mo17-chr",1:10, sep = "")
    pos = "Mo17-pos"
    chr.parent = "Mo17-chr"
  }

  l.bQTLs <- l.methylation <- vector(mode = "list", length = n.chromosomes)
  names(l.bQTLs) <- names(l.methylation) <- chrs
  for(i in 1:n.chromosomes){
    l.bQTLs[[i]] <- subset(df.bQTLs, df.bQTLs[,chr.parent] == chrs[i])
    l.bQTLs[[i]] <- l.bQTLs[[i]][order(l.bQTLs[[i]][,pos]),]
    
    l.methylation[[i]] <- subset(df.methylation, df.methylation$chr == chrs[i])
    l.methylation[[i]] <- l.methylation[[i]][order(l.methylation[[i]]$start),]
    l.methylation[[i]] <- subset(l.methylation[[i]], l.methylation[[i]]$start >= min(l.bQTLs[[i]][,pos] - th.bp_offset) 
                                                   & l.methylation[[i]]$start <= max(l.bQTLs[[i]][,pos]) + th.bp_offset)
  }
  
  strt<-Sys.time() 
  cl<-makeCluster(n.cpus)
  registerDoParallel(cl)
  
  l.vals.methylation <- foreach(i = 1:n.chromosomes, .packages=c("seqinr", "VariantAnnotation", "Biostrings")) %dopar% {   

    v.pos.methylation <- l.methylation[[i]]$start
    v.pos.bQTLs <- l.bQTLs[[i]][,pos]

    v.vals.methylation <- numeric(length(v.pos.bQTLs))
  
    p.start <- v.pos.bQTLs[1] - th.bp_offset
    p.end <- v.pos.bQTLs[min(length(v.pos.bQTLs), n.snps.split)]
    
    idx.split <- which(p.start <= v.pos.methylation & v.pos.methylation <= p.end)
    v.pos.methylation.split <- v.pos.methylation[idx.split]
    v.vals.methylation.split <- l.methylation[[i]]$val[idx.split]
    
    for(j in 1:length(v.pos.bQTLs)){
      
      if(j %% n.snps.split == 0){
        p.start <- v.pos.bQTLs[j] - th.bp_offset
        p.end <- v.pos.bQTLs[min(length(v.pos.bQTLs), (j + n.snps.split))]
        
        idx.split <- which(p.start <= v.pos.methylation & v.pos.methylation <= p.end)
        
        v.pos.methylation.split <- v.pos.methylation[idx.split]
        v.vals.methylation.split <- l.methylation[[i]]$val[idx.split]
      }  
      
      idx = which(abs(v.pos.methylation.split - v.pos.bQTLs[j]) <= th.bp_offset)
      if(length(idx) > 0){
        v.vals.methylation[j] = mean(v.vals.methylation.split[idx])
      }
    }
    
    names(v.vals.methylation) <- l.bQTLs[[i]]$id
    
    v.vals.methylation
  }
  stopCluster(cl)
  print(Sys.time()-strt)
  
  v.vals.methylation = c()
  for(i in 1:n.chromosomes){
    v.vals.methylation <- c(v.vals.methylation, l.vals.methylation[[i]])
  }
  
  v.vals.methylation
}


#' A result statistics
#'
#' This function analysis the results
#' @param df.cluster_annotations inferred cluster condition annotations 
#' @param m.functionality gene cluster acitivity heatmap
#' @param th.bp_offset (default = 20)
#' @keywords 
#' @export
#' @examples
#' cat_function()
methylation_occupancy = function(df.ASBs, 
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
                                 th.bp_offset = 20,
                                 n.chromosomes = 10,
                                 n.cpus = 1,
                                 folder = "tmp"){

  message("(differential) methylation values versus post frequency")
  
  for(k in 1:length(v.filenames)){ # TODO: change back
    message("Processing, ", v.filenames[k])
    chrs <- paste(v.parents[k], "-chr", 1:n.chromosomes, sep = "")

    if(v.formats[k] == "bedGraph"){
      df.methylation <- read.table(v.filenames[k], header = FALSE, sep ="\t", quote = "", stringsAsFactors = FALSE)
    }else if(v.formats[k] == "bigWig"){
      df.methylation <- import.bw(v.filenames[k])
      df.methylation <- as.data.frame(df.methylation)
      df.methylation <- df.methylation[,c(1,2,6)]
    }
    names(df.methylation) <- c("chr",  "start", "val")
    df.methylation <- subset(df.methylation, df.methylation$chr %in% chrs)
    df.methylation = subset(df.methylation, df.methylation$val > 0) 
    gc()
    
    message("Checking ", nrow(df.methylation), " entries of methylation for ", nrow(df.ASBs), " ASBs")
    v.vals.methylation <- check_methylation(df.ASBs, 
                                            df.methylation,
                                            th.bp_offset = th.bp_offset,
                                            parent = v.parents[k],
                                            has_ranges = v.has_ranges[k],
                                            n.cpus = n.cpus, 
                                            n.chromosomes = n.chromosomes,
                                            n.snps.split = 50)
    
    gc()
    saveRDS(v.vals.methylation, paste(folder, paste("v.vals.methylation.ASBs.",v.parents[k],".", v.variants[k], ".rds", sep = ""), sep ="/"))
    
    message("Checking ", nrow(df.methylation), " entries of methylation for ", nrow(df.bgSNPs), " bgSNPs")
    v.vals.methylation <- check_methylation(df.bgSNPs, 
                                            df.methylation,
                                            th.bp_offset = th.bp_offset,
                                            parent = v.parents[k],
                                            has_ranges = v.has_ranges[k],
                                            n.cpus = n.cpus, 
                                            n.chromosomes = n.chromosomes,
                                            n.snps.split = 1000)
    
    saveRDS(v.vals.methylation, paste(folder, paste("v.vals.methylation.bgSNPs.",v.parents[k],".", v.variants[k], ".rds", sep = ""), sep ="/"))
    
    rm(df.methylation)
    gc()
  }
  
}


load_methylation = function(df.ASBs, 
                           df.bgSNPs,
                           v.parents = c("B73", "Mo17", "B73", "Mo17", "B73", "Mo17"),
                           v.variants = c("CG", "CG", "CHH", "CHH", "CHG", "CHG"),
                           folder = "tmp"){
  
  v.fields <-paste(v.variants,v.parents, sep ="-")
  df.ASBs[,v.fields] <- NA
  df.bgSNPs[,v.fields] <- NA
  
  for(k in 1:length(v.variants)){
    
    field <- v.fields[k]
    parent <- v.parents[k]

    v.vals.methylation <- readRDS(paste(folder, paste("v.vals.methylation.ASBs.",v.parents[k],".", v.variants[k], ".rds", sep = ""), sep ="/"))
    v.snp_ids <- as.numeric(names(v.vals.methylation))
    df.ASBs[v.snp_ids, field] <- v.vals.methylation
    
    v.vals.methylation <- readRDS(paste(folder, paste("v.vals.methylation.bgSNPs.",v.parents[k],".", v.variants[k], ".rds", sep = ""), sep ="/"))
    v.snp_ids <- as.numeric(names(v.vals.methylation))
    df.bgSNPs[v.snp_ids, field] <- v.vals.methylation
    
  }
  
  return(list(df.ASBs=df.ASBs, df.bgSNPs=df.bgSNPs))
}




# 
# plot_methylation_vs_allelic_bias = function(df.ASBs,
#                                             l.SNPs_w_chromatin,
#                                             
#                                             th.quantile = 0.95,
#                                             
#                                             th.bp_offset = 20
# ){
#   
#   # Test plot background or ASBs
#   rows = c("bp_window", "number_ASBs", "high_bias", "with_methylation_data", "have_methylation", "have_motifs", 
#            "have_motif_and_methylation", "have_motif_and_no_methylation", 
#            "have_no_motif_and_methylation", "have_no_motif_and_no_methylation",
#            "have_methylation_and_high_bias", "have_motifs_and_high_bias", 
#            "have_motif_and_methylation_and_high_bias", "have_motif_and_no_methylation_and_high_bias",
#            "have_no_motif_and_methylation_and_high_bias", "have_no_motif_and_no_methylation_and_high_bias")
#   
#   df.methylation_summary <- data.frame(matrix(nrow = length(rows), ncol = 1))
#   rownames(df.methylation_summary) = rows
#   
#   l.SNP_methylated <- readRDS(paste("l.ASBs_methylated.rds", sep ="")) # load!!
#   df.ASBs[,"has_motif"] <- df.ASBs_w_motifs$has_motif # needs the motif ASBs
#   
#   df.SNPs <- l.SNPs[[1]]
#   
#   idx_hb = which(df.SNPs$POSTfreq >= 0.85 | df.SNPs$POSTfreq <= 0.15)
#   high_bias = length(idx_hb)
#   
#   for(k in 1:length(offsets)){ # number of available ASBs, yes and no, subset by motifs... different offsets ... => create a dataframe
#     
#     idx1 = (l.SNP_methylated[[k]]$MO17_CpG >= 0.7 & l.SNP_methylated[[k]]$B73_CpG <= 0.3) 
#     idx2 = (l.SNP_methylated[[k]]$MO17_CpG <= 0.3 & l.SNP_methylated[[k]]$B73_CpG >= 0.7)
#     idx3 = (l.SNP_methylated[[k]]$MO17_CHG >= 0.7 & l.SNP_methylated[[k]]$B73_CHG <= 0.3)
#     idx4 = (l.SNP_methylated[[k]]$MO17_CHG <= 0.3 & l.SNP_methylated[[k]]$B73_CHG >= 0.7)
#     idx5 = (l.SNP_methylated[[k]]$MO17_CHH >= 0.7 & l.SNP_methylated[[k]]$B73_CHH <= 0.3)
#     idx6 = (l.SNP_methylated[[k]]$MO17_CHH <= 0.3 & l.SNP_methylated[[k]]$B73_CHH >= 0.7)
#     idx = idx1 | idx2 | idx3 | idx4 | idx4 | idx5 | idx6
#     
#     message("------------------")
#     message("results for offset, ",  offsets[k], "bp distance methylation position to ASB")
#     
#     perc = round(length(idx[!is.na(idx)]) / length(idx) * 100, 1)
#     
#     w_methylation_data = length(idx[!is.na(idx)])
#     message("methylation data available for ", perc, "% (", w_methylation_data, ") ASBs")
#     
#     # message("methylated (for available):")
#     df.ASBs[,"is_methylated"] = idx
#     have_methylation = as.numeric(table(idx)["TRUE"])
#     have_motifs = as.numeric(table(df.ASBs$has_motif)["TRUE"])
#     
#     res = df.ASBs$has_motif & df.ASBs$is_methylated
#     have_motif_and_methylation = as.numeric(table(res)["TRUE"])
#     
#     res = df.ASBs$has_motif & !df.ASBs$is_methylated
#     have_motif_and_no_methylation = as.numeric(table(res)["TRUE"])
#     
#     res = !df.ASBs$has_motif & df.ASBs$is_methylated
#     have_no_motif_and_methylation = as.numeric(table(res)["TRUE"])
#     
#     res = !df.ASBs$has_motif & !df.ASBs$is_methylated
#     have_no_motif_and_no_methylation = as.numeric(table(res)["TRUE"])
#     
#     # check with high bias
#     have_methylation_and_high_bias = as.numeric(table(idx[idx_hb])["TRUE"])
#     have_motifs_and_high_bias = as.numeric(table(df.ASBs[idx_hb,]$has_motif)["TRUE"])
#     
#     res = df.ASBs[idx_hb,]$has_motif & df.ASBs[idx_hb,]$is_methylated
#     have_motif_and_methylation_and_high_bias = as.numeric(table(res)["TRUE"])
#     
#     res = df.ASBs[idx_hb,]$has_motif & !df.ASBs[idx_hb,]$is_methylated
#     have_motif_and_no_methylation_and_high_bias = as.numeric(table(res)["TRUE"])
#     
#     res = !df.ASBs[idx_hb,]$has_motif & df.ASBs[idx_hb,]$is_methylated
#     have_no_motif_and_methylation_and_high_bias = as.numeric(table(res)["TRUE"])
#     
#     res = !df.ASBs[idx_hb,]$has_motif & !df.ASBs[idx_hb,]$is_methylated
#     have_no_motif_and_no_methylation_and_high_bias = as.numeric(table(res)["TRUE"])
#     
#     vals <- c(offsets[k], nrow(df.ASBs), high_bias, w_methylation_data, have_methylation, have_motifs, 
#               have_motif_and_methylation, have_motif_and_no_methylation, 
#               have_no_motif_and_methylation, have_no_motif_and_no_methylation,
#               have_methylation_and_high_bias, have_motifs_and_high_bias, 
#               have_motif_and_methylation_and_high_bias, have_motif_and_no_methylation_and_high_bias, 
#               have_no_motif_and_methylation_and_high_bias, have_no_motif_and_no_methylation_and_high_bias)
#     
#     df <- data.frame(vals)
#     rownames(df) <- as.character(rows)
#     
#     df.methylation_summary <- cbind(df.methylation_summary, df)
#     
#     message("---------------")
#   }
#   
#   df.methylation_summary = df.methylation_summary[,-1]
#   colnames(df.methylation_summary) <- as.character(offsets)
#   
#   write.csv(df.methylation_summary, "df.methylation_summary.csv")
#   
#   # methylierung scatter plots?
#   
#   dev.off()
#   k = 2
#   title = paste("allelic biases vs methylation levels (averaged per SNP using ",  offsets[k], "bp window)", sep = "")
#   par(mfrow=c(3,2))
#   scatter_plot(df.SNPs$POSTfreq, l.SNP_methylated[[k]]$B73_CpG, xlab = "allelic bias", ylab = "methylation level (B73 CpG)")
#   scatter_plot(df.SNPs$POSTfreq, l.SNP_methylated[[k]]$MO17_CpG, xlab = "allelic bias", ylab = "methylation level (MO17 CpG)")
#   scatter_plot(df.SNPs$POSTfreq, l.SNP_methylated[[k]]$B73_CHG, xlab = "allelic bias", ylab = "methylation level (B73 CHG)")
#   scatter_plot(df.SNPs$POSTfreq, l.SNP_methylated[[k]]$MO17_CHG, xlab = "allelic bias", ylab = "methylation level (MO17 CHG)")
#   scatter_plot(df.SNPs$POSTfreq, l.SNP_methylated[[k]]$B73_CHH, xlab = "allelic bias", ylab = "methylation level (B73 CHH)")
#   scatter_plot(df.SNPs$POSTfreq, l.SNP_methylated[[k]]$MO17_CHH, xlab = "allelic bias", ylab = "methylation level (MO17 CHH)")
#   mtext(title, side = 3, line = -2, outer = TRUE, cex = 0.8)
#   
#   # repeat with displaying only methylation explained ASBs
#   dev.off()
#   title = paste("allelic biases vs methylation levels (averaged per SNP using ",  offsets[k], "bp window)", sep = "")
#   par(mfrow=c(2,2))
#   idx = which(l.SNP_methylated[[k]]$MO17_CpG <= 0.3 & l.SNP_methylated[[k]]$B73_CpG >= 0.7)
#   scatter_plot(df.SNPs$POSTfreq[idx], l.SNP_methylated[[k]]$B73_CpG[idx], xlab = "allelic bias", ylab = "methylation level (B73 CpG)")
#   
#   idx = which(l.SNP_methylated[[k]]$MO17_CpG >= 0.7 & l.SNP_methylated[[k]]$B73_CpG <= 0.3)
#   scatter_plot(df.SNPs$POSTfreq[idx], l.SNP_methylated[[k]]$MO17_CpG[idx], xlab = "allelic bias", ylab = "methylation level (MO17 CpG)")
#   
#   idx = which(l.SNP_methylated[[k]]$MO17_CHG <= 0.3 & l.SNP_methylated[[k]]$B73_CHG >= 0.7)
#   scatter_plot(df.SNPs$POSTfreq[idx], l.SNP_methylated[[k]]$B73_CHG[idx], xlab = "allelic bias", ylab = "methylation level (B73 CHG)")
#   
#   idx = which(l.SNP_methylated[[k]]$MO17_CHG >= 0.7 & l.SNP_methylated[[k]]$B73_CHG <= 0.3)
#   scatter_plot(df.SNPs$POSTfreq[idx], l.SNP_methylated[[k]]$MO17_CHG[idx], xlab = "allelic bias", ylab = "methylation level (MO17 CHG)")
#   # 
#   # idx = which(l.SNP_methylated[[k]]$MO17_CHH <= 0.3 & l.SNP_methylated[[k]]$B73_CHH >= 0.7)
#   # scatter_plot(df.SNPs$POSTfreq[idx], l.SNP_methylated[[k]]$B73_CHH[idx], xlab = "allelic bias", ylab = "methylation level (B73 CHH)")
#   # 
#   # idx = which(l.SNP_methylated[[k]]$MO17_CHH >= 0.7 & l.SNP_methylated[[k]]$B73_CHH <= 0.3)
#   # scatter_plot(df.SNPs$POSTfreq[idx], l.SNP_methylated[[k]]$MO17_CHH[idx], xlab = "allelic bias", ylab = "methylation level (MO17 CHH)")
#   mtext(title, side = 3, line = -2, outer = TRUE, cex = 0.8)
#   
#   
#   ### 
#   
#   dev.off()
#   k = 2
#   title = paste("allelic biases vs methylation levels (averaged per SNP using ",  offsets[k], "bp window)", sep = "")
#   par(mfrow=c(3,1))
#   
#   idx1 = (l.SNP_methylated[[k]]$MO17_CpG <= 0.1 & l.SNP_methylated[[k]]$B73_CpG >= 0.7)
#   idx2 = (l.SNP_methylated[[k]]$MO17_CpG >= 0.7 & l.SNP_methylated[[k]]$B73_CpG <= 0.1)
#   idx = which(idx1 | idx2)
#   col_vals = rep("black", nrow(df.SNPs)) 
#   col_vals[idx] = "red"
#   
#   scatter_plot(df.SNPs$POSTfreq, l.SNP_methylated[[k]]$MO17_CpG - l.SNP_methylated[[k]]$B73_CpG, xlab = "allelic bias", ylab = "methylation level difference (CpG)", 
#                x_limit = c(0,1), y_limit = c(-1,1), col_vals = col_vals)
#   
#   ### 
#   idx1 = (l.SNP_methylated[[k]]$MO17_CHG <= 0.1 & l.SNP_methylated[[k]]$B73_CHG >= 0.7)
#   idx2 = (l.SNP_methylated[[k]]$MO17_CHG >= 0.7 & l.SNP_methylated[[k]]$B73_CHG <= 0.1)
#   idx = which(idx1 | idx2)
#   col_vals2 = rep("black", nrow(df.SNPs)) 
#   col_vals2[idx] = "red"
#   
#   scatter_plot(c(df.SNPs$POSTfreq, df.SNPs$POSTfreq), 
#                c(l.SNP_methylated[[k]]$MO17_CpG - l.SNP_methylated[[k]]$B73_CpG,
#                  l.SNP_methylated[[k]]$MO17_CHG - l.SNP_methylated[[k]]$B73_CHG), 
#                xlab = "allelic bias", ylab = "methylation level difference - CpG (red), CHG (blue))", 
#                x_limit = c(0,1), y_limit = c(-1,1), col_vals = c(col_vals, col_vals2))
#   
#   
#   ###
#   idx1 = (l.SNP_methylated[[k]]$MO17_CHG <= 0.1 & l.SNP_methylated[[k]]$B73_CHG >= 0.7)
#   idx2 = (l.SNP_methylated[[k]]$MO17_CHG >= 0.7 & l.SNP_methylated[[k]]$B73_CHG <= 0.1)
#   idx = which(idx1 | idx2)
#   col_vals = rep("black", nrow(df.SNPs)) 
#   col_vals[idx] = "red"
#   
#   scatter_plot(df.SNPs$POSTfreq, l.SNP_methylated[[k]]$MO17_CHG - l.SNP_methylated[[k]]$B73_CHG, xlab = "allelic bias", ylab = "methylation level difference (CHG)", 
#                x_limit = c(0,1), y_limit = c(-1,1), col_vals = col_vals)
#   
#   ###
#   
#   idx1 = (l.SNP_methylated[[k]]$MO17_CHH <= 0.1 & l.SNP_methylated[[k]]$B73_CHH >= 0.7)
#   idx2 = (l.SNP_methylated[[k]]$MO17_CHH >= 0.7 & l.SNP_methylated[[k]]$B73_CHH <= 0.1)
#   idx = which(idx1 | idx2)
#   col_vals = rep("black", nrow(df.SNPs)) 
#   col_vals[idx] = "red"
#   
#   scatter_plot(df.SNPs$POSTfreq, l.SNP_methylated[[k]]$MO17_CHH - l.SNP_methylated[[k]]$B73_CHH, xlab = "allelic bias", ylab = "methylation level difference (CHH)", 
#                x_limit = c(0,1), y_limit = c(-1,1), col_vals = col_vals)
#   
#   
#   
#   scatter_plot(df.SNPs$POSTfreq, l.SNP_methylated[[k]]$MO17_CHG - l.SNP_methylated[[k]]$B73_CHG, xlab = "allelic bias", ylab = "methylation level difference (CHG)")
#   scatter_plot(df.SNPs$POSTfreq, l.SNP_methylated[[k]]$MO17_CHH - l.SNP_methylated[[k]]$B73_CHH, xlab = "allelic bias", ylab = "methylation level difference (CHH)")
#   mtext(title, side = 3, line = -2, outer = TRUE, cex = 0.8)
#   
#   
# }
# 

                                                    
                                                    
# differential_methylation_vs_allelic_bias = function(df.ASBs,
#                                                     l.SNPs_w_chromatin,
#                                                     th.quantile = 0.95,
#                                                     th.bp_offset = 20,
#                                                     th.distance_to_ASB = 2000,
#                                                     n.chromosomes = 10,
#                                                     n.cpus = 1,
#                                                     b.val = TRUE,
#                                                     v.species = c("Mo17", "B73"),
#                                                     degrees = 15,
#                                                     th.padding = 50,
#                                                     width = 6,
#                                                     height = 4,
#                                                     group = "genic_and_nongenic",
#                                                     folder_tmp = "",
# ){ 
# 
#   l.diff_data = vector(mode = "list", length = 3)
#   l.diff_data[[1]] = l.SNPs_w_chromatin[[1]][[3]] - l.SNPs_w_chromatin[[1]][[2]]
#   
#   l.diff_SNPs_w_chromatin <- vector(mode = "list", length = 2)
#   l.diff_SNPs_w_chromatin[[1]] <- vector(mode = "list", length = 4)
#   l.diff_SNPs_w_chromatin[[2]] <- vector(mode = "list", length = 4)
#   
#   l.diff_SNPs_w_chromatin[[1]][[1]] = l.SNPs_w_chromatin[[1]][[1]]$DNase
#   l.diff_SNPs_w_chromatin[[2]][[1]] = l.SNPs_w_chromatin[[2]][[1]]$DNase
#   
#   l.diff_SNPs_w_chromatin[[1]][[2]] = l.SNPs_w_chromatin[[1]][[3]]$MO17_CpG - l.SNPs_w_chromatin[[1]][[2]]$B73_CpG
#   l.diff_SNPs_w_chromatin[[2]][[2]] = l.SNPs_w_chromatin[[2]][[3]]$MO17_CpG - l.SNPs_w_chromatin[[2]][[2]]$B73_CpG
#   
#   l.diff_SNPs_w_chromatin[[1]][[3]] = l.SNPs_w_chromatin[[1]][[5]]$MO17_CHG - l.SNPs_w_chromatin[[1]][[4]]$B73_CHG
#   l.diff_SNPs_w_chromatin[[2]][[3]] = l.SNPs_w_chromatin[[2]][[5]]$MO17_CHG - l.SNPs_w_chromatin[[2]][[4]]$B73_CHG
#   
#   l.diff_SNPs_w_chromatin[[1]][[4]] = l.SNPs_w_chromatin[[1]][[7]]$MO17_CHH - l.SNPs_w_chromatin[[1]][[6]]$B73_CHH
#   l.diff_SNPs_w_chromatin[[2]][[4]] = l.SNPs_w_chromatin[[2]][[7]]$MO17_CHH - l.SNPs_w_chromatin[[2]][[6]]$B73_CHH
#   
#   data = paste(c("DNase  >=", "CpG (MO17 - B73) >=", "CHG  (MO17 - B73)  >=", "CHH  (MO17 - B73)  >="), th.quantile)
#   for(i in 1:4){
#     th = quantile(l.diff_SNPs_w_chromatin[[2]][[i]], th.quantile, na.rm = T)
#     idx = which(l.diff_SNPs_w_chromatin[[1]][[i]] >= th)
#     df.ASBs[,data[i]] = "no"
#     df.ASBs[idx,data[i]] = "yes"
#   }
#   
#   data = paste(c("CpG (MO17 - B73) <=", "CHG  (MO17 - B73)  <=", "CHH  (MO17 - B73)  <="), 1 - th.quantile)
#   for(i in 2:4){
#     th = quantile(l.diff_SNPs_w_chromatin[[2]][[i]], 1 - th.quantile, na.rm = T)
#     idx = which(l.diff_SNPs_w_chromatin[[1]][[i]] <= th)
#     df.ASBs[,data[i-1]] = "no"
#     df.ASBs[idx,data[i-1]] = "yes"
#   }
# 
#   
#   
#   v.datasets <- c("DNase", "CpG", "CHG", "CHH")
#   for(k in 1:length(v.datasets)){
#     
#     # differential methylation grouped density (ASBs, bgSNPs)
#     df = data.frame(val = c(l.diff_SNPs_w_chromatin[[k]][[1]], l.diff_SNPs_w_chromatin[[k]][[2]]),
#                     group = c(rep("ASBs", length(l.diff_SNPs_w_chromatin[[k]][[1]])), 
#                               rep("bg SNPs", length(l.diff_SNPs_w_chromatin[[k]][[2]]))))
#     
#     ggplot(df, aes(x=val, fill=group)) + geom_density(alpha=0.4) + theme_bw()
#     file = paste(folder_tmp, paste(v.datasets[k], ".pdf", sep = ""), sep = "/")
#     ggsave(file, width = 20, height = 20, units = "cm")
#     
#     
#     # differential methylation vs postfrequency (ABSs only)
#     
#     
#    # scatter_plot(df.ASBs$POSTfreq, l.diff_SNPs_w_chromatin[[k]][[1]], y_limit = c(0, 1.5))
#   }
#   
#   
#   
#   
#   
#   
# 
#   
#   
#   
#   
#   
#   
#   
#   l.ASB_SNPs_w_chromatin = readRDS("l.SNPs_w_chromatin_ASBs.rds")
#   l.bgSNPs_w_chromatin = readRDS("l.bgSNPs_w_chromatin.rds")
#   
#   diff_ASB_MO17_minus_B73 = l.ASB_SNPs_w_chromatin[[3]]$MO17_CpG - l.ASB_SNPs_w_chromatin[[2]]$B73_CpG
#   diff_bgSNPs_MO17_minus_B73 = l.bgSNPs_w_chromatin[[3]]$MO17_CpG - l.bgSNPs_w_chromatin[[2]]$B73_CpG
#   
#   th.under = quantile(diff_bgSNPs_MO17_minus_B73, 0.05, na.rm = T)
#   th.over = quantile(diff_bgSNPs_MO17_minus_B73, 0.95, na.rm = T)
#   
#   df.ASBs =  l.ASB_SNPs_w_chromatin[[1]]
#   df.ASBs["delta_MO17_minus_B73"] = diff_ASB_MO17_minus_B73
#   df.ASBs["significance_MO17_minus_B73"] = NA
#   
#   idx = which(df.ASBs$delta_MO17_minus_B73 >= th.over)
#   df.ASBs$significance_MO17_minus_B73[idx] = "upper_95"
#   
#   idx = which(df.ASBs$delta_MO17_minus_B73 <= th.under)
#   df.ASBs$significance_MO17_minus_B73[idx] = "lower_005"
#   
#   
#   
#   df.ASBs["DNase"] = l.ASB_SNPs_w_chromatin[[1]]$DNase
#   df.ASBs["significance_DNase"] = NA
#   
#   th.over = quantile((l.bgSNPs_w_chromatin[[1]]$DNase), 0.95, na.rm = T)
#   
#   idx = which(df.ASBs$DNase >= th.over)
#   df.ASBs$significance_DNase[idx] = "upper_95"
#   
#   
#   
#   hist(diff_ASB_MO17_minus_B73, 100)
#   
#   # for(k in 1:length(v.datasets)){
#   #   write.csv2(l.results[[k]][[1]], paste("D:/HashSeq/tmp/", v.datasets[k], "_ASBs.csv", sep = ""), row.names = F)
#   #   write.csv2(l.results[[k]][[2]], paste("D:/HashSeq/tmp/", v.datasets[k], "_bgSNPs.csv", sep = ""), row.names = F)
#   # }
#   
#   
#   
#   
#   
#   ### 
#   df.ASBs <- l.SNPs[[1]]
#   for(j in 1:length(v.datasets)){
#     df.ASBs[v.datasets[j]] <- 1
#   }
#   
#   
#   l.SNPs_w_chromatin = readRDS("l.SNPs_w_chromatin_ASBs.rds")
#   
#   
#   for(k in 1:length(v.datasets)){
#     
#     df.bg_data = read.csv2(paste("D:/HashSeq/tmp/", v.datasets[k], "_bgSNPs.csv", sep = ""))
#     
#     th = quantile(df.bg_data$val, 0.95)
#     
#     l.SNPs_w_chromatin[[k]][">95BG"] = ifelse(l.SNPs_w_chromatin[[k]][,v.datasets[k]] >= th, "yes", "no")
#     
#     
#     df.ASBs[,v.datasets[k]] =  ifelse(l.SNPs_w_chromatin[[k]][,v.datasets[k]] >= th, "yes", "no")
#     # 
#     
#     # df.ASBs[v.datasets[k]] = 
#   }
#   
#   write.csv2(df.ASBs, "A:/junkDNA.ai/df.ASBs.csv",row.names = F)
#   
#   
#   table(l.SNPs_w_chromatin[[3]]$`>95BG`)
#   
#   
#   ### with RNASEQ?
#   
#   # df.bQTL_gene_partitioning <- df.ASBs
#   # df.ASB_w_RNASeq <-  merge(df.bQTL_gene_partitioning, df.rnaseq.B73vsMo17, by.x = "gene.ID", by.y = "Gene_ID_AGPv4", all = FALSE)
#   # 
#   # df.ASBs["RNASeq <= 0.05"] = NA
#   # df.ASBs["RNASeq (foldchange)"] = NA
#   # idx = which(df.ASBs$gene.ID %in% df.rnaseq.B73vsMo17$Gene_ID_AGPv4)
#   # for(i in 1:nrow(df.ASBs)){
#   #   idx = which(df.rnaseq.B73vsMo17$Gene_ID_AGPv4 == df.ASBs$gene.ID[i])[1]
#   #   if(length(idx) > 0 & !is.na(idx)){
#   #     
#   #     pval = t.test(c(df.rnaseq.B73vsMo17$B73_control_rep1_batch1[idx],df.rnaseq.B73vsMo17$B73_control_rep2_batch1[idx],df.rnaseq.B73vsMo17$B73_control_rep2_batch1[idx]),
#   #                   c(df.rnaseq.B73vsMo17$Mo17_control_rep1_batch2[idx],df.rnaseq.B73vsMo17$Mo17_control_rep2_batch2[idx],df.rnaseq.B73vsMo17$Mo17_control_rep3_batch2[idx]))$p.val
#   #     
#   #     if(pval <= 0.05 & !is.na(pval)){
#   #       df.ASBs$`RNASeq <= 0.05`[i] = "yes"
#   #     }else{
#   #       df.ASBs$`RNASeq <= 0.05`[i] = "no"
#   #     } 
#   #     
#   #     df.ASBs$`RNASeq (foldchange)`[i] = df.rnaseq.B73vsMo17$B73_vs_Mo17[idx]
#   #   }
#   # }
#   # 
#   # data = c("DNase  >= 0.95", "CpG (MO17 - B73) >= 0.95", "CHG  (MO17 - B73)  >= 0.95", "CHH  (MO17 - B73)  >= 0.95", 
#   #          "CpG (MO17 - B73) <= 0.05", "CHG  (MO17 - B73)  <= 0.05", "CHH  (MO17 - B73)  <= 0.05", "RNASeq <= 0.05")
#   # 
#   # df.ASBs["any"] = ifelse(apply(df.ASBs[,data], 1, function(m) { any(m == "yes") }) == TRUE, "yes", "no")
#   # write.table(df.ASBs, "df.ASBs_DNaseMethylationRNASeq.txt", row.names = F, sep = "\t")
#   
# }
# 


