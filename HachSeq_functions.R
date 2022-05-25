save.xlsx <- function (file, ...)
{
  require(xlsx, quietly = TRUE)
  objects <- list(...)
  fargs <- as.list(match.call(expand.dots = TRUE))
  objnames <- as.character(fargs)[-c(1, 2)]
  nobjects <- length(objects)
  for (i in 1:nobjects) {
    if (i == 1)
      write.xlsx(objects[[i]], file, sheetName = objnames[i], row.names = FALSE)
    else write.xlsx(objects[[i]], file, sheetName = objnames[i],
                    append = TRUE, row.names = FALSE)
  }
  print(paste("Workbook", file, "has", nobjects, "worksheets."))
}

install_and_load_libraries <- function(){
  
  # CRAN
  list.of.packages <- c("ggplot2", "reshape2","doParallel", "pheatmap", "igraph", "seqinr", "foreach", "plotly", "rtracklayer")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  # bioconductor
  # source("https://bioconductor.org/biocLite.R")
  # install.packages("BiocManager")
  # BiocManager::install(version = '3.12')
  
  list.of.packages <- c("Biostrings", "VariantAnnotation","BSgenome")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) BiocManager::install(new.packages)
  
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
}



v.filenames <- c("datasets_paper/Enhancer_HM/GSE94251_H3K9ac_ist.bedGraph", 
                 "datasets_paper/Enhancer_HM/GSE94291_DNase_ist.bedGraph", 
                 "datasets_paper/Methylation/B73_CHG.bw", 
                 "datasets_paper/Methylation/B73_CHH.bw", 
                 "datasets_paper/Methylation/B73_CpG.bw", 
                 "datasets_paper/Methylation/MO17_CHG.bw", 
                 "datasets_paper/Methylation/MO17_CHH.bw", 
                 "datasets_paper/Methylation/MO17_CpG.bw")

v.datasets <- c("H3K9", "DNase",  "B73_CHG", "B73_CHH", "B73_CpG", "MO17_CHG", "MO17_CHH", "MO17_CpG")
v.formats <- c("bedGraph", "bedGraph", "bigWig", "bigWig", "bigWig", "bigWig", "bigWig", "bigWig")
v.has_ranges = c("yes", "yes", "no", "no", "no", "no", "no", "no")


l.SNPs <- readRDS(paste("tmp/l.bQTL_gene_partitioning_withGeneDistances_backgroundSampled_", timeStamp, ".rds", sep = ""))


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
analyse_asb_methylation_vs_postfrequency = function(l.SNPs=l.SNPs,
                                                    th.bp_offset = 20){
  
  
  
  th.distance_to_ASB <-  2000
  n.chromosomes <- 10
  
  n.cpus <- 5
  # n.binWidth = 20
  b.val = TRUE
  v.species = c("Mo17", "B73")

  v.binWidth = 
    degrees = 15
  th.distance_to_ASB = 5000
  th.padding = 50
  width = 6
  height = 4
  
  group = "genic_and_nongenic"
  

  l.SNPs <- readRDS(paste("D:/HashSeq/tmp/l.bQTL_gene_partitioning_withGeneDistances_backgroundSampled_117.rds", sep = ""))
  
  
  message("plot asb methylation values versus post frequency")
  
  df.ASBs <- l.SNPs[[1]]
  v.gns_ASBs = unique(df.ASBs$gene.ID)
  df.bpSNPs <- l.SNPs[[2]]
  df.bpSNPs <- subset(df.bpSNPs, !df.bpSNPs$gene.ID %in% v.gns_ASBs) # 
  
  l.SNP_selection = vector(mode = "list", length = 2)
  l.SNP_selection[[1]] = df.ASBs
  l.SNP_selection[[2]] = df.bpSNPs
  
  v.filenames <- c("D:/HashSeq/datasets_paper/Enhancer_HM/GSE94291_DNase_ist.bedGraph", 
                   "D:/HashSeq/datasets_paper/Methylation/B73_CpG.bw", 
                   "D:/HashSeq/datasets_paper/Methylation/MO17_CpG.bw")
  
  v.datasets <- c("DNase", "B73_CpG",  "MO17_CpG" ) 
  v.formats <- c("bedGraph", "bigWig", "bigWig")
  v.has_ranges = c("yes",  "no", "no")

  
  l.results <- vector(mode = "list", length = length(v.datasets))
    
  for(k in 1:length(v.datasets)){
  
    l.results[[k]] = vector(mode = "list", length = 2)
    
    message("Processing, ", v.datasets[k])
    if(v.formats[k] == "bedGraph"){
      df.dataset <- read.table(v.filenames[k], header = FALSE, sep ="\t", quote = "", stringsAsFactors = FALSE)
    }else if(v.formats[k] == "bigWig"){
      df.dataset <- import.bw(v.filenames[k])
      df.dataset <- as.data.frame(df.dataset)
      df.dataset <- df.dataset[,c(1,2,3,6)]
    }
    names(df.dataset) <- c("chr",  "start", "end", "val")
    
    if(v.datasets[k] == "DNase"){
      df.dataset = subset(df.dataset, df.dataset$val > 0) # remove lack in coverage 
    }
    
    for(s in 1:length(l.SNP_selection)){
      
      message(k, " / ", s)
      
      df.SNPs = l.SNP_selection[[s]]
      v.vals = numeric(nrow(df.SNPs))
      v.postfreqs = numeric(nrow(df.SNPs))
      idx_SNP = 1
      for(i in 1:n.chromosomes){
        df.SNPs.i <- subset(df.SNPs, df.SNPs$contig == i)
        pb <- txtProgressBar(min = 0, max = nrow(df.SNPs.i), style = 3)
        for(j in 1:nrow(df.SNPs.i)){
          setTxtProgressBar(pb, j)
          if(v.has_ranges[k] == "yes"){
            idx = which(df.dataset$start <= df.SNPs.i$position[j] & df.SNPs.i$position[j] <= df.dataset$end)  
          }else{ # what is the closest 
            idx = which(abs(df.dataset$start-df.SNPs.i$position[j]) <= th.bp_offset) # min(abs(df.dataset$start-df.SNPs.i$position[j])))
          }
          if(length(idx) > 0){
            v.vals[idx_SNP] = mean(df.dataset$val[idx])
            if(s == 1){
              v.postfreqs[idx_SNP] = df.SNPs.i$POSTfreq[j]
            }else{
              v.postfreqs[idx_SNP] = 0.5
            }
            idx_SNP = idx_SNP + 1
          }
        }
        close(pb)
      }
      
      v.vals = v.vals[1:idx_SNP]
      v.postfreqs = v.postfreqs[1:idx_SNP]
 
      df = data.frame(val = v.vals, postfreq = v.postfreqs)  
       
      l.results[[k]][[s]] = df 
      
    }
  }

  # saveRDS(l.results, "l.results.rds")
  
  
  
  for(k in 1:length(v.datasets)){
    
    df = data.frame(val = c(l.results[[k]][[1]]$val, l.results[[k]][[2]]$val),
                    group = c(rep("ASBs", nrow(l.results[[k]][[1]])), rep("bg SNPs", nrow(l.results[[k]][[2]]))))
    
    ggplot(df, aes(x=val, fill=group)) + geom_density(alpha=0.4) + theme_bw()

    file = paste("D:/HashSeq/tmp/", v.datasets[k], ".pdf", sep = "")
    ggsave(file, width = 20, height = 20, units = "cm")
    
  }
  
  for(k in 1:length(v.datasets)){
    write.csv2(l.results[[k]][[1]], paste("D:/HashSeq/tmp/", v.datasets[k], "_ASBs.csv", sep = ""), row.names = F)
    write.csv2(l.results[[k]][[2]], paste("D:/HashSeq/tmp/", v.datasets[k], "_bgSNPs.csv", sep = ""), row.names = F)
  }
  
  
}


#' A result statistics
#'
#' This function analysis the results
#' @param df.cluster_annotations inferred cluster condition annotations 
#' @param m.functionality gene cluster acitivity heatmap
#' @keywords 
#' @export
#' @examples
#' cat_function()
distances_Methylation_to_ASB  = function(l.bQTL_gene_partitioning, 
                                         th.distance_to_ASB = 5000,
                                         th.distance_to_ASB_plot = 2000,
                                         n.chromosomes = 10,
                                         n.cpus = 5,
                                         n.binWidth = 50,
                                         b.val = TRUE,
                                         v.species = c("Mo17", "B73"),
                                         v.levels = c(0.15,0.85),
                                         v.stringency = c( "< 0.15 | > 0.85",  "< 0.15 | > 0.85"),
                                         degrees = 15,
                                         th.padding = 50,
                                         width = 6,
                                         height = 4,
                                         group = "genic_and_nongenic"){
  
  
  l.bQTL_gene_partitioning <- readRDS(paste("D:/HashSeq/tmp/l.bQTL_gene_partitioning_withGeneDistances_backgroundSampled_117.rds", sep = ""))
  
 
  df.SNPs <-  l.bQTL_gene_partitioning[[1]]
  
  # separate by directionality - only for ASBs
  l.df.SNPs.directionality <- vector(mode = "list", length = 2)
  l.df.SNPs.directionality[[1]]  <- subset(df.SNPs, df.SNPs$POSTfreq <= v.levels[1])
  l.df.SNPs.directionality[[2]]  <- subset(df.SNPs, df.SNPs$POSTfreq >= v.levels[2])
  
  names(v.species) <- v.species
  names(v.stringency) <- v.species
  
  l.res_dists <- vector(mode = "list", length = length(v.datasets))
  
  v.filenames <- c("D:/HashSeq/datasets_paper/Enhancer_HM/GSE94291_DNase_ist.bedGraph", 
                   "D:/HashSeq/datasets_paper/Methylation/B73_CpG.bw", 
                   "D:/HashSeq/datasets_paper/Methylation/MO17_CpG.bw")
  
  v.datasets <- c("DNase", "B73_CpG",  "MO17_CpG" ) 
  v.formats <- c("bedGraph", "bigWig", "bigWig")
  v.has_ranges = c("yes",  "no", "no")
  
  l.results <- vector(mode = "list", length = length(v.datasets))
  
  for(k in 1:length(v.datasets)){

    message("Processing, ", v.datasets[k])
    if(v.formats[k] == "bedGraph"){
      df.dataset <- read.table(v.filenames[k], header = FALSE, sep ="\t", quote = "", stringsAsFactors = FALSE)
    }else if(v.formats[k] == "bigWig"){
      df.dataset <- import.bw(v.filenames[k])
      df.dataset <- as.data.frame(df.dataset)
      df.dataset <- df.dataset[,c(1,2,3,6)]
    }
    names(df.dataset) <- c("chr",  "start", "end", "val")
    
    if(v.datasets[k] == "DNase"){
      df.dataset = subset(df.dataset, df.dataset$val > 0) # remove lack in coverage 
    }
  
    # separating, all, genic, non genic
    df.distSet <- c()
    for(s in 1:length(l.df.SNPs.directionality)){
      df.SNPs <- l.df.SNPs.directionality[[s]]
      v.dist_to_ASB <- c()
      v.vals.dist_to_ASB <- c()
      
      for(i in 1:n.chromosomes){
        
        df.SNPs.i <- subset(df.SNPs, df.SNPs$contig == i)
        df.distanceDataset.i <- subset(df.dataset, df.dataset$chr == i)
        
        v.dist_to_ASB.i <- c()
        v.vals.dist_to_ASB.i <- c()
        
        for(j in 1:nrow(df.SNPs.i)){
          v.dist_to_ASB.j <- (df.distanceDataset.i$start - df.SNPs.i$position[j])
          # the histone dataset (H3K9) must be treated strand sensitive
          if(df.SNPs.i$strand[j] == "-"){# & k == 1){
            v.dist_to_ASB.j <- v.dist_to_ASB.j * -1
          }
          
          v.vals.dist_to_ASB.j <- df.distanceDataset.i$val
          
          idx <- which(abs(v.dist_to_ASB.j) <= th.distance_to_ASB)
          
          v.vals.dist_to_ASB.j <- v.vals.dist_to_ASB.j[idx]
          
          v.dist_to_ASB.j <- v.dist_to_ASB.j[idx]
          v.dist_to_ASB.i <- c(v.dist_to_ASB.i, v.dist_to_ASB.j)
          v.vals.dist_to_ASB.i <- c(v.vals.dist_to_ASB.i, v.vals.dist_to_ASB.j)
        }
        v.dist_to_ASB <- c(v.dist_to_ASB, v.dist_to_ASB.i)
        v.vals.dist_to_ASB <- c(v.vals.dist_to_ASB, v.vals.dist_to_ASB.i)
      }
    
      df.distSet <- rbind(df.distSet, data.frame(dist = as.numeric(v.dist_to_ASB), val = as.numeric(v.vals.dist_to_ASB), color = v.species[s], line_type = v.stringency[s]))
    }
    l.results[[k]] <- df.distSet
    
    rm(df.dataset)
  }
  
  saveRDS(l.results, paste("A:/junkDNA.ai/HaschSeq/datasets/Methylation_ASB_distances/l.data.rds", sep = ""))
  
  # store datasets 
  for(k in 1:length(v.datasets)){
    df.distSet <- l.results[[k]]
    write.table(df.distSet, paste("A:/junkDNA.ai/HaschSeq/results/Methylation_ASB_distances/df.distance_plots_",v.datasets[k], ".txt", sep = ""), sep = "\t", row.names = FALSE)
  }

 
  l.res_plots <- vector(mode = "list", length = length(v.datasets))
  for(k in 1:length(v.datasets)){
    df.distPlot <- c()
    df.distSet <- l.results[[k]]
    for(s in 1:length(l.df.SNPs.directionality)){
      df.SNPs <- l.df.SNPs.directionality[[s]]
      df.distSet.s <- subset(df.distSet, df.distSet$color == v.species[s] & df.distSet$line_type == v.stringency[s])
      v.dist_to_ASB <- df.distSet.s$dist
      v.vals.dist_to_ASB <- df.distSet.s$val
      v.breaks_sets <- seq(-th.distance_to_ASB, th.distance_to_ASB, n.binWidth)
      v.mean_per_break <- numeric(length(v.breaks_sets))
      names(v.mean_per_break) <- v.breaks_sets
      for(j in 1:(length(v.breaks_sets) - 1)){
        idx <- which(v.breaks_sets[j] <= v.dist_to_ASB & v.dist_to_ASB <= v.breaks_sets[j + 1])
        v.mean_per_break[j] <- mean(v.vals.dist_to_ASB[idx])
      }
      v.mean_per_break <- v.mean_per_break[-length(v.mean_per_break)]
      df.distPlot <- rbind(df.distPlot, data.frame(bin = as.numeric(names(v.mean_per_break)), val = v.mean_per_break, color = v.species[s], line_type = v.stringency[s], stringsAsFactors = FALSE))
    }
    l.res_plots[[k]] <- df.distPlot
  }
  names(l.res_plots) <- v.datasets 
  saveRDS(l.res_plots, paste("A:/junkDNA.ai/HaschSeq/datasets/Methylation_ASB_distances/l.chromatin_opening_distance_plots_", group,"_", "binwidth_",n.binWidth, ".rds", sep = ""))

  for(k in 1:length(v.datasets)){
    df.distPlot <- l.res_plots[[k]]
    write.table(df.distPlot, paste("A:/junkDNA.ai/HaschSeq/results/Methylation_ASB_distances/df.distance_plots_binned_",v.datasets[k], "_", group,"_", "binwidth_", n.binWidth, ".txt", sep = ""), sep = "\t", row.names = FALSE)
  }
  
  for(k in 1:length(v.datasets)){
    p <- ggplot(l.res_plots[[k]], aes(x=bin, y=val, colour = color, linetype = line_type)) + stat_smooth(aes(x = bin, y = val), method = "lm",  formula = y ~ poly(x, degree= degrees), se = FALSE, n = th.distance_to_ASB_plot) + scale_colour_manual(values = c("green", "blue")) + theme_bw() + xlim(-th.distance_to_ASB_plot + th.padding , th.distance_to_ASB_plot - th.padding)
    pdf(paste("A:/junkDNA.ai/HaschSeq/results/Methylation_ASB_distances/df.distance_plots_binned_",v.datasets[k], "_", group,"_", "binwidth_",n.binWidth, ".pdf", sep = ""), width = width, height = height)
    print(p)
    dev.off()
    pg = ggplot_build(p)
    df = data.frame(bin = pg$data[[1]]$x, val = pg$data[[1]]$y, color = pg$data[[1]]$colour, linetype = pg$data[[1]]$linetype)
    write.table(df, paste("A:/junkDNA.ai/HaschSeq/results/Methylation_ASB_distances/df.trendline_binned_",v.datasets[k], "_", group,"_", "binwidth_", n.binWidth, ".txt", sep = ""), sep = "\t", row.names = FALSE)
  }
  
  
    
  

  message("Distance plots with motif information")
  
  df.ASB.MOTIFS.only <- l.motif_analysis[[1]]
  df.ASB.MOTIFS.only <- unique(df.ASB.MOTIFS.only[,c(1,2)])
  
  # separate by directionality - only for ASBs
  df.SNPs <-  l.bQTL_gene_partitioning[[1]]
  l.df.SNPs.directionality <- vector(mode = "list", length = 6)
  l.df.SNPs.directionality[[1]] <- subset(df.SNPs, df.SNPs$POSTfreq <= 0.15)
  l.df.SNPs.directionality[[2]]  <- subset(df.SNPs, df.SNPs$POSTfreq <= 0.15)
  l.df.SNPs.directionality[[3]]  <- subset(df.SNPs, df.SNPs$POSTfreq <= 0.15)
  l.df.SNPs.directionality[[4]] <- subset(df.SNPs, df.SNPs$POSTfreq >= 0.85)
  l.df.SNPs.directionality[[5]]  <- subset(df.SNPs, df.SNPs$POSTfreq >= 0.85)
  l.df.SNPs.directionality[[6]]  <- subset(df.SNPs, df.SNPs$POSTfreq >= 0.85)
  
  v.species <- c("Mo17", "Mo17", "Mo17" , "B73", "B73", "B73")
  v.stringency <- c("< 0.15 | > 0.85", "< 0.15 | > 0.85 in motifs", "< 0.15 | > 0.85 not in motifs", "< 0.15 | > 0.85", "< 0.15 | > 0.85 in motifs", "< 0.15 | > 0.85 not in motifs")
  
  names(v.species) <- v.species
  names(v.stringency) <- v.species
  
  l.res_plots <- vector(mode = "list", length = length(v.datasets))
  for(k in 1:length(v.datasets)){
    
    message("Processing, ", v.datasets[k])
    
    if(v.formats[k] == "bedGraph"){
      df.dataset <- read.table(v.filenames[k], header = FALSE, sep ="\t", quote = "", stringsAsFactors = FALSE)
    }else if(v.formats[k] == "bigWig"){
      df.dataset <- import.bw(v.filenames[k])
      df.dataset <- as.data.frame(df.dataset)
      df.dataset <- df.dataset[,c(1,2,3,6)]
    }
    names(df.dataset) <- c("chr",  "start", "end", "val")
  
    df.distPlot <- c()
    for(s in 1:length(l.df.SNPs.directionality)){
      
      df.SNPs <- l.df.SNPs.directionality[[s]]
      
      v.dist_to_ASB <- c()
      v.vals.dist_to_ASB <- c()
      
      for(i in 1:n.chromosomes){
        
        # consider strand sensitivity of the H3K9 histone dataset
        df.SNPs.i <- subset(df.SNPs, df.SNPs$contig == i)
        
        # motif
        df.ASB.MOTIFS.only.i <- subset(df.ASB.MOTIFS.only, df.ASB.MOTIFS.only$contig == i)
        
        if(s == 2 | s == 5){
          df.SNPs.i <- subset(df.SNPs.i, df.SNPs.i$position %in% df.ASB.MOTIFS.only.i$position) # directionality filter
        }else if(s == 3 | s == 6){
          df.SNPs.i <- subset(df.SNPs.i, !df.SNPs.i$position %in% df.ASB.MOTIFS.only.i$position)
        }
        
        df.distanceDataset.i <- subset(df.dataset, df.dataset$chr == i)
        
        v.dist_to_ASB.i <- c()
        v.vals.dist_to_ASB.i <- c()
        
          for(j in 1:nrow(df.SNPs.i)){ 
            
            v.dist_to_ASB.j <- (df.distanceDataset.i$start - df.SNPs.i$position[j])
            
            # all datasets must be treated strand sensitive
            
            if(df.SNPs.i$strand[j] == "-"){#  & k == 1){
              v.dist_to_ASB.j <- v.dist_to_ASB.j * -1
            }
            
            # if(df.SNPs.i$strand[j] == "-" & k == 1){
            #   v.dist_to_ASB.j <- v.dist_to_ASB.j * -1
            # }
            
            v.vals.dist_to_ASB.j <- df.distanceDataset.i$val
            idx <- which(abs(v.dist_to_ASB.j) <= th.distance_to_ASB)
            
            if(b.val){
              v.vals.dist_to_ASB.j <- v.vals.dist_to_ASB.j[idx]
            }
            
            v.dist_to_ASB.j <- v.dist_to_ASB.j[idx]
            v.dist_to_ASB.i <- c(v.dist_to_ASB.i, v.dist_to_ASB.j)
            v.vals.dist_to_ASB.i <- c(v.vals.dist_to_ASB.i, v.vals.dist_to_ASB.j)
            
          }
          
        v.dist_to_ASB <- c(v.dist_to_ASB, v.dist_to_ASB.i)
        v.vals.dist_to_ASB <- c(v.vals.dist_to_ASB, v.vals.dist_to_ASB.i)
        
      }
      
      # + / - 1000
      # v.dist_to_ASB <- v.dist_to_ASB[abs(v.dist_to_ASB) <= th.distance_to_ASB]
      
      
      v.breaks_sets <- seq(-th.distance_to_ASB, (th.distance_to_ASB+1), n.binWidth)
      v.break_means <- numeric(length(v.breaks_sets))
      names(v.break_means) <- v.breaks_sets
      
      for(j in 1:(length(v.breaks_sets) - 1)){
        idx <- which(v.breaks_sets[j] <= v.dist_to_ASB & v.dist_to_ASB <= v.breaks_sets[j + 1])
        v.break_means[j] <- mean(v.vals.dist_to_ASB[idx])
      }
      
      v.break_means <- v.break_means[-length(v.break_means)]
      
      # maker larger bins - 20 
      # names(v.vals.dist_to_ASB) <- v.dist_to_ASB
      # v.break_means <- numeric(length(unique(names(v.vals.dist_to_ASB))))
      # names(v.break_means) <- unique(names(v.vals.dist_to_ASB))
      # 
      # for(j in 1:length(v.break_means)){
      #   idx.j <- which(names(v.vals.dist_to_ASB) == names(v.break_means)[j])
      #   v.break_means[j] = mean(v.vals.dist_to_ASB[idx.j])
      # }
      
      # plot(as.numeric(names(v.break_means)), v.break_means, col = "red")
      
      df.distPlot <- rbind(df.distPlot, data.frame(bin = as.numeric(names(v.break_means)), val = v.break_means, color = v.species[s], line_type = v.stringency[s]))
    }
    
    l.res_plots[[k]] <- df.distPlot
    rm(df.dataset)
    
  }
  names(l.res_plots) <- v.datasets 
  
  saveRDS(l.res_plots, paste("tmp/l.chromatin_opening_distance_plots_With_motifs_", timeStamp, ".rds", sep = ""))
  
  
  
  
  
  #l.res_plots <- readRDS(paste("tmp/l.chromatin_opening_distance_plots_With_motifs_", timeStamp, ".rds", sep = ""))
  
  i <- 1
  
  df.distPlot <- l.res_plots[[i]] 
  
  
  ggplot(df.distPlot, aes(x=bin, y=val, colour = color, linetype = line_type)) + stat_smooth(aes(x = bin, y = val), method = "lm",  formula = y ~ poly(x, 21), se = FALSE, n = 1000) + scale_colour_manual(values = c("green", "blue")) + theme_bw() + xlim(-950, 950)
  ggplot(df.distPlot, aes(x=bin, y=val, colour = color, linetype = line_type)) + stat_smooth(aes(x = bin, y = val), method = "lm",  formula = y ~ poly(x, 21), se = FALSE) + scale_colour_manual(values = c("green", "blue")) + theme_bw()
  
  
  # 6 x 4
  for(i in 1:length(l.res_plots)){
    df.distPlot <- l.res_plots[[k]]
    write.table(df.distPlot, paste("figures_paper/df.distance_plots_binned_",v.datasets[k], "_", timeStamp, "_with_motif_information.txt", sep = ""), sep = "\t", row.names = FALSE)
    p <- ggplot(l.res_plots[[i]], aes(x=bin, y=val, colour = color, linetype = line_type)) + stat_smooth(aes(x = bin, y = val), method = "lm",  formula = y ~ poly(x, 21), se = FALSE, n = 1000) + scale_colour_manual(values = c("blue", "green")) + theme_bw() + xlim(-950, 950)
    #ggplot(df.distPlot, aes(x=bin, y=val, colour = color, linetype = line_type)) + stat_smooth(aes(x = bin, y = val), method = "lm",  formula = y ~ poly(x, 21), se = FALSE, n = 1000) + scale_colour_manual(values = c("green", "blue")) + theme_bw() + xlim(-950, 950)
    pdf(paste("figures_paper/DistancePlots/", names(l.res_plots)[i],"_with_motif_information.pdf", sep = ""), width = 6, height = 4)
    print(p)
    dev.off()
  }
  
  


}



# Figure 1 F
AT_ZM_overlap = function(){
  
  df.ZM_Ath_overlaps <- read.csv("datasets_paper/ZmvsAth_1_to_1.csv")
  df.ZM_Ath_overlaps["p.val"] <- 1
  df.ZM_Ath_overlaps["fc"] <- 0
  
  # overla analysis 
  for(i in 1:nrow(df.ZM_Ath_overlaps)){
    
    ## domain genes in target genes
    hitInSample = n_A_B = df.ZM_Ath_overlaps$Overlap[i]
    sampleSize = n_A = df.ZM_Ath_overlaps$unique_Ath[i] +  hitInSample
    ## domain genes in genome
    hitInPop = n_B =  df.ZM_Ath_overlaps$unique_Zm[i] + hitInSample
    popSize = n_C =  df.ZM_Ath_overlaps$total[1]
    
    failInPop = n_C-n_B
    
    fc <- (n_A_B / n_A) / (n_B / n_C)
    pval = fisher.test(matrix(c(hitInSample, hitInPop-hitInSample, sampleSize-hitInSample, failInPop-sampleSize +hitInSample), 2, 2), alternative='greater')$p.value; 
    
    df.ZM_Ath_overlaps$fc[i] <- fc
    df.ZM_Ath_overlaps$p.val[i] <- pval
    
  }
  
  ####
  
  v.categories <- c("BR targets but not regulated", "BR targets and regulated", "no BR targets but regulated", "no BR targets and not regulated")
  
  m.heatmap.fc <- matrix(0, nrow = 4, ncol = 4, dimnames = list(v.categories, v.categories))
  m.heatmap.pval <- matrix(1, nrow = 4, ncol = 4, dimnames = list(v.categories, v.categories))
  
  m.heatmap.fc[1:16] <- df.ZM_Ath_overlaps$fc
  m.heatmap.fc <- t(m.heatmap.fc)
  
  m.heatmap.pval[1:16] <- df.ZM_Ath_overlaps$p.val
  m.heatmap.pval <- t(m.heatmap.pval)
  
  library(pheatmap)
  p <- pheatmap((m.heatmap.fc), cluster_rows = FALSE,cluster_cols = FALSE, labels_row = NULL, labels_col = NULL)
  
}






