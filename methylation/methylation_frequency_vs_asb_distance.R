methylation_frequency_vs_asb_distance = function(l.bQTL_gene_partitioning){
  
  # Figures: 2 i,j,k
  
  library(ggplot2)
  
  n.binWidth = 20
  degrees = 15 #21
  th.distance_to_ASB = 2000
  th.padding = 50
  width = 6
  height = 4
  group = "all"
  timeStamp = 117
  
  # l.bQTL_gene_partitioning <- readRDS(paste("tmp/l.bQTL_gene_partitioning_withGeneDistances_backgroundSampled_", timeStamp, ".rds", sep = ""))
  df.SNPs <-  l.bQTL_gene_partitioning[[1]]
  
  
  # separate by directionality - only for ASBs
  l.df.SNPs.directionality <- vector(mode = "list", length = 4)
  l.df.SNPs.directionality[[1]] <- subset(df.SNPs, df.SNPs$POSTfreq <= 0.5)
  l.df.SNPs.directionality[[2]]  <- subset(df.SNPs, df.SNPs$POSTfreq <= 0.15)
  l.df.SNPs.directionality[[3]] <- subset(df.SNPs, df.SNPs$POSTfreq >= 0.5)
  l.df.SNPs.directionality[[4]]  <- subset(df.SNPs, df.SNPs$POSTfreq >= 0.85)
  
  # v.species = c("Mo17 < 0.5", "Mo17 < 0.15", "B73 > 0.5", "B73 > 0.85")
  
  v.species <- c("Mo17", "Mo17" , "B73", "B73")
  v.stringency <- c("< 0.5 | > 0.5", "< 0.15 | > 0.85", "< 0.5 | > 0.5", "< 0.15 | > 0.85")
  
  v.datasets <- c("H3K9", "DNase", "MartBS123_CHG", "MartBS123_CG","B73_CHG", "B73_CHH", "B73_CpG", "MO17_CHG", "MO17_CHH", "MO17_CpG")
  l.res_dists <- c(readRDS(paste("tmp/l.res_dists_general_", group,"_", timeStamp, ".rds", sep = "")), readRDS("tmp/l.res_dists_novel_methylation_datasets.rds"))
  
  l.res_plots <- vector(mode = "list", length = length(v.datasets))
  for(k in 1:length(v.datasets)){
    df.distPlot <- c()
    df.distSet <- l.res_dists[[k]]
    for(s in 1:length(l.df.SNPs.directionality)){
      
      df.SNPs <- l.df.SNPs.directionality[[s]]
      df.distSet.s <- subset(df.distSet, df.distSet$color == v.species[s] & df.distSet$line_type == v.stringency[s])
      
      v.dist_to_ASB <- df.distSet.s$dist
      v.vals.dist_to_ASB <- df.distSet.s$val
      # df.distPlot <- rbind(df.distPlot, data.frame(bin = as.numeric(names(v.mean_per_break)), val = v.mean_per_break, color = v.species[s], line_type = v.stringency[s]))
      
      v.breaks_sets <- seq(-th.distance_to_ASB, th.distance_to_ASB, n.binWidth)
      v.mean_per_break <- numeric(length(v.breaks_sets))
      names(v.mean_per_break) <- v.breaks_sets
      
      for(j in 1:(length(v.breaks_sets) - 1)){
        idx <- which(v.breaks_sets[j] <= v.dist_to_ASB & v.dist_to_ASB <= v.breaks_sets[j + 1])
        v.mean_per_break[j] <- mean(v.vals.dist_to_ASB[idx])
      }
      
      v.mean_per_break <- v.mean_per_break[-length(v.mean_per_break)]
      
      # names(v.vals.dist_to_ASB) <- v.dist_to_ASB
      # v.break_means <- numeric(length(unique(names(v.vals.dist_to_ASB))))
      # names(v.break_means) <- unique(names(v.vals.dist_to_ASB))
      # 
      # for(j in 1:length(v.break_means)){
      #   idx.j <- which(names(v.vals.dist_to_ASB) == names(v.break_means)[j])
      #   v.break_means[j] = mean(v.vals.dist_to_ASB[idx.j])
      # }
      # plot(as.numeric(names(v.break_means)), v.break_means, col = "red")
      
      df.distPlot <- rbind(df.distPlot, data.frame(bin = as.numeric(names(v.mean_per_break)), val = v.mean_per_break, color = v.species[s], line_type = v.stringency[s], stringsAsFactors = FALSE))
      
    }
    l.res_plots[[k]] <- df.distPlot
  }
  names(l.res_plots) <- v.datasets 
  saveRDS(l.res_plots, paste("../tmp/l.chromatin_opening_distance_plots_", group,"_", "binwidth_",n.binWidth, "_", timeStamp, ".rds", sep = ""))
  
  
  l.res_plots = readRDS(paste("../tmp/l.chromatin_opening_distance_plots_", group,"_", "binwidth_",n.binWidth, "_", timeStamp, ".rds", sep = ""))   
  for(k in 1:length(v.datasets)){
    df.distPlot <- l.res_plots[[k]]
    write.table(df.distPlot, paste("DistancePlots/df.distance_plots_binned_",v.datasets[k], "_", group,"_", "binwidth_", n.binWidth, "_", timeStamp, ".txt", sep = ""), sep = "\t", row.names = FALSE)
  }
  
  for(k in 1:length(v.datasets)){
    p <- ggplot(l.res_plots[[k]], aes(x=bin, y=val, colour = color, linetype = line_type)) + stat_smooth(aes(x = bin, y = val), method = "lm",  formula = y ~ poly(x, degree= degrees), se = FALSE, n = th.distance_to_ASB) + scale_colour_manual(values = c("green", "blue")) + theme_bw() + xlim(-th.distance_to_ASB + th.padding , th.distance_to_ASB - th.padding)
    #ggplot(df.distPlot, aes(x=bin, y=val, colour = color, linetype = line_type)) + stat_smooth(aes(x = bin, y = val), method = "lm",  formula = y ~ poly(x, 21), se = FALSE, n = 1000) + scale_colour_manual(values = c("green", "blue")) + theme_bw() + xlim(-950, 950)
    pdf(paste("figures_paper/DistancePlots/df.distance_plots_binned_",v.datasets[k], "_", group,"_", "binwidth_",n.binWidth, "_",  timeStamp,".pdf", sep = ""), width = width, height = height)
    print(p)
    dev.off()
  }
  
}