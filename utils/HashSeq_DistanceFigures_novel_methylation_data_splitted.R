

# do a for loop 

for (thing in ls()) {
  message(thing)
  print(object.size(get(thing)), units='auto')
}



# parameter sets ... read from csv 


th.distance_to_ASB <-  2000
n.chromosomes <- 10

n.cpus <- 4
v.binWidth = c(20,40,60,75,100)
b.val = TRUE
v.species = c("Mo17", "B73")
v.groups <- c("all", "genic", "non_genic")

for(b in 1:length(v.binWidth)){
  n.binWidth = v.binWidth[b]
  for(g in 1:length(v.groups)){
    
        group = v.groups[g]
       
        df.SNPs <-  l.bQTL_gene_partitioning[[1]]
        # df.SNPs <- readRDS(paste("tmp/l.bQTL_gene_partitioning_withGeneDistances_backgroundSampled_", timeStamp, ".rds", sep = ""))[[1]]
        
        if(group == "genic"){
          df.SNPs <- subset(df.SNPs, df.SNPs$non_genic == "no")
        }else if(group == "non_genic"){
          df.SNPs <- subset(df.SNPs, df.SNPs$non_genic == "yes") 
        }
        
        # separate by directionality - only for ASBs
        l.df.SNPs.directionality <- vector(mode = "list", length = 4)
        l.df.SNPs.directionality[[1]] <- subset(df.SNPs, df.SNPs$POSTfreq <= 0.5)
        l.df.SNPs.directionality[[2]]  <- subset(df.SNPs, df.SNPs$POSTfreq <= 0.15)
        l.df.SNPs.directionality[[3]] <- subset(df.SNPs, df.SNPs$POSTfreq >= 0.5)
        l.df.SNPs.directionality[[4]]  <- subset(df.SNPs, df.SNPs$POSTfreq >= 0.85)
        
        # v.species = c("Mo17 < 0.5", "Mo17 < 0.15", "B73 > 0.5", "B73 > 0.85")
        
        v.species <- c("Mo17", "Mo17" , "B73", "B73")
        v.stringency <- c("< 0.5 | > 0.5", "< 0.15 | > 0.85", "< 0.5 | > 0.5", "< 0.15 | > 0.85")
        
        names(v.species) <- v.species 
        names(v.stringency) <- v.species
        
        v.filenames <- c("data/methylation/B73_CHG.bw", 
                         "data/methylation/B73_CHH.bw", 
                         "data/methylation/B73_CpG.bw", 
                         "data/methylation/MO17_CHG.bw", 
                         "data/methylation/MO17_CHH.bw", 
                         "data/methylation/MO17_CpG.bw")
        
        
        v.datasets <- c("B73_CHG", "B73_CHH", "B73_CpG", "MO17_CHG", "MO17_CHH", "MO17_CpG")
        
        v.formats <- c("bigWig", "bigWig", "bigWig", "bigWig", "bigWig", "bigWig")
        
        l.res_dists <- vector(mode = "list", length = length(v.datasets))
        
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
        
        
          # using the B73 
          if(FALSE){
            # distance plot to ASBs, binding peaks, etc. 
            df.distPlot <- c()
            strt<-Sys.time() 
            cl<-makeCluster(min(n.chromosomes, n.cpus))
            registerDoParallel(cl)
            l.res <- foreach(i = 1:n.chromosomes) %dopar% { 
            # for(i in 1:n.chromosomes){
              v.dist_to_ASB <- c()
              v.vals.dist_to_ASB <- c()
              
              df.SNPs.i <- subset(df.SNPs, df.SNPs$contig == i)
              df.distanceDataset.i <- subset(df.distanceDataset, df.distanceDataset$chr == i)
              
              v.dist_to_ASB.i <- c()
              v.vals.dist_to_ASB.i <- c()
              
              for(j in 1:nrow(df.SNPs.i)){
                
                v.dist_to_ASB.j <- (df.distanceDataset.i$start - df.SNPs.i$position[j])
                v.vals.dist_to_ASB.j <- df.distanceDataset.i$val
                
                idx <- which(abs(v.dist_to_ASB.j) <= th.distance_to_ASB)
                
                if(b.val){
                  v.vals.dist_to_ASB.j <- v.vals.dist_to_ASB.j[idx]
                }
                v.dist_to_ASB.j <- v.dist_to_ASB.j[idx]
                
                v.dist_to_ASB.i <- c(v.dist_to_ASB.i, v.dist_to_ASB.j)
                v.vals.dist_to_ASB.i <- c(v.vals.dist_to_ASB.i, v.vals.dist_to_ASB.j)
              }
              
              list(v.dist_to_ASB = v.dist_to_ASB, v.vals.dist_to_ASB=v.vals.dist_to_ASB)
              # v.dist_to_ASB <- c(v.dist_to_ASB, v.dist_to_ASB.i)
              # v.vals.dist_to_ASB <- c(v.vals.dist_to_ASB, v.vals.dist_to_ASB.i)
            }
            stopCluster(cl)
            print(Sys.time()-strt)
          }
        
          
          # + / - 1000
          # v.dist_to_ASB <- v.dist_to_ASB[abs(v.dist_to_ASB) <= th.distance_to_ASB]
          # 
          # names(v.vals.dist_to_ASB) <- v.dist_to_ASB
          # v.break_means <- numeric(length(unique(names(v.vals.dist_to_ASB))))
          # names(v.break_means) <- unique(names(v.vals.dist_to_ASB))
          # 
          # for(j in 1:length(v.break_means)){
          #   idx.j <- which(names(v.vals.dist_to_ASB) == names(v.break_means)[j])
          #   v.break_means[j] = mean(v.vals.dist_to_ASB[idx.j])aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
          # }
          # 
          #   # plot(as.numeric(names(v.break_means)), v.break_means, col = "red")
          #   
          #   df.distPlot <- rbind(df.distPlot, data.frame(bin = as.numeric(names(v.break_means)), val = v.break_means, species = v.species[s]))
          #   
          # }
          # 
          # ggplot(df.distPlot, aes(x=bin, y=val, col = species)) + stat_smooth(aes(x = bin, y = val), method = "lm",  formula = y ~ poly(x, 21), se = FALSE) + theme_bw()
        
        
          # v.colors_species <- c("green", "green", "blue", "blue")
          # names(v.colors_species) <- v.species
          # 
          # v.linetype_species <- c("solid", "dashed", "solid", "dashed")
          # names(v.linetype_species) <- v.species
          # 
          #####
          
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
            df.distSet <- rbind(df.distSet, data.frame(dist = as.numeric(v.dist_to_ASB), val = as.numeric(v.vals.dist_to_ASB), color = v.species[s], line_type = v.stringency[s]))
            
          }
          
          l.res_dists[[k]] <- df.distSet
          
          
          rm(df.dataset)
        
        }
         
        
        saveRDS(l.res_dists, paste("tmp/l.res_dists_novel_methylation_",v.datasets[k], "_", group,"_", n.binWidth, "_", timeStamp, ".rds", sep = ""))
        #saveRDS(l.res_dists, "tmp/l.res_dists_novel_methylation_datasets.rds")
        
        
        
        # l.res_dists <- readRDS("tmp/l.res_dists_novel_methylation_datasets.rds")
        
        
        
        # store datasets 
        for(k in 1:length(v.datasets)){
        
          df.distSet <- l.res_dists[[k]]
          # write.csv(df.distSet, paste("output/df.distance_plots_",v.datasets[k], "_", timeStamp, ".csv", sep = ""))
          
         # write.table(df.distSet, paste("output/methylation/df.distance_plots_",v.datasets[k], "_", timeStamp, ".txt", sep = ""), sep = "\t", row.names = FALSE)
          write.table(df.distSet, paste("figures_paper/DistancePlots/df.distance_plots_novel_methylation_",v.datasets[k], "_", group,"_", n.binWidth,"_", timeStamp, ".txt", sep = ""), sep = "\t", row.names = FALSE) 
        }
        
        
        
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
        #saveRDS(l.res_plots, paste("tmp/l.novel_methylation_datasets_distance_plots_", timeStamp, ".rds", sep = ""))
        saveRDS(l.res_plots, paste("tmp/l.novel_methylation_datasets_distance_plots_",v.datasets[k], "_", group,"_", n.binWidth , "_", timeStamp, ".rds", sep = ""))
        
        
        for(k in 1:length(v.datasets)){
          
          df.distPlot <- l.res_plots[[k]]
          #write.table(df.distPlot, paste("output/df.novel_methylation_datasets_distance_plots_binned_",v.datasets[k], "_", timeStamp, ".txt", sep = ""), sep = "\t", row.names = FALSE)
          write.table(df.distPlot, paste("figures_paper/DistancePlots/df.novel_methylation_datasets_distance_plots_binned_",v.datasets[k], "_", group,"_", n.binWidth, "_", timeStamp, ".txt", sep = ""), sep = "\t", row.names = FALSE)
          
        }
        
        
        
        
        # i <- 1
        
        #df.distPlot <- l.res_plots[[i]] 
        
        
        # 6 x 4
        
        # l.res_plots <- readRDS(paste("tmp/l.chromatin_opening_distance_plots_", timeStamp, ".rds", sep = ""))
        
        # ggplot(df.distPlot, aes(x=bin, y=val, colour = color, linetype = line_type)) + geom_line()
        
        #for(i in 1:length(l.res_plots)){
          
        
         # p <- ggplot(l.res_plots[[i]], aes(x=bin, y=val, colour = color, linetype = line_type)) + stat_smooth(aes(x = bin, y = val), method = "lm",  formula = y ~ poly(x, 21), se = FALSE, n = 1000) + scale_colour_manual(values = c("green", "blue")) + theme_bw() + xlim(-950, 950)
          
          #ggplot(df.distPlot, aes(x=bin, y=val, colour = color, linetype = line_type)) + stat_smooth(aes(x = bin, y = val), method = "lm",  formula = y ~ poly(x, 21), se = FALSE, n = 1000) + scale_colour_manual(values = c("green", "blue")) + theme_bw() + xlim(-950, 950)
          
          #pdf(paste("figures_paper/DistancePlots/Novel_Methylation/", names(l.res_plots)[i],".pdf", sep = ""), width = 6, height = 4)
          #print(p)
          #dev.off()
          
        #}
          
  }
}

# 
# + scale_linetype_manual(values = )
# 
# ggplot(df.distPlot, aes(x=bin, y=val,  fill = directionality, linetype = directionality)) + stat_smooth(aes(x = bin, y = val), method = "lm",  formula = y ~ poly(x, 21), se = FALSE) + theme_bw() 
# 
# 
# + scale_colour_manual(values = v.colors_species) + scale_linetype_manual(values = )
# 
# 
# 
# 
# 
# ggplot(mort3, aes(x = year, y = BCmort, col = State, linetype = State)) +
#   geom_line(lwd = 1) +
#   scale_linetype_manual(values = c(rep("solid", 10), rep("dashed", 6))) +
#   scale_color_manual(values = c(brewer.pal(10, "Set3"), brewer.pal(6, "Set3"))) +
#   opts(title = "BC mortality") +
#   theme_bw()
# 




# ggplot2 figures 


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




# gene expression directionality heatmap

