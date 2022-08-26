scatter_plot <- function(x, y, xlab = "allelic bias", ylab = "dnase level", x_limit = NULL, y_limit = NULL){
  
  df <- data.frame(x,y)
  
  ## Use densCols() output to get density at each point
  v <- densCols(x,y, colramp=colorRampPalette(c("black", "white")))
  df$dens <- col2rgb(v)[1,] + 1L
  
  ## Map densities to colors
  cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", 
                              "#FCFF00", "#FF9400", "#FF3100"))(256)
  df$col <- cols[df$dens]

  ## Plot it, reordering rows so that densest points are plotted on top
  plot(y~x, data=df[order(df$dens),], pch=20, col=col, cex=0.7, xlab = xlab, ylab = ylab, xlim=x_limit, ylim=y_limit)
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
methylation_occupancy = function(l.SNPs,
                                n.chromosomes = 10,
                                n.cpus = 1,
                                v.filenames = c("data/dnase/GSE94291_DNase_ist.bedGraph", 
                                                 "data/methylation/B73_CpG.bw", 
                                                 "data/methylation/MO17_CpG.bw",
                                                 
                                                 "data/methylation/B73_CHG.bw",
                                                 "data/methylation/MO17_CHG.bw",
                                                 
                                                 "data/methylation/B73_CHH.bw",
                                                 "data/methylation/MO17_CHH.bw"),
                                
                                v.datasets = c("DNase", "B73_CpG",  "MO17_CpG", "B73_CHG", "MO17_CHG", "B73_CHH", "MO17_CHH"),
                                v.formats = c("bedGraph", "bigWig", "bigWig", "bigWig", "bigWig", "bigWig", "bigWig"),
                                v.has_ranges = c("yes",  "no", "no",  "no", "no",  "no", "no")
                                ){
  
  message("(differential) methylation values versus post frequency")
  
  df.ASBs <- l.SNPs[[1]]
  v.gns_ASBs = unique(df.ASBs$gene.ID)
  df.bpSNPs <- l.SNPs[[2]]
  df.bpSNPs <- subset(df.bpSNPs, !df.bpSNPs$gene.ID %in% v.gns_ASBs) # 
  
  l.SNP_selection = vector(mode = "list", length = 2)
  l.SNP_selection[[1]] = df.ASBs
  l.SNP_selection[[2]] = df.bpSNPs
  
  for(s in 1:length(l.SNP_selection)){
    for(j in 1:length(v.datasets)){
      l.SNP_selection[[s]][v.datasets[j]] <- NA
    }
  }

  l.SNPs_w_chromatin <- vector(mode = "list", length = 2)
  l.SNPs_w_chromatin[[1]] <- vector(mode = "list", length = length(v.datasets))
  l.SNPs_w_chromatin[[2]] <- vector(mode = "list", length = length(v.datasets))
  
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
    
    for(s in 1:length(l.SNP_selection)){
      
      message(k, " / ", s)
      
      df.SNPs = l.SNP_selection[[s]]
      
      df.SNPs_w_chromatin = c()
      
      strt<-Sys.time() 
      cl<-makeCluster(n.cpus)
      registerDoParallel(cl)
      
      l.SNPs_w_data <- foreach(i = 1:n.chromosomes, .packages=c("seqinr", "VariantAnnotation", "Biostrings")) %dopar% {   
        
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
            df.SNPs.i[j,v.datasets[k]] = mean(df.dataset$val[idx])
          }
        }
        close(pb)
        df.SNPs.i
      }
      
      stopCluster(cl)
      print(Sys.time()-strt)
      
      df.SNPs_w_chromatin = c()
      for(i in 1:n.chromosomes){
        df.SNPs_w_chromatin <- rbind(df.SNPs_w_chromatin, l.SNPs_w_data[[i]])
      }
      
      l.SNPs_w_chromatin[[s]][[k]] = df.SNPs_w_chromatin
      
      rm(l.SNPs_w_data)
      rm(df.SNPs_w_chromatin)
    }
    
    rm(df.dataset)
  }
  
  saveRDS(l.SNPs_w_chromatin, paste(folder_tmp, "l.SNPs_w_chromatin_ASBs_and_bgSNPs.rds", sep = "/"))
  
}



plot_methylation_vs_allelic_bias = function(df.ASBs,
                                            l.SNPs_w_chromatin,
                                            
                                            th.quantile = 0.95,
                                            
                                            th.bp_offset = 20
){
  
  # Test plot background or ASBs
  rows = c("bp_window", "number_ASBs", "high_bias", "with_methylation_data", "have_methylation", "have_motifs", 
           "have_motif_and_methylation", "have_motif_and_no_methylation", 
           "have_no_motif_and_methylation", "have_no_motif_and_no_methylation",
           "have_methylation_and_high_bias", "have_motifs_and_high_bias", 
           "have_motif_and_methylation_and_high_bias", "have_motif_and_no_methylation_and_high_bias",
           "have_no_motif_and_methylation_and_high_bias", "have_no_motif_and_no_methylation_and_high_bias")
  
  df.methylation_summary <- data.frame(matrix(nrow = length(rows), ncol = 1))
  rownames(df.methylation_summary) = rows
  
  l.SNP_methylated <- readRDS(paste("l.ASBs_methylated.rds", sep ="")) # load!!
  df.ASBs[,"has_motif"] <- df.ASBs_w_motifs$has_motif # needs the motif ASBs
  
  df.SNPs <- l.SNPs[[1]]
  
  idx_hb = which(df.SNPs$POSTfreq >= 0.85 | df.SNPs$POSTfreq <= 0.15)
  high_bias = length(idx_hb)
  
  for(k in 1:length(offsets)){ # number of available ASBs, yes and no, subset by motifs... different offsets ... => create a dataframe
    
    idx1 = (l.SNP_methylated[[k]]$MO17_CpG >= 0.7 & l.SNP_methylated[[k]]$B73_CpG <= 0.3) 
    idx2 = (l.SNP_methylated[[k]]$MO17_CpG <= 0.3 & l.SNP_methylated[[k]]$B73_CpG >= 0.7)
    idx3 = (l.SNP_methylated[[k]]$MO17_CHG >= 0.7 & l.SNP_methylated[[k]]$B73_CHG <= 0.3)
    idx4 = (l.SNP_methylated[[k]]$MO17_CHG <= 0.3 & l.SNP_methylated[[k]]$B73_CHG >= 0.7)
    idx5 = (l.SNP_methylated[[k]]$MO17_CHH >= 0.7 & l.SNP_methylated[[k]]$B73_CHH <= 0.3)
    idx6 = (l.SNP_methylated[[k]]$MO17_CHH <= 0.3 & l.SNP_methylated[[k]]$B73_CHH >= 0.7)
    idx = idx1 | idx2 | idx3 | idx4 | idx4 | idx5 | idx6
    
    message("------------------")
    message("results for offset, ",  offsets[k], "bp distance methylation position to ASB")
    
    perc = round(length(idx[!is.na(idx)]) / length(idx) * 100, 1)
    
    w_methylation_data = length(idx[!is.na(idx)])
    message("methylation data available for ", perc, "% (", w_methylation_data, ") ASBs")
    
    # message("methylated (for available):")
    df.ASBs[,"is_methylated"] = idx
    have_methylation = as.numeric(table(idx)["TRUE"])
    have_motifs = as.numeric(table(df.ASBs$has_motif)["TRUE"])
    
    res = df.ASBs$has_motif & df.ASBs$is_methylated
    have_motif_and_methylation = as.numeric(table(res)["TRUE"])
    
    res = df.ASBs$has_motif & !df.ASBs$is_methylated
    have_motif_and_no_methylation = as.numeric(table(res)["TRUE"])
    
    res = !df.ASBs$has_motif & df.ASBs$is_methylated
    have_no_motif_and_methylation = as.numeric(table(res)["TRUE"])
    
    res = !df.ASBs$has_motif & !df.ASBs$is_methylated
    have_no_motif_and_no_methylation = as.numeric(table(res)["TRUE"])
    
    # check with high bias
    have_methylation_and_high_bias = as.numeric(table(idx[idx_hb])["TRUE"])
    have_motifs_and_high_bias = as.numeric(table(df.ASBs[idx_hb,]$has_motif)["TRUE"])
    
    res = df.ASBs[idx_hb,]$has_motif & df.ASBs[idx_hb,]$is_methylated
    have_motif_and_methylation_and_high_bias = as.numeric(table(res)["TRUE"])
    
    res = df.ASBs[idx_hb,]$has_motif & !df.ASBs[idx_hb,]$is_methylated
    have_motif_and_no_methylation_and_high_bias = as.numeric(table(res)["TRUE"])
    
    res = !df.ASBs[idx_hb,]$has_motif & df.ASBs[idx_hb,]$is_methylated
    have_no_motif_and_methylation_and_high_bias = as.numeric(table(res)["TRUE"])
    
    res = !df.ASBs[idx_hb,]$has_motif & !df.ASBs[idx_hb,]$is_methylated
    have_no_motif_and_no_methylation_and_high_bias = as.numeric(table(res)["TRUE"])
    
    vals <- c(offsets[k], nrow(df.ASBs), high_bias, w_methylation_data, have_methylation, have_motifs, 
              have_motif_and_methylation, have_motif_and_no_methylation, 
              have_no_motif_and_methylation, have_no_motif_and_no_methylation,
              have_methylation_and_high_bias, have_motifs_and_high_bias, 
              have_motif_and_methylation_and_high_bias, have_motif_and_no_methylation_and_high_bias, 
              have_no_motif_and_methylation_and_high_bias, have_no_motif_and_no_methylation_and_high_bias)
    
    df <- data.frame(vals)
    rownames(df) <- as.character(rows)
    
    df.methylation_summary <- cbind(df.methylation_summary, df)
    
    message("---------------")
  }
  
  df.methylation_summary = df.methylation_summary[,-1]
  colnames(df.methylation_summary) <- as.character(offsets)
  
  write.csv(df.methylation_summary, "df.methylation_summary.csv")
  
  # methylierung scatter plots?
  
  dev.off()
  k = 2
  title = paste("allelic biases vs methylation levels (averaged per SNP using ",  offsets[k], "bp window)", sep = "")
  par(mfrow=c(3,2))
  scatter_plot(df.SNPs$POSTfreq, l.SNP_methylated[[k]]$B73_CpG, xlab = "allelic bias", ylab = "methylation level (B73 CpG)")
  scatter_plot(df.SNPs$POSTfreq, l.SNP_methylated[[k]]$MO17_CpG, xlab = "allelic bias", ylab = "methylation level (MO17 CpG)")
  scatter_plot(df.SNPs$POSTfreq, l.SNP_methylated[[k]]$B73_CHG, xlab = "allelic bias", ylab = "methylation level (B73 CHG)")
  scatter_plot(df.SNPs$POSTfreq, l.SNP_methylated[[k]]$MO17_CHG, xlab = "allelic bias", ylab = "methylation level (MO17 CHG)")
  scatter_plot(df.SNPs$POSTfreq, l.SNP_methylated[[k]]$B73_CHH, xlab = "allelic bias", ylab = "methylation level (B73 CHH)")
  scatter_plot(df.SNPs$POSTfreq, l.SNP_methylated[[k]]$MO17_CHH, xlab = "allelic bias", ylab = "methylation level (MO17 CHH)")
  mtext(title, side = 3, line = -2, outer = TRUE, cex = 0.8)
  
  # 
  # repeat with displaying only methylation explained ASBs
  dev.off()
  title = paste("allelic biases vs methylation levels (averaged per SNP using ",  offsets[k], "bp window)", sep = "")
  par(mfrow=c(2,2))
  idx = which(l.SNP_methylated[[k]]$MO17_CpG <= 0.3 & l.SNP_methylated[[k]]$B73_CpG >= 0.7)
  scatter_plot(df.SNPs$POSTfreq[idx], l.SNP_methylated[[k]]$B73_CpG[idx], xlab = "allelic bias", ylab = "methylation level (B73 CpG)")
  
  idx = which(l.SNP_methylated[[k]]$MO17_CpG >= 0.7 & l.SNP_methylated[[k]]$B73_CpG <= 0.3)
  scatter_plot(df.SNPs$POSTfreq[idx], l.SNP_methylated[[k]]$MO17_CpG[idx], xlab = "allelic bias", ylab = "methylation level (MO17 CpG)")
  
  idx = which(l.SNP_methylated[[k]]$MO17_CHG <= 0.3 & l.SNP_methylated[[k]]$B73_CHG >= 0.7)
  scatter_plot(df.SNPs$POSTfreq[idx], l.SNP_methylated[[k]]$B73_CHG[idx], xlab = "allelic bias", ylab = "methylation level (B73 CHG)")
  
  idx = which(l.SNP_methylated[[k]]$MO17_CHG >= 0.7 & l.SNP_methylated[[k]]$B73_CHG <= 0.3)
  scatter_plot(df.SNPs$POSTfreq[idx], l.SNP_methylated[[k]]$MO17_CHG[idx], xlab = "allelic bias", ylab = "methylation level (MO17 CHG)")
  # 
  # idx = which(l.SNP_methylated[[k]]$MO17_CHH <= 0.3 & l.SNP_methylated[[k]]$B73_CHH >= 0.7)
  # scatter_plot(df.SNPs$POSTfreq[idx], l.SNP_methylated[[k]]$B73_CHH[idx], xlab = "allelic bias", ylab = "methylation level (B73 CHH)")
  # 
  # idx = which(l.SNP_methylated[[k]]$MO17_CHH >= 0.7 & l.SNP_methylated[[k]]$B73_CHH <= 0.3)
  # scatter_plot(df.SNPs$POSTfreq[idx], l.SNP_methylated[[k]]$MO17_CHH[idx], xlab = "allelic bias", ylab = "methylation level (MO17 CHH)")
  mtext(title, side = 3, line = -2, outer = TRUE, cex = 0.8)
  
  
  ### 
  
  dev.off()
  k = 2
  title = paste("allelic biases vs methylation levels (averaged per SNP using ",  offsets[k], "bp window)", sep = "")
  par(mfrow=c(3,1))
  
  idx1 = (l.SNP_methylated[[k]]$MO17_CpG <= 0.1 & l.SNP_methylated[[k]]$B73_CpG >= 0.7)
  idx2 = (l.SNP_methylated[[k]]$MO17_CpG >= 0.7 & l.SNP_methylated[[k]]$B73_CpG <= 0.1)
  idx = which(idx1 | idx2)
  col_vals = rep("black", nrow(df.SNPs)) 
  col_vals[idx] = "red"
  
  scatter_plot(df.SNPs$POSTfreq, l.SNP_methylated[[k]]$MO17_CpG - l.SNP_methylated[[k]]$B73_CpG, xlab = "allelic bias", ylab = "methylation level difference (CpG)", 
               x_limit = c(0,1), y_limit = c(-1,1), col_vals = col_vals)
  
  ### 
  idx1 = (l.SNP_methylated[[k]]$MO17_CHG <= 0.1 & l.SNP_methylated[[k]]$B73_CHG >= 0.7)
  idx2 = (l.SNP_methylated[[k]]$MO17_CHG >= 0.7 & l.SNP_methylated[[k]]$B73_CHG <= 0.1)
  idx = which(idx1 | idx2)
  col_vals2 = rep("black", nrow(df.SNPs)) 
  col_vals2[idx] = "red"
  
  scatter_plot(c(df.SNPs$POSTfreq, df.SNPs$POSTfreq), 
               c(l.SNP_methylated[[k]]$MO17_CpG - l.SNP_methylated[[k]]$B73_CpG,
                 l.SNP_methylated[[k]]$MO17_CHG - l.SNP_methylated[[k]]$B73_CHG), 
               xlab = "allelic bias", ylab = "methylation level difference - CpG (red), CHG (blue))", 
               x_limit = c(0,1), y_limit = c(-1,1), col_vals = c(col_vals, col_vals2))
  
  
  ###
  idx1 = (l.SNP_methylated[[k]]$MO17_CHG <= 0.1 & l.SNP_methylated[[k]]$B73_CHG >= 0.7)
  idx2 = (l.SNP_methylated[[k]]$MO17_CHG >= 0.7 & l.SNP_methylated[[k]]$B73_CHG <= 0.1)
  idx = which(idx1 | idx2)
  col_vals = rep("black", nrow(df.SNPs)) 
  col_vals[idx] = "red"
  
  scatter_plot(df.SNPs$POSTfreq, l.SNP_methylated[[k]]$MO17_CHG - l.SNP_methylated[[k]]$B73_CHG, xlab = "allelic bias", ylab = "methylation level difference (CHG)", 
               x_limit = c(0,1), y_limit = c(-1,1), col_vals = col_vals)
  
  ###
  
  idx1 = (l.SNP_methylated[[k]]$MO17_CHH <= 0.1 & l.SNP_methylated[[k]]$B73_CHH >= 0.7)
  idx2 = (l.SNP_methylated[[k]]$MO17_CHH >= 0.7 & l.SNP_methylated[[k]]$B73_CHH <= 0.1)
  idx = which(idx1 | idx2)
  col_vals = rep("black", nrow(df.SNPs)) 
  col_vals[idx] = "red"
  
  scatter_plot(df.SNPs$POSTfreq, l.SNP_methylated[[k]]$MO17_CHH - l.SNP_methylated[[k]]$B73_CHH, xlab = "allelic bias", ylab = "methylation level difference (CHH)", 
               x_limit = c(0,1), y_limit = c(-1,1), col_vals = col_vals)
  
  
  
  scatter_plot(df.SNPs$POSTfreq, l.SNP_methylated[[k]]$MO17_CHG - l.SNP_methylated[[k]]$B73_CHG, xlab = "allelic bias", ylab = "methylation level difference (CHG)")
  scatter_plot(df.SNPs$POSTfreq, l.SNP_methylated[[k]]$MO17_CHH - l.SNP_methylated[[k]]$B73_CHH, xlab = "allelic bias", ylab = "methylation level difference (CHH)")
  mtext(title, side = 3, line = -2, outer = TRUE, cex = 0.8)
  
  
}


                                                    
                                                    
 


differential_methylation_vs_allelic_bias = function(df.ASBs,
                                                    l.SNPs_w_chromatin,
                                                    
                                                    th.quantile = 0.95,
                                                    
                                                    th.bp_offset = 20,
                                                    th.distance_to_ASB = 2000,
                                                    n.chromosomes = 10,
                                                    n.cpus = 1,
                                                    b.val = TRUE,
                                                    v.species = c("Mo17", "B73"),
                                                    degrees = 15,
                                                    th.distance_to_ASB = 5000,
                                                    th.padding = 50,
                                                    width = 6,
                                                    height = 4,
                                                    group = "genic_and_nongenic",
                                                    folder_tmp = "",
){ 

  l.diff_data = vector(mode = "list", length = 3)
  l.diff_data[[1]] = l.SNPs_w_chromatin[[1]][[3]] - l.SNPs_w_chromatin[[1]][[2]]
  
  l.diff_SNPs_w_chromatin <- vector(mode = "list", length = 2)
  l.diff_SNPs_w_chromatin[[1]] <- vector(mode = "list", length = 4)
  l.diff_SNPs_w_chromatin[[2]] <- vector(mode = "list", length = 4)
  
  l.diff_SNPs_w_chromatin[[1]][[1]] = l.SNPs_w_chromatin[[1]][[1]]$DNase
  l.diff_SNPs_w_chromatin[[2]][[1]] = l.SNPs_w_chromatin[[2]][[1]]$DNase
  
  l.diff_SNPs_w_chromatin[[1]][[2]] = l.SNPs_w_chromatin[[1]][[3]]$MO17_CpG - l.SNPs_w_chromatin[[1]][[2]]$B73_CpG
  l.diff_SNPs_w_chromatin[[2]][[2]] = l.SNPs_w_chromatin[[2]][[3]]$MO17_CpG - l.SNPs_w_chromatin[[2]][[2]]$B73_CpG
  
  l.diff_SNPs_w_chromatin[[1]][[3]] = l.SNPs_w_chromatin[[1]][[5]]$MO17_CHG - l.SNPs_w_chromatin[[1]][[4]]$B73_CHG
  l.diff_SNPs_w_chromatin[[2]][[3]] = l.SNPs_w_chromatin[[2]][[5]]$MO17_CHG - l.SNPs_w_chromatin[[2]][[4]]$B73_CHG
  
  l.diff_SNPs_w_chromatin[[1]][[4]] = l.SNPs_w_chromatin[[1]][[7]]$MO17_CHH - l.SNPs_w_chromatin[[1]][[6]]$B73_CHH
  l.diff_SNPs_w_chromatin[[2]][[4]] = l.SNPs_w_chromatin[[2]][[7]]$MO17_CHH - l.SNPs_w_chromatin[[2]][[6]]$B73_CHH
  
  data = paste(c("DNase  >=", "CpG (MO17 - B73) >=", "CHG  (MO17 - B73)  >=", "CHH  (MO17 - B73)  >="), th.quantile)
  for(i in 1:4){
    th = quantile(l.diff_SNPs_w_chromatin[[2]][[i]], th.quantile, na.rm = T)
    idx = which(l.diff_SNPs_w_chromatin[[1]][[i]] >= th)
    df.ASBs[,data[i]] = "no"
    df.ASBs[idx,data[i]] = "yes"
  }
  
  data = paste(c("CpG (MO17 - B73) <=", "CHG  (MO17 - B73)  <=", "CHH  (MO17 - B73)  <="), 1 - th.quantile)
  for(i in 2:4){
    th = quantile(l.diff_SNPs_w_chromatin[[2]][[i]], 1 - th.quantile, na.rm = T)
    idx = which(l.diff_SNPs_w_chromatin[[1]][[i]] <= th)
    df.ASBs[,data[i-1]] = "no"
    df.ASBs[idx,data[i-1]] = "yes"
  }

  
  
  v.datasets <- c("DNase", "CpG", "CHG", "CHH")
  for(k in 1:length(v.datasets)){
    
    # differential methylation grouped density (ASBs, bgSNPs)
    df = data.frame(val = c(l.diff_SNPs_w_chromatin[[k]][[1]], l.diff_SNPs_w_chromatin[[k]][[2]]),
                    group = c(rep("ASBs", length(l.diff_SNPs_w_chromatin[[k]][[1]])), rep("bg SNPs", length(l.diff_SNPs_w_chromatin[[k]][[2]]))))
    
    ggplot(df, aes(x=val, fill=group)) + geom_density(alpha=0.4) + theme_bw()
    file = paste(folder_tmp, paste(v.datasets[k], ".pdf", sep = ""), sep = "/")
    ggsave(file, width = 20, height = 20, units = "cm")
    
    
    # differential methylation vs postfrequency (ABSs only)
    
    
   # scatter_plot(df.ASBs$POSTfreq, l.diff_SNPs_w_chromatin[[k]][[1]], y_limit = c(0, 1.5))
  }
  
  
  
  
  
  

  
  
  
  
  
  
  
  l.ASB_SNPs_w_chromatin = readRDS("l.SNPs_w_chromatin_ASBs.rds")
  l.bgSNPs_w_chromatin = readRDS("l.bgSNPs_w_chromatin.rds")
  
  diff_ASB_MO17_minus_B73 = l.ASB_SNPs_w_chromatin[[3]]$MO17_CpG - l.ASB_SNPs_w_chromatin[[2]]$B73_CpG
  diff_bgSNPs_MO17_minus_B73 = l.bgSNPs_w_chromatin[[3]]$MO17_CpG - l.bgSNPs_w_chromatin[[2]]$B73_CpG
  
  th.under = quantile(diff_bgSNPs_MO17_minus_B73, 0.05, na.rm = T)
  th.over = quantile(diff_bgSNPs_MO17_minus_B73, 0.95, na.rm = T)
  
  df.ASBs =  l.ASB_SNPs_w_chromatin[[1]]
  df.ASBs["delta_MO17_minus_B73"] = diff_ASB_MO17_minus_B73
  df.ASBs["significance_MO17_minus_B73"] = NA
  
  idx = which(df.ASBs$delta_MO17_minus_B73 >= th.over)
  df.ASBs$significance_MO17_minus_B73[idx] = "upper_95"
  
  idx = which(df.ASBs$delta_MO17_minus_B73 <= th.under)
  df.ASBs$significance_MO17_minus_B73[idx] = "lower_005"
  
  
  
  df.ASBs["DNase"] = l.ASB_SNPs_w_chromatin[[1]]$DNase
  df.ASBs["significance_DNase"] = NA
  
  th.over = quantile((l.bgSNPs_w_chromatin[[1]]$DNase), 0.95, na.rm = T)
  
  idx = which(df.ASBs$DNase >= th.over)
  df.ASBs$significance_DNase[idx] = "upper_95"
  
  
  
  hist(diff_ASB_MO17_minus_B73, 100)
  
  # for(k in 1:length(v.datasets)){
  #   write.csv2(l.results[[k]][[1]], paste("D:/HashSeq/tmp/", v.datasets[k], "_ASBs.csv", sep = ""), row.names = F)
  #   write.csv2(l.results[[k]][[2]], paste("D:/HashSeq/tmp/", v.datasets[k], "_bgSNPs.csv", sep = ""), row.names = F)
  # }
  
  
  
  
  
  ### 
  df.ASBs <- l.SNPs[[1]]
  for(j in 1:length(v.datasets)){
    df.ASBs[v.datasets[j]] <- 1
  }
  
  
  l.SNPs_w_chromatin = readRDS("l.SNPs_w_chromatin_ASBs.rds")
  
  
  for(k in 1:length(v.datasets)){
    
    df.bg_data = read.csv2(paste("D:/HashSeq/tmp/", v.datasets[k], "_bgSNPs.csv", sep = ""))
    
    th = quantile(df.bg_data$val, 0.95)
    
    l.SNPs_w_chromatin[[k]][">95BG"] = ifelse(l.SNPs_w_chromatin[[k]][,v.datasets[k]] >= th, "yes", "no")
    
    
    df.ASBs[,v.datasets[k]] =  ifelse(l.SNPs_w_chromatin[[k]][,v.datasets[k]] >= th, "yes", "no")
    # 
    
    # df.ASBs[v.datasets[k]] = 
  }
  
  write.csv2(df.ASBs, "A:/junkDNA.ai/df.ASBs.csv",row.names = F)
  
  
  table(l.SNPs_w_chromatin[[3]]$`>95BG`)
  
  
  ### with RNASEQ?
  
  # df.bQTL_gene_partitioning <- df.ASBs
  # df.ASB_w_RNASeq <-  merge(df.bQTL_gene_partitioning, df.rnaseq.B73vsMo17, by.x = "gene.ID", by.y = "Gene_ID_AGPv4", all = FALSE)
  # 
  # df.ASBs["RNASeq <= 0.05"] = NA
  # df.ASBs["RNASeq (foldchange)"] = NA
  # idx = which(df.ASBs$gene.ID %in% df.rnaseq.B73vsMo17$Gene_ID_AGPv4)
  # for(i in 1:nrow(df.ASBs)){
  #   idx = which(df.rnaseq.B73vsMo17$Gene_ID_AGPv4 == df.ASBs$gene.ID[i])[1]
  #   if(length(idx) > 0 & !is.na(idx)){
  #     
  #     pval = t.test(c(df.rnaseq.B73vsMo17$B73_control_rep1_batch1[idx],df.rnaseq.B73vsMo17$B73_control_rep2_batch1[idx],df.rnaseq.B73vsMo17$B73_control_rep2_batch1[idx]),
  #                   c(df.rnaseq.B73vsMo17$Mo17_control_rep1_batch2[idx],df.rnaseq.B73vsMo17$Mo17_control_rep2_batch2[idx],df.rnaseq.B73vsMo17$Mo17_control_rep3_batch2[idx]))$p.val
  #     
  #     if(pval <= 0.05 & !is.na(pval)){
  #       df.ASBs$`RNASeq <= 0.05`[i] = "yes"
  #     }else{
  #       df.ASBs$`RNASeq <= 0.05`[i] = "no"
  #     } 
  #     
  #     df.ASBs$`RNASeq (foldchange)`[i] = df.rnaseq.B73vsMo17$B73_vs_Mo17[idx]
  #   }
  # }
  # 
  # data = c("DNase  >= 0.95", "CpG (MO17 - B73) >= 0.95", "CHG  (MO17 - B73)  >= 0.95", "CHH  (MO17 - B73)  >= 0.95", 
  #          "CpG (MO17 - B73) <= 0.05", "CHG  (MO17 - B73)  <= 0.05", "CHH  (MO17 - B73)  <= 0.05", "RNASeq <= 0.05")
  # 
  # df.ASBs["any"] = ifelse(apply(df.ASBs[,data], 1, function(m) { any(m == "yes") }) == TRUE, "yes", "no")
  # write.table(df.ASBs, "df.ASBs_DNaseMethylationRNASeq.txt", row.names = F, sep = "\t")
  
}




methylation_frequency_vs_distance_to_asb <- function(df.ASBs, 
                                                    th.distance_to_ASB = 2000,
                                                    n.chromosomes = 10,
                                                    n.cpus = 4,
                                                    n.binWidth = 20,
                                                    b.val = TRUE,
                                                    v.species = c("Mo17", "B73"),
                                                    v.groups = c("all", "genic", "non_genic"),
                                  
                                                    v.species = c("Mo17", "Mo17" , "B73", "B73"),
                                                    v.stringency = c("< 0.5 | > 0.5", "< 0.15 | > 0.85", "< 0.5 | > 0.5", "< 0.15 | > 0.85"),
                                                    
                                                    v.filenames = c("data/methylation/B73_CHG.bw", 
                                                                     "data/methylation/B73_CHH.bw", 
                                                                     "data/methylation/B73_CpG.bw", 
                                                                     "data/methylation/MO17_CHG.bw", 
                                                                     "data/methylation/MO17_CHH.bw", 
                                                                     "data/methylation/MO17_CpG.bw"),
                                                    
                                                    v.datasets = c("B73_CHG", "B73_CHH", "B73_CpG", "MO17_CHG", "MO17_CHH", "MO17_CpG"),
                                                    v.formats = c("bigWig", "bigWig", "bigWig", "bigWig", "bigWig", "bigWig"),
                                                    
                                                    folder_tmp = "tmp",
                                                    folder_output = "results",
                                                    
                                                    do.plot = F){
  
  # Genomic feature profiling of ASBs (Methylation)
  
  # Methylation levels for CG, CHG, and CHH for B73 and Mo17 were extracted from Regulski et al. 2013
  # Methylation frequency versus distance (up to +/- 2 kbp) around each ASB were averaged over 20bp bins, 
  # and visualized by regions bound by BZR1 with
  # either high or low affinity levels depending on the inbred line. For B73, high and low
  # affinity bound regions were defined by a post frequency of >= 0.85 or <= 0.15, respectively
  # and oppositely for Mo17 by a post frequency <= 0.15 and >= 0.85), respectively.
  
  names(v.species) <- v.species 
  names(v.stringency) <- v.species
  
  for(g in 1:length(v.groups)){
      
      group = v.groups[g]
      df.SNPs <- df.ASBs 
      
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

          }
          stopCluster(cl)
          print(Sys.time()-strt)
        }
        
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
      
      # saveRDS(l.res_dists, paste(folder_tmp, paste("l.res_dists_methylation_",v.datasets[k], "_", group,"_", n.binWidth, "_", timeStamp, ".rds", sep = ""), sep ="/"))
      
      for(k in 1:length(l.res_dists)){
        write.table(l.res_dists[[k]], paste(folder_output, paste("df.methylation_frequency_vs_distance_to_asb_",v.datasets[k], "_", group,"_", n.binWidth,"_", timeStamp, ".txt", sep = ""), sep ="/"), sep = "\t", row.names = FALSE)
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
          
          df.distPlot <- rbind(df.distPlot, data.frame(bin = as.numeric(names(v.mean_per_break)), val = v.mean_per_break, 
                                                       color = v.species[s], line_type = v.stringency[s], stringsAsFactors = FALSE))
        }
        l.res_plots[[k]] <- df.distPlot
      }
      
      names(l.res_plots) <- v.datasets 
      # saveRDS(l.res_plots, paste(folder_tmp, paste("l.novel_methylation_datasets_distance_plots_",v.datasets[k], "_", group,"_", n.binWidth , "_", timeStamp, ".rds", sep = ""), sep ="/"))
      

      for(k in 1:length(l.res_plots)){
        write.table(l.res_plots[[k]], paste(folder_output, paste("df.methylation_frequency_vs_distance_to_asb_binned_",v.datasets[k], "_", group,"_", n.binWidth, "_", timeStamp, ".txt", sep = ""), sep ="/"), sep = "\t", row.names = FALSE)
      }
      
      # 6 x 4
      if(do.plot){
        # Fig. 4 i)-k) Average CpG (h), CHG (i) and CHH(j) methylation frequency in B73 (green) and Mo17 (blue) over ASB loci with a least 85% binding
        # bias towards B73 or Mo17. High affinity ( ___ ) and low affinity ( …. ) alleles are considered separately for each genotype.
        for(i in 1:length(l.res_plots)){
            p <- ggplot(l.res_plots[[i]], aes(x=bin, y=val, colour = color, linetype = line_type)) + stat_smooth(aes(x = bin, y = val), method = "lm",  formula = y ~ poly(x, 21), se = FALSE, n = 1000) + scale_colour_manual(values = c("green", "blue")) + theme_bw() + xlim(-950, 950)
            pdf(paste(folder_output, paste("Novel_Methylation/", names(l.res_plots)[i],".pdf", sep = ""), sep ="/"), width = 6, height = 4)
            print(p)
            dev.off()
          }
      }
      
    }
  
  
}

