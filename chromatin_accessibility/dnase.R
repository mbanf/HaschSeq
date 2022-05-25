asb_dnase_allelic_bias = function(l.SNPs=l.SNPs, th.bp_offset = 20,   n.cpus = 5){
  # figures: 2h 
  
  th.distance_to_ASB <-  2000
  n.chromosomes <- 10
  
  # n.binWidth = 20
  b.val = TRUE
  v.species = c("Mo17", "B73")
  
  # v.binWidth = 
  degrees = 15
  th.distance_to_ASB = 5000
  th.padding = 50
  width = 6
  height = 4
  
  group = "genic_and_nongenic"
  
  message("plot asb methylation values versus post frequency")
  
  df.ASBs <- l.SNPs[[1]]
  df.bpSNPs <- l.SNPs[[2]]
  
  l.SNP_selection = vector(mode = "list", length = 2)
  l.SNP_selection[[1]] = df.ASBs
  l.SNP_selection[[2]] = df.bpSNPs
  

  v.datasets <- c("DNase", "B73_CpG",  "MO17_CpG", "B73_CHG", "MO17_CHG", "B73_CHH", "MO17_CHH") 
  v.formats <- c("bedGraph", "bigWig", "bigWig", "bigWig", "bigWig", "bigWig", "bigWig")
  v.has_ranges = c("yes",  "no", "no",  "no", "no",  "no", "no")
  

  
  # l.results <- vector(mode = "list", length = length(v.datasets))
  
  l.SNPs_w_chromatin <- vector(mode = "list", length = 2)
  l.SNPs_w_chromatin[[1]] <- vector(mode = "list", length = length(v.datasets))
  l.SNPs_w_chromatin[[2]] <- vector(mode = "list", length = length(v.datasets))
  
  
  
  v.filename <- "../data/methylation/GSE94291_DNase_ist.bedGraph"
  df.dataset <- read.table(v.filename, header = FALSE, sep ="\t", quote = "", stringsAsFactors = FALSE)
  names(df.dataset) <- c("chr",  "start", "end", "val")
  df.dataset = subset(df.dataset, df.dataset$val > 0) # remove lack in coverage 
  
  df.SNPs <- l.SNPs[[1]][,c("contig", "position", "gene.ID", "POSTfreq")]
  
  df.SNPs_w_chromatin = c()
  
  strt<-Sys.time() 
  cl<-makeCluster(n.cpus)
  registerDoParallel(cl)
  
  l.SNPs_w_data <- foreach(i = 1:n.chromosomes, .packages=c("seqinr", "VariantAnnotation", "Biostrings")) %dopar% {   
    df.SNPs.i <- subset(df.SNPs, df.SNPs$contig == i)
    df.dataset.i = subset(df.dataset, df.dataset$chr == i)
    pb <- txtProgressBar(min = 0, max = nrow(df.SNPs.i), style = 3)
    for(j in 1:nrow(df.SNPs.i)){
      setTxtProgressBar(pb, j)
      idx = which(df.dataset.i$start <= df.SNPs.i$position[j] & df.SNPs.i$position[j] <= df.dataset.i$end)  
      if(length(idx) > 0){
        df.SNPs.i[j,"dnase"] = mean(df.dataset.i$val[idx])
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
  
  plot(df.SNPs_w_chromatin$POSTfreq, df.SNPs_w_chromatin$dnase)
  
  
  ## Data in a data.frame
  x1 <- df.SNPs_w_chromatin$POSTfreq
  x2 <- df.SNPs_w_chromatin$dnase
  scatter_plot(x1, x2, y_limit = c(0, 1.5))
  
  
  
  
  
    for(s in 1:length(l.SNP_selection)){
      
      message(k, " / ", s)
      
      df.SNPs = l.SNP_selection[[s]]
      
      df.SNPs_w_chromatin = c()
      
      strt<-Sys.time() 
      cl<-makeCluster(n.cpus)
      registerDoParallel(cl)
      
      l.SNPs_w_data <- foreach(i = 1:n.chromosomes, .packages=c("seqinr", "VariantAnnotation", "Biostrings")) %dopar% {   
        df.SNPs.i <- subset(df.SNPs, df.SNPs$contig == i)
        df.dataset.i = subset(df.dataset, df.dataset$chr == i)
        pb <- txtProgressBar(min = 0, max = nrow(df.SNPs.i), style = 3)
        for(j in 1:nrow(df.SNPs.i)){
          setTxtProgressBar(pb, j)
          if(v.has_ranges[k] == "yes"){
            idx = which(df.dataset.i$start <= df.SNPs.i$position[j] & df.SNPs.i$position[j] <= df.dataset.i$end)  
          }else{ # what is the closest 
            idx = which(abs(df.dataset.i$start-df.SNPs.i$position[j]) <= th.bp_offset) # min(abs(df.dataset$start-df.SNPs.i$position[j])))
          }
          if(length(idx) > 0){
            df.SNPs.i[j,v.datasets[k]] = mean(df.dataset.i$val[idx])
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
      
      # saveRDS(l.SNPs_w_chromatin, "l.SNPs_w_chromatin_ASBs_and_bgSNPs.rds")
      rm(l.SNPs_w_data)
      rm(df.SNPs_w_chromatin)
    }
    rm(df.dataset)
  }
  
  l.SNPs_w_chromatin = readRDS("l.SNPs_w_chromatin_ASBs_and_bgSNPs.rds")
  
  l.diff_data = vector(mode = "list", length = 3)
  l.diff_data[[1]] = l.SNPs_w_chromatin[[1]][[3]] - l.SNPs_w_chromatin[[1]][[2]]
  
  l.diff_SNPs_w_chromatin <- vector(mode = "list", length = 2)
  l.diff_SNPs_w_chromatin[[1]] <- vector(mode = "list", length = 4)
  l.diff_SNPs_w_chromatin[[2]] <- vector(mode = "list", length = 4)
  
  
  # estimate the number of genes... 
  l.diff_SNPs_w_chromatin[[1]][[1]] = l.SNPs_w_chromatin[[1]][[1]]$DNase
  l.diff_SNPs_w_chromatin[[2]][[1]] = l.SNPs_w_chromatin[[2]][[1]]$DNase
  
  l.diff_SNPs_w_chromatin[[1]][[2]] = l.SNPs_w_chromatin[[1]][[3]]$MO17_CpG - l.SNPs_w_chromatin[[1]][[2]]$B73_CpG
  l.diff_SNPs_w_chromatin[[2]][[2]] = l.SNPs_w_chromatin[[2]][[3]]$MO17_CpG - l.SNPs_w_chromatin[[2]][[2]]$B73_CpG
  
  l.diff_SNPs_w_chromatin[[1]][[3]] = l.SNPs_w_chromatin[[1]][[5]]$MO17_CHG - l.SNPs_w_chromatin[[1]][[4]]$B73_CHG
  l.diff_SNPs_w_chromatin[[2]][[3]] = l.SNPs_w_chromatin[[2]][[5]]$MO17_CHG - l.SNPs_w_chromatin[[2]][[4]]$B73_CHG
  
  l.diff_SNPs_w_chromatin[[1]][[4]] = l.SNPs_w_chromatin[[1]][[7]]$MO17_CHH - l.SNPs_w_chromatin[[1]][[6]]$B73_CHH
  l.diff_SNPs_w_chromatin[[2]][[4]] = l.SNPs_w_chromatin[[2]][[7]]$MO17_CHH - l.SNPs_w_chromatin[[2]][[6]]$B73_CHH
  
  
  # filter the background snps by asb associated genes
  # filter background snps 
  
  df.ASBs <- l.SNPs[[1]][,c("contig", "position", "gene.ID")]
  df.bpSNPs <- l.SNPs[[2]]
  
  v.gns_ASBs = unique(df.ASBs$gene.ID)
  idx.bgSNPs = which(!df.bpSNPs$gene.ID %in% v.gns_ASBs)
  # df.bpSNPs <- subset(df.bpSNPs, !df.bpSNPs$gene.ID %in% v.gns_ASBs) # remove bgSNPs associated with ASB linked genes ?
  
  # th = quantile(l.diff_SNPs_w_chromatin[[2]][[1]][idx.bgSNPs], 0.95, na.rm = T)
  # length(which( l.diff_SNPs_w_chromatin[[1]][[1]]> th))
  
  data = c("DNase  >= 0.95", "CpG (MO17 - B73) >= 0.95", "CHG  (MO17 - B73)  >= 0.95", "CHH  (MO17 - B73)  >= 0.95")
  
  for(i in 1:4){
    th = quantile(l.diff_SNPs_w_chromatin[[2]][[i]][idx.bgSNPs], 0.95, na.rm = T)
    idx = which(l.diff_SNPs_w_chromatin[[1]][[i]] >= th)
    
    df.ASBs[,data[i]] = "no"
    df.ASBs[idx,data[i]] = "yes"
  }
  
  
  data = c("CpG (MO17 - B73) <= 0.05", "CHG  (MO17 - B73)  <= 0.05", "CHH  (MO17 - B73)  <= 0.05")
  
  for(i in 2:4){
    th = quantile(l.diff_SNPs_w_chromatin[[2]][[i]][idx.bgSNPs], 0.05, na.rm = T)
    idx = which(l.diff_SNPs_w_chromatin[[1]][[i]] <= th)
    df.ASBs[,data[i-1]] = "no"
    df.ASBs[idx,data[i-1]] = "yes"
    
  }
  
  ### 
  idx = names(df.ASBs)[!names(df.ASBs) %in% c("contig", "position", "gene.ID")]
  res = apply(df.ASBs[,idx], 1, function(a) if("yes" %in% a){TRUE}else{FALSE})
  
  df.ASBs[,"is_methylated"] = res
  
  # saveRDS(df.ASBs, "df.ASBs_methylation.rds")
  # saveRDS(l.bgSNPs_w_chromatin, "l.bgSNPs_w_chromatin.rds")
  
  # 
  
  
  # saveRDS(l.results, "l.results.rds")
  l.ASB_SNPs_w_chromatin = readRDS("../tmp/l.SNPs_w_chromatin_ASBs.rds")
  l.bgSNPs_w_chromatin = readRDS("../tmp/l.bgSNPs_w_chromatin.rds")
  
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
  
  for(k in 1:length(v.datasets)){
    
    df = data.frame(val = c(l.results[[k]][[1]]$val, l.results[[k]][[2]]$val),
                    group = c(rep("ASBs", nrow(l.results[[k]][[1]])), rep("bg SNPs", nrow(l.results[[k]][[2]]))))
    
    ggplot(df, aes(x=val, fill=group)) + geom_density(alpha=0.4) + theme_bw()
    
    file = paste("/tmp/", v.datasets[k], ".pdf", sep = "")
    ggsave(file, width = 20, height = 20, units = "cm")
    
  }
  
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
  
  
  ###
  
  df.bQTL_gene_partitioning <- df.ASBs
  df.ASB_w_RNASeq <-  merge(df.bQTL_gene_partitioning, df.rnaseq.B73vsMo17, by.x = "gene.ID", by.y = "Gene_ID_AGPv4", all = FALSE)
  
  df.ASBs["RNASeq <= 0.05"] = NA
  df.ASBs["RNASeq (foldchange)"] = NA
  idx = which(df.ASBs$gene.ID %in% df.rnaseq.B73vsMo17$Gene_ID_AGPv4)
  for(i in 1:nrow(df.ASBs)){
    idx = which(df.rnaseq.B73vsMo17$Gene_ID_AGPv4 == df.ASBs$gene.ID[i])[1]
    if(length(idx) > 0 & !is.na(idx)){
      
      pval = t.test(c(df.rnaseq.B73vsMo17$B73_control_rep1_batch1[idx],df.rnaseq.B73vsMo17$B73_control_rep2_batch1[idx],df.rnaseq.B73vsMo17$B73_control_rep2_batch1[idx]),
                    c(df.rnaseq.B73vsMo17$Mo17_control_rep1_batch2[idx],df.rnaseq.B73vsMo17$Mo17_control_rep2_batch2[idx],df.rnaseq.B73vsMo17$Mo17_control_rep3_batch2[idx]))$p.val
      
      if(pval <= 0.05 & !is.na(pval)){
        df.ASBs$`RNASeq <= 0.05`[i] = "yes"
      }else{
        df.ASBs$`RNASeq <= 0.05`[i] = "no"
      } 
      
      df.ASBs$`RNASeq (foldchange)`[i] = df.rnaseq.B73vsMo17$B73_vs_Mo17[idx]
    }
  }
  
  data = c("DNase  >= 0.95", "CpG (MO17 - B73) >= 0.95", "CHG  (MO17 - B73)  >= 0.95", "CHH  (MO17 - B73)  >= 0.95", 
           "CpG (MO17 - B73) <= 0.05", "CHG  (MO17 - B73)  <= 0.05", "CHH  (MO17 - B73)  <= 0.05", "RNASeq <= 0.05")
  
  df.ASBs["any"] = ifelse(apply(df.ASBs[,data], 1, function(m) { any(m == "yes") }) == TRUE, "yes", "no")
  
  # t-test and fold change 
  write.table(df.ASBs, "df.ASBs_DNaseMethylationRNASeq.txt", row.names = F, sep = "\t")
  
}