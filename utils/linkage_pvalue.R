postTotal.significant <- l.selection[[1]] # 

# initial linkage disequilibrium filtering 
postTotal.significant.tmp <- postTotal.significant
postTotal.significant <- postTotal.significant[-(1:nrow(postTotal.significant)),]

n.flips <- 0

strt<-Sys.time() 
cl<-makeCluster(min(n.chromosomes, n.cpus))
registerDoParallel(cl)

l.postTotal.significant <- foreach(i = 1:n.chromosomes) %dopar% {     
  
  postTotal.significant.i <- subset(postTotal.significant.tmp, postTotal.significant.tmp$contig == i)
  postTotal.significant.i <- postTotal.significant.i[order(postTotal.significant.i$position),]
  
  finished <- FALSE
  df.keep <- postTotal.significant[-(1:nrow(postTotal.significant)),]
  
  while(!finished){
    
    # estimate distances of neighboring SNPs
    dist <- postTotal.significant.i$position[2:nrow(postTotal.significant.i)] - postTotal.significant.i$position[1:(nrow(postTotal.significant.i) - 1)]
    
    if(!any(dist < s.disequilibriumDistance)){
      finished <- TRUE
    }
    
    for(j in 1:length(dist)){
      
      if(dist[j] < s.disequilibriumDistance){
        
        # pairwise evaluation
        set <- postTotal.significant.i[c(j,(j+1)),]
        
        # using p-values to guide selection 
        #if(s == 1) # in case of significant
        i.max <- which(set$`p-value (corrected)` == max(set$`p-value (corrected)`))
        # else
        #   i.max <- which(set$`p-value (corrected)` == min(set$`p-value (corrected)`))
        # 
        # different directions - keep both
        if(set$POSTfreq[1] < 0.5 & set$POSTfreq[2] > 0.5 | set$POSTfreq[1] > 0.5 & set$POSTfreq[2] < 0.5){
          
          # keep both
          n.flips <- n.flips + 1
          df.keep <- rbind(df.keep, set)
          
        }
        
        # if s = 1 : remove SNP within proximity with higher p-value
        idx <- ifelse(i.max == 1, j, j + 1)
        postTotal.significant.i <- postTotal.significant.i[-idx,]
        
        break
        
      }
    }
  }
  
  postTotal.significant.i <- rbind(postTotal.significant.i, df.keep)
  postTotal.significant.i <- unique(postTotal.significant.i)
  
  # postTotal.significant <- rbind(postTotal.significant, postTotal.significant.i)
  postTotal.significant.i
  
}

stopCluster(cl)
print(Sys.time()-strt)

# filtered significant bQTLs
postTotal.significant <- c()
for(i in 1:n.chromosomes){
  postTotal.significant <- rbind(postTotal.significant, l.postTotal.significant[[i]])
}

l.selection[[1]] <- postTotal.significant

hist(l.selection[[1]]$POSTfreq, breaks = 100) # histogram of significant peaks

if(b.save_intermediate_results){
  write.table(l.selection[[1]], paste(folder_output, "/bQTL_after_input_blacklisting_pvalue_filter_in_peaks_after_LinkageDisequilibrium.txt", sep = ""),  row.names = FALSE, quote = FALSE, sep ="\t")
}