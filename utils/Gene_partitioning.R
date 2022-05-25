


# removes duplicates in ASBs
create_background_distribution_with_similar_gene_partitioning <- function(l.bQTL_gene_partitioning = l.bQTL_gene_partitioning, multiplyer = 10, v.partitions = v.partitions, b.duplicateRemoval = FALSE, seed.randomGenerator = 1234){
  
  set.seed(seed.randomGenerator)
  
  df.bQTL_gene_partitioning = l.bQTL_gene_partitioning[[1]]
  
  v.partitions.duplicates <- v.partitions[c(1,2,8)]
  
  if(b.duplicateRemoval){
    # identify duplicate partitions and remove => enter into general pipeline
    for(i in 1:2){
      i.set <- which(apply(l.bQTL_gene_partitioning[[i]][,v.partitions.duplicates], 1, function(m) length(which(m == "yes"))) == 2)
      i.set <- (!1:nrow(l.bQTL_gene_partitioning[[i]]) %in% i.set)
      l.bQTL_gene_partitioning[[i]] <- l.bQTL_gene_partitioning[[i]][i.set, ]
    }
  }
  # v.distribution <- apply(l.bQTL_gene_partitioning[[1]][,v.partitions.duplicates], 1, table)["yes",]
  df.bQTL_gene_partitioning_bg_equivalent <- c()
  
  for(i in 1:n.chromosomes){
    
    df.bQTL_gene_partitioning.asb.i <- subset(l.bQTL_gene_partitioning[[1]],  l.bQTL_gene_partitioning[[1]]$contig == i)
    df.bQTL_gene_partitioning_bg.i  <- subset(l.bQTL_gene_partitioning[[2]],  l.bQTL_gene_partitioning[[2]]$contig == i)
    # df.bg_snps.i <- subset(l.postTotal[[2]], l.postTotal[[2]]$contig == i)
    
    # <- df.bQTL_gene_partitioning.i[,v.partitions
    v.distribution <- apply(df.bQTL_gene_partitioning.asb.i[,v.partitions], 2, table)["yes",]
    v.distribution <- v.distribution * multiplyer
    
    df.bQTL_gene_partitioning_bg_equivalent.i <- c()
    
    for(j in 1:length(v.partitions)){
      if(v.partitions[j] != "gene"){
        df.bQTL_gene_partitioning_bg.i.set <- df.bQTL_gene_partitioning_bg.i[which(as.character(df.bQTL_gene_partitioning_bg.i[,v.partitions[j]]) == "yes"),]
        # df.bQTL_gene_partitioning_bg.i.set <- df.bQTL_gene_partitioning_bg.i.set[which(apply(df.bQTL_gene_partitioning_bg.i.set[,v.partitions[-j]], 1, function(m) all(as.character(m) == "no")) == TRUE),]
        i.set <- sample(nrow(df.bQTL_gene_partitioning_bg.i.set), v.distribution[j])
        df.bQTL_gene_partitioning_bg_equivalent.i <- rbind(df.bQTL_gene_partitioning_bg_equivalent.i, df.bQTL_gene_partitioning_bg.i.set[i.set,])
      }
    }
    
    df.bQTL_gene_partitioning_bg_equivalent <- rbind(df.bQTL_gene_partitioning_bg_equivalent, df.bQTL_gene_partitioning_bg_equivalent.i)
    
  }
  
  v.distribution <- apply(l.bQTL_gene_partitioning[[1]][,v.partitions], 2, table)["yes",]
  v.distribution.bg_equivalent <- apply(df.bQTL_gene_partitioning_bg_equivalent[,v.partitions], 2, table)["yes",]
  #df.bQTL_gene_partitioning_bg_equivalent <- l.bQTL_gene_partitioning[[2]]
  # print(v.distribution)
  # print(v.distribution.bg_equivalent)
  # print(v.distribution / v.distribution.bg_equivalent)
  # 
  df.distribution <- data.frame(ASBs = v.distribution, bgSNPs = v.distribution.bg_equivalent, multiplyer = (v.distribution.bg_equivalent / v.distribution))
  print(df.distribution)
  
  write.table(df.distribution, "df.distribution.txt", sep ="\t")
  
  l.bQTL_gene_partitioning[[2]] <- df.bQTL_gene_partitioning_bg_equivalent
  
  
  return(list(l.bQTL_gene_partitioning = l.bQTL_gene_partitioning, df.distribution = df.distribution))
  
}





perform_SNP_to_gene_partitioning <- function(df.gene_annotation=df.gene_annotation,
                                             v.genePartitions=v.genePartitions
                                             ){
  
  message("perform SNP partitioning")
  
  message("bQTL to gene partitioning - in eQTL, GWAS, RNASeq")
  
  l.bQTL_gene_partitioning <- vector(mode = "list", length = 2)
  df.partitions <- data.frame(partition = c("promoter_5kb", "promoter_1kb", "gene", "five_prime_UTR", "exon", "intron", "three_prime_UTR", "post_gene_1kb",  "non_genic", "total_bQTL"),  # "post_gene_5kb",
                              significant = numeric(10), 
                              non_significant = numeric(10))
  
  if(!b.loadGenePartitioning){
    
    for(s in 1:2){
      
      message("identifying gene partitions - positions of ASBs with ")
      
      if(TRUE){
        if(s == 1){
          postTotal.significant <- l.selection[[s]] # select the significant subset 
        }else{
          postTotal.significant <- l.postTotal[[s]]  
        }
      }else{
        postTotal.significant <- l.selection[[s]] # select the significant subset 
      }
      
      
      df.bQTL_gene_partitioning <- c()
      
      strt<-Sys.time() 
      cl<-makeCluster(min(n.chromosomes, n.cpus))
      registerDoParallel(cl)
      
      l.bQTL_gene_partitioning_per_chromosome <-  foreach(i = 1:n.chromosomes, .packages=c("seqinr", "VariantAnnotation", "Biostrings")) %dopar% { 
        
        message(paste("processing chromosome", i))
        
        postTotal.significant.i <- subset(postTotal.significant, postTotal.significant$contig == i)
        
        postTotal.significant.i["promoter_5kb"] <- "no"
        postTotal.significant.i["promoter_1kb"] <- "no"
        postTotal.significant.i["gene"] <- "no"
        postTotal.significant.i["five_prime_UTR"] <- "no"
        postTotal.significant.i["exon"] <- "no"
        postTotal.significant.i["intron"] <- "no"
        postTotal.significant.i["three_prime_UTR"] <- "no"
        postTotal.significant.i["post_gene_1kb"] <- "no"
        
        postTotal.significant.i["non_genic"] <- "no"
        
        postTotal.significant.i["gene.ID"] <- NA
        postTotal.significant.i["transcript.ID"] <- NA
        postTotal.significant.i["strand"] <- NA
        
        postTotal.significant.i["distance_to_gene"] <- NA
        
      
        df.gff.i <- subset(df.gene_annotation, df.gene_annotation$chr == i)
        
        l.df.gff.i <- vector(mode = "list", length = length(v.genePartitions))
        l.df.gff.i[[1]] <- subset(df.gff.i, df.gff.i$partition %in% v.genePartitions[1])
        l.df.gff.i[[2]] <- subset(df.gff.i, df.gff.i$partition %in% v.genePartitions[2])
        l.df.gff.i[[3]] <- subset(df.gff.i, df.gff.i$partition %in% v.genePartitions[3])
        l.df.gff.i[[4]] <- subset(df.gff.i, df.gff.i$partition %in% v.genePartitions[4])
        l.df.gff.i[[5]] <- subset(df.gff.i, df.gff.i$partition %in% v.genePartitions[5])
        
        # nearest gene distance 
        
        # preselect all in gene 
        pb <- txtProgressBar(min = 0, max = nrow(postTotal.significant.i), style = 3)
        for(j in 1:nrow(postTotal.significant.i)){
          
          setTxtProgressBar(pb, j)
      
          # non strand sensitive part #
          
          # in gene
          i.set <- which(l.df.gff.i[[1]]$pos.start < postTotal.significant.i$position[j] & postTotal.significant.i$position[j] < l.df.gff.i[[1]]$pos.stop)
          
          if(length(i.set) > 0){
            
            postTotal.significant.i$gene[j] <- "yes"
            postTotal.significant.i$gene.ID[j] <- l.df.gff.i[[1]]$gene.ID[i.set[1]]
            postTotal.significant.i$strand[j] <- l.df.gff.i[[1]]$strand[i.set[1]]
            
            
            postTotal.significant.i$distance_to_gene[j] <- 0
            
            # in 5 prime utr
            i.set <- which(l.df.gff.i[[2]]$pos.start < postTotal.significant.i$position[j] & postTotal.significant.i$position[j] < l.df.gff.i[[2]]$pos.stop)
            if(length(i.set) > 0){
              postTotal.significant.i$five_prime_UTR[j] <- "yes"
            }else{
              # in 3 prime utr
              i.set <- which(l.df.gff.i[[4]]$pos.start < postTotal.significant.i$position[j] & postTotal.significant.i$position[j] < l.df.gff.i[[4]]$pos.stop)
              if(length(i.set) > 0){
                postTotal.significant.i$three_prime_UTR[j] <- "yes"
              }else{
                # in CDS (exon)
                i.set <- which(l.df.gff.i[[5]]$pos.start < postTotal.significant.i$position[j] & postTotal.significant.i$position[j] < l.df.gff.i[[5]]$pos.stop)
                if(length(i.set) > 0){
                  postTotal.significant.i$exon[j] <- "yes"
                }else{
                  # i.set <- which(l.df.gff.i[[3]]$pos.start < postTotal.significant.i$position[j] & postTotal.significant.i$position[j] < l.df.gff.i[[3]]$pos.stop)
                  #if(length(i.set) > 0){  
                  postTotal.significant.i$intron[j] <- "yes"
                  #}
                }
              } 
            }
          }
        }
        close(pb)
        
        # strand sensitive part #
        
        # identify bQTL in 5kb and 1kb distance 
        l.df.gff.i[[1]]["up_5kb"] <- 0
        l.df.gff.i[[1]]["up_1kb"] <- 0
        l.df.gff.i[[1]]["start_gene"] <- 0
        l.df.gff.i[[1]]["end_gene"] <- 0
        l.df.gff.i[[1]]["down_1kb"] <- 0
        l.df.gff.i[[1]]["down_5kb"] <- 0
        
        for(k in 1:nrow(l.df.gff.i[[1]])){
          
          # strand sensitive 
          if(l.df.gff.i[[1]]$strand[k] == "+"){# promoter before (-) the gene
            
            l.df.gff.i[[1]]$up_5kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k] - 5000)
            l.df.gff.i[[1]]$up_1kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k] - 1000)
            l.df.gff.i[[1]]$start_gene[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k])
            
            l.df.gff.i[[1]]$end_gene[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k])
            l.df.gff.i[[1]]$down_1kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k] + 1000)
            l.df.gff.i[[1]]$down_5kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k] + 5000)
            
          }else if(l.df.gff.i[[1]]$strand[k] == "-"){# promoter after (+) the gene - perspective from positive strand 
            
            l.df.gff.i[[1]]$up_1kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k] + 1000)
            l.df.gff.i[[1]]$up_5kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k] + 5000)
            
            l.df.gff.i[[1]]$start_gene[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k])
            
            l.df.gff.i[[1]]$end_gene[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k])
            l.df.gff.i[[1]]$down_1kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k] - 1000)
            l.df.gff.i[[1]]$down_5kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k] - 5000)
            
            
          }
        }
        
        
        # search only for genes?
        i.set_upDown <- which(postTotal.significant.i$gene == "no")
        
        pb <- txtProgressBar(min = 0, max = length(i.set_upDown), style = 3)
        
        for(j in 1:length(i.set_upDown)){
          
          setTxtProgressBar(pb, j)
          
          # in 1 KB upstream - both strands 
          i.set <- which(l.df.gff.i[[1]]$up_1kb < postTotal.significant.i$position[i.set_upDown[j]] & postTotal.significant.i$position[i.set_upDown[j]] < l.df.gff.i[[1]]$start_gene 
                         | l.df.gff.i[[1]]$start_gene < postTotal.significant.i$position[i.set_upDown[j]] & postTotal.significant.i$position[i.set_upDown[j]] < l.df.gff.i[[1]]$up_1kb)
          
          if(length(i.set) > 0){
            
            if(postTotal.significant.i$gene[i.set_upDown[j]] == "no"){
              postTotal.significant.i$promoter_1kb[i.set_upDown[j]] <- "yes"
            }
            
            postTotal.significant.i$gene.ID[i.set_upDown[j]] <- l.df.gff.i[[1]]$gene.ID[i.set[1]]
            postTotal.significant.i$strand[i.set_upDown[j]] <- l.df.gff.i[[1]]$strand[i.set[1]]
           
            if(postTotal.significant.i$strand[i.set_upDown[j]] == "+"){
              postTotal.significant.i$distance_to_gene[i.set_upDown[j]] <- l.df.gff.i[[1]]$pos.start[i.set[1]] - postTotal.significant.i$position[i.set_upDown[j]]
            }else{
              postTotal.significant.i$distance_to_gene[i.set_upDown[j]] <- postTotal.significant.i$position[i.set_upDown[j]] - l.df.gff.i[[1]]$pos.stop[i.set[1]]
            }
            
            
            # take the closest
            # postTotal.significant.i$distance_to_gene[i.set_upDown[j]] <- 
            
          }else{
            
            # in 5 KB upstream - consider both strands
            i.set <- which(l.df.gff.i[[1]]$up_5kb < postTotal.significant.i$position[i.set_upDown[j]] & postTotal.significant.i$position[i.set_upDown[j]] < l.df.gff.i[[1]]$start_gene
                           | l.df.gff.i[[1]]$start_gene < postTotal.significant.i$position[i.set_upDown[j]] & postTotal.significant.i$position[i.set_upDown[j]] < l.df.gff.i[[1]]$up_5kb)
            
            if(length(i.set)  > 0){
              
              if(postTotal.significant.i$gene[i.set_upDown[j]] == "no"){
                postTotal.significant.i$promoter_5kb[i.set_upDown[j]] <- "yes"  
              }
              
              postTotal.significant.i$gene.ID[i.set_upDown[j]] <- l.df.gff.i[[1]]$gene.ID[i.set[1]]
              postTotal.significant.i$strand[i.set_upDown[j]] <- l.df.gff.i[[1]]$strand[i.set[1]]
              
              if(postTotal.significant.i$strand[i.set_upDown[j]] == "+"){
                postTotal.significant.i$distance_to_gene[i.set_upDown[j]] <- l.df.gff.i[[1]]$pos.start[i.set[1]] - postTotal.significant.i$position[i.set_upDown[j]]
              }else{
                postTotal.significant.i$distance_to_gene[i.set_upDown[j]] <- postTotal.significant.i$position[i.set_upDown[j]] - l.df.gff.i[[1]]$pos.stop[i.set[1]]
              }
              
            }
            
          }
          
          
          # now also duoble partition assignments 
          i.set <- which(l.df.gff.i[[1]]$end_gene < postTotal.significant.i$position[i.set_upDown[j]] & postTotal.significant.i$position[i.set_upDown[j]] < l.df.gff.i[[1]]$down_1kb
                         | l.df.gff.i[[1]]$down_1kb < postTotal.significant.i$position[i.set_upDown[j]] & postTotal.significant.i$position[i.set_upDown[j]] < l.df.gff.i[[1]]$end_gene)
          
          if(length(i.set)  > 0){
            
            if(postTotal.significant.i$gene[i.set_upDown[j]] == "no"){
              postTotal.significant.i$post_gene_1kb[i.set_upDown[j]] <- "yes"
            }
            
            postTotal.significant.i$gene.ID[i.set_upDown[j]] <- l.df.gff.i[[1]]$gene.ID[i.set[1]]
            postTotal.significant.i$strand[i.set_upDown[j]] <- l.df.gff.i[[1]]$strand[i.set[1]]
            
            if(postTotal.significant.i$strand[i.set_upDown[j]] == "-"){ # inverse because post gene
              postTotal.significant.i$distance_to_gene[i.set_upDown[j]] <- l.df.gff.i[[1]]$pos.start[i.set[1]] - postTotal.significant.i$position[i.set_upDown[j]]
            }else{
              postTotal.significant.i$distance_to_gene[i.set_upDown[j]] <- postTotal.significant.i$position[i.set_upDown[j]] - l.df.gff.i[[1]]$pos.stop[i.set[1]]
            }
          }
          # }else{
          #   
          #   i.set <- which(l.df.gff.i[[1]]$end_gene < postTotal.significant.i$position[i.set_upDown[j]] & postTotal.significant.i$position[i.set_upDown[j]] < l.df.gff.i[[1]]$down_5kb
          #                  | l.df.gff.i[[1]]$down_5kb < postTotal.significant.i$position[i.set_upDown[j]] & postTotal.significant.i$position[i.set_upDown[j]] < l.df.gff.i[[1]]$end_gene)
          #   
          #   if(length(i.set) > 0){
          #     
          #     if(postTotal.significant.i$gene[i.set_upDown[j]] == "no"){
          #       postTotal.significant.i$post_gene_5kb[i.set_upDown[j]] <- "yes"
          #     }
          #     
          #     postTotal.significant.i$gene.ID[i.set_upDown[j]] <- l.df.gff.i[[1]]$gene.ID[i.set[1]]
          #     
          #   }
          # }
        }
        close(pb)
        
        postTotal.significant.i
        
      }
      
      stopCluster(cl)
      print(Sys.time()-strt)
      
      
      df.bQTL_gene_partitioning <- c()
      for(i in 1:n.chromosomes){
        df.bQTL_gene_partitioning <- rbind(df.bQTL_gene_partitioning, l.bQTL_gene_partitioning_per_chromosome[[i]])
      }
      
  
      # if not found anywehere near or in gene
      i.set <- apply(df.bQTL_gene_partitioning, 1, function(m) {all(m[v.partitions] == "no")})
      i.set <- which(i.set == TRUE)
      df.bQTL_gene_partitioning$non_genic[i.set] <- "yes"
      
      # find the nearest gene for those too
      # df.bQTL_gene_partitioning <- l.bQTL_gene_partitioning[[s]]
  
      df.gene_annotation.set <- subset(df.gene_annotation, df.gene_annotation$partition == "gene")
      df.gene_annotation.plus <- subset(df.gene_annotation.set, df.gene_annotation.set$strand == "+")
      df.gene_annotation.minus <- subset(df.gene_annotation.set, df.gene_annotation.set$strand == "-")
      
       
      for(i in 1:n.chromosomes){
        
        df.gene_annotation.plus.i <- subset(df.gene_annotation.plus, df.gene_annotation.plus$chr  == i)
        df.gene_annotation.minus.i <- subset(df.gene_annotation.minus, df.gene_annotation.minus$chr  == i)
        
        for(j in 1:length(i.set)){
          
          chr <- df.bQTL_gene_partitioning$contig[i.set[j]]
          pos <- df.bQTL_gene_partitioning$position[i.set[j]]
    
          if(chr == i){
          
            v.genes <- c(df.gene_annotation.plus.i$gene.ID,  df.gene_annotation.minus.i$gene.ID)
            v.distances <- c((df.gene_annotation.plus.i$pos.start - pos), (pos - df.gene_annotation.minus.i$pos.stop))
            v.strands <- c(rep("+", nrow(df.gene_annotation.plus.i)), rep("-", nrow(df.gene_annotation.minus.i)))
            
            idx <- which(v.distances > 0)
            v.genes <- v.genes[idx]
            v.distances <- v.distances[idx]
            v.strands <- v.strands[idx]
            
            idx <- which(v.distances == min(v.distances))
           
            df.bQTL_gene_partitioning$distance_to_gene[i.set[j]] <-  v.distances[idx[1]]
            df.bQTL_gene_partitioning$gene.ID[i.set[j]] <-  v.genes[idx[1]]
            df.bQTL_gene_partitioning$strand[i.set[j]] <-  v.strands[idx[1]]
            
          }
          
        }
        
      }
      
    
      for(j in 1:length(v.partitions)){
        print(v.partitions[j])
        tmp <- subset(df.bQTL_gene_partitioning, df.bQTL_gene_partitioning[,v.partitions[j]] == "yes")
        print(mean(abs(tmp$POSTfreq - 0.5)))
      }
      
      # numbers - summary
      df.partitions[1,s+1] <- table(df.bQTL_gene_partitioning$promoter_5kb)[2]
      df.partitions[2,s+1] <- table(df.bQTL_gene_partitioning$promoter_1kb)[2]
      df.partitions[3,s+1] <- table(df.bQTL_gene_partitioning$gene)[2]
      df.partitions[4,s+1] <- table(df.bQTL_gene_partitioning$five_prime_UTR)[2]
      df.partitions[5,s+1] <- table(df.bQTL_gene_partitioning$exon)[2]
      df.partitions[6,s+1] <- table(df.bQTL_gene_partitioning$intron)[2]
      df.partitions[7,s+1] <- table(df.bQTL_gene_partitioning$three_prime_UTR)[2]
      df.partitions[8,s+1] <- table(df.bQTL_gene_partitioning$post_gene_1kb)[2]
      # df.partitions[9,s+1] <- table(df.bQTL_gene_partitioning$post_gene_5kb)[2]
      
      # adapted - no post 5 kb
      df.partitions[9,s+1] <- nrow(postTotal.significant) - sum(df.partitions[1:8,s+1])
      df.partitions[10,s+1] <- nrow(postTotal.significant)
      
      
      # saveRDS(df.partitions, "df.partitions.rds")
      # test <- apply(df.bQTL_gene_partitioning, 1, function(m) {length(which(m[c(17:18,20:24)] == "yes"))})
      
      # i.set <- apply( l.bQTL_gene_partitioning[[2]],1, function(m) {all(m[4:12] == "no")})
      
      l.bQTL_gene_partitioning[[s]] <- df.bQTL_gene_partitioning
      
    } # 27.04. closed
    
    saveRDS(l.bQTL_gene_partitioning, paste("tmp/l.bQTL_gene_partitioning_withGeneDistances_", timeStamp, ".rds", sep = ""))
    
  }else{
    l.bQTL_gene_partitioning <- readRDS(paste("tmp/l.bQTL_gene_partitioning_", timeStamp, ".rds", sep = ""))  
    for(s in 1:length(l.bQTL_gene_partitioning)){
      
      df.bQTL_gene_partitioning <- l.bQTL_gene_partitioning[[s]]
      
      df.partitions[1,s+1] <- table(df.bQTL_gene_partitioning$promoter_5kb)[2]
      df.partitions[2,s+1] <- table(df.bQTL_gene_partitioning$promoter_1kb)[2]
      df.partitions[3,s+1] <- table(df.bQTL_gene_partitioning$gene)[2]
      df.partitions[4,s+1] <- table(df.bQTL_gene_partitioning$five_prime_UTR)[2]
      df.partitions[5,s+1] <- table(df.bQTL_gene_partitioning$exon)[2]
      df.partitions[6,s+1] <- table(df.bQTL_gene_partitioning$intron)[2]
      df.partitions[7,s+1] <- table(df.bQTL_gene_partitioning$three_prime_UTR)[2]
      df.partitions[8,s+1] <- table(df.bQTL_gene_partitioning$post_gene_1kb)[2]
      
      df.partitions[9,s+1] <- nrow(l.bQTL_gene_partitioning[[s]]) - sum(df.partitions[c(1,2,4,5,6,7,8),s+1])
      
      df.partitions[10,s+1] <- nrow(l.bQTL_gene_partitioning[[s]])
      
    }
  }
  
  
  message("Figure 2 - Relative comparative partition analysis ")
  
  df.partitions["percentage_significant"] <- df.partitions[,2] / df.partitions[10,2]
  df.partitions["percentage_non_significant"] <- df.partitions[,3] / df.partitions[10,3]
  
  write.csv(df.partitions, paste("output/df.partition_plus_", timeStamp, ".csv", sep = ""))
  
  # length(unique(l.bQTL_gene_partitioning[[1]]$gene.ID))
  # post-frequency and expression level levels
  
  
  
  write.csv(df.partitions, paste("output/df.partition_plus_", timeStamp, ".csv", sep = ""))
  
  # duplicate removal
  
  # if(FALSE){
  #   df.bQTL_gene_partitioning_bg_equivalent <- read.csv("output/df.bQTL_gene_partitioning_bg_equivalent.csv")
  #   l.bQTL_gene_partitioning[[2]] <- df.bQTL_gene_partitioning_bg_equivalent
  #   
  #   
  #   v.partitions.duplicates <- v.partitions[c(1,2,8)]
  #   
  #   for(i in 1:2){
  #     i.set <- which(apply(l.bQTL_gene_partitioning[[i]][,v.partitions.duplicates], 1, function(m) length(which(m == "yes"))) == 2)
  #     i.set <- (!1:nrow(l.bQTL_gene_partitioning[[i]]) %in% i.set)
  #     l.bQTL_gene_partitioning[[i]] <- l.bQTL_gene_partitioning[[i]][i.set, ]
  #   }
  #   
  #   
  #   saveRDS(l.bQTL_gene_partitioning, paste("tmp/l.bQTL_gene_partitioning_noDublicates",timeStamp, ".rds", sep = ""))  
  #   
  #   
  # }  
  

}















# if(FALSE){
#   
#   
#   
#   df.partitions["relative_distribution_significant"] <- df.partitions$percentage_significant / c(4,1,3,0.411,0.865, 2.1,0.417,1,4,20.24, 1000)
#   df.partitions["relative_distribution_nonsignificant"] <- df.partitions$percentage_non_significant / c(4,1,3,0.411,0.865, 2.1,0.417,1,4,20.24,1000)
#   
#   # figure 2 - bargraphs
#   
#   
#   
#   # with duplicates 
#   
# }


if(FALSE){
  
  message("perform peak partitioning")
  
  if(!b.loadGenePartitioning){
  
    strt<-Sys.time() 
    cl<-makeCluster(min(n.chromosomes, n.cpus))
    registerDoParallel(cl)
    
    l.peak_gene_partitioning_per_chromosome <-  foreach(i = 1:n.chromosomes, .packages=c("seqinr", "VariantAnnotation", "Biostrings")) %dopar% { 
      #for(i in 1:n.chromosomes){
      
      message(paste("processing chromosome", i))
      
      df.peaks.i <- subset(df.peaks, df.peaks$seqnames == i)
      
      df.peaks.i["promoter_5kb"] <- "no"
      df.peaks.i["promoter_1kb"] <- "no"
      df.peaks.i["gene"] <- "no"
      df.peaks.i["five_prime_UTR"] <- "no"
      df.peaks.i["exon"] <- "no"
      df.peaks.i["intron"] <- "no"
      df.peaks.i["three_prime_UTR"] <- "no"
      df.peaks.i["post_gene_1kb"] <- "no"
      df.peaks.i["non_genic"] <- "no"
      
      df.peaks.i["gene.ID"] <- NA
      
      
      df.gff.i <- subset(df.gene_annotation, df.gene_annotation$chr == i)
      
      l.df.gff.i <- vector(mode = "list", length = length(v.genePartitions))
      l.df.gff.i[[1]] <- subset(df.gff.i, df.gff.i$partition %in% v.genePartitions[1])
      l.df.gff.i[[2]] <- subset(df.gff.i, df.gff.i$partition %in% v.genePartitions[2])
      l.df.gff.i[[3]] <- subset(df.gff.i, df.gff.i$partition %in% v.genePartitions[3])
      l.df.gff.i[[4]] <- subset(df.gff.i, df.gff.i$partition %in% v.genePartitions[4])
      l.df.gff.i[[5]] <- subset(df.gff.i, df.gff.i$partition %in% v.genePartitions[5])
      
      # preselect all in gene 
      pb <- txtProgressBar(min = 0, max = nrow(df.peaks.i), style = 3)
      
      for(j in 1:nrow(df.peaks.i)){
        
        setTxtProgressBar(pb, j)
        
        # in gene 
        i.set <- which(l.df.gff.i[[1]]$pos.start < df.peaks.i$posPeak[j] & df.peaks.i$posPeak[j] < l.df.gff.i[[1]]$pos.stop)
        
        if(length(i.set) > 0){
          
          df.peaks.i$gene[j] <- "yes"
          df.peaks.i$gene.ID[j] <- l.df.gff.i[[1]]$gene.ID[i.set[1]]
          
          # in 5 prime utr
          i.set <- which(l.df.gff.i[[2]]$pos.start < df.peaks.i$posPeak[j] & df.peaks.i$posPeak[j] < l.df.gff.i[[2]]$pos.stop)
          if(length(i.set) > 0){
            df.peaks.i$five_prime_UTR[j] <- "yes"
          }else{
            # in 3 prime utr
            i.set <- which(l.df.gff.i[[4]]$pos.start < df.peaks.i$posPeak[j] & df.peaks.i$posPeak[j] < l.df.gff.i[[4]]$pos.stop)
            if(length(i.set) > 0){
              df.peaks.i$three_prime_UTR[j] <- "yes"
            }else{
              # in CDS (exon)
              i.set <- which(l.df.gff.i[[5]]$pos.start < df.peaks.i$posPeak[j] & df.peaks.i$posPeak[j] < l.df.gff.i[[5]]$pos.stop)
              if(length(i.set) > 0){
                df.peaks.i$exon[j] <- "yes"
              }else{
                # i.set <- which(l.df.gff.i[[3]]$pos.start < postTotal.significant.i$position[j] & postTotal.significant.i$position[j] < l.df.gff.i[[3]]$pos.stop)
                #if(length(i.set) > 0){  
                df.peaks.i$intron[j] <- "yes"
                #}
              }
            } 
          }
        }
      }
      close(pb)
      
      # identify bQTL in 5kb and 1kb distance 
      l.df.gff.i[[1]]["up_5kb"] <- 0
      l.df.gff.i[[1]]["up_1kb"] <- 0
      l.df.gff.i[[1]]["start_gene"] <- 0
      l.df.gff.i[[1]]["end_gene"] <- 0
      l.df.gff.i[[1]]["down_1kb"] <- 0
      l.df.gff.i[[1]]["down_5kb"] <- 0
      
      for(k in 1:nrow(l.df.gff.i[[1]])){
        
        if(l.df.gff.i[[1]]$strand[k] == "+"){
          
          l.df.gff.i[[1]]$up_5kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k] - 5000)
          l.df.gff.i[[1]]$up_1kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k] - 1000)
          l.df.gff.i[[1]]$start_gene[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k])
          
          l.df.gff.i[[1]]$end_gene[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k])
          l.df.gff.i[[1]]$down_1kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k] + 1000)
          l.df.gff.i[[1]]$down_5kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k] + 5000)
          
        }else if(l.df.gff.i[[1]]$strand[k] == "-"){
          
          l.df.gff.i[[1]]$up_1kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k] + 1000)
          l.df.gff.i[[1]]$up_5kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k] + 5000)
          
          l.df.gff.i[[1]]$start_gene[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k])
          
          l.df.gff.i[[1]]$end_gene[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k])
          l.df.gff.i[[1]]$down_1kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k] - 1000)
          l.df.gff.i[[1]]$down_5kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k] - 5000)
          
        }
      }
      
      
      # search only for genes?
      i.set_upDown <- which(df.peaks.i$gene == "no")
      
      pb <- txtProgressBar(min = 0, max = length(i.set_upDown), style = 3)
      
      for(j in 1:length(i.set_upDown)){
        
        setTxtProgressBar(pb, j)
        
        # in 1 KB upstream 
        i.set <- which(l.df.gff.i[[1]]$up_1kb < df.peaks.i$posPeak[i.set_upDown[j]] & df.peaks.i$posPeak[i.set_upDown[j]] < l.df.gff.i[[1]]$start_gene 
                       | l.df.gff.i[[1]]$start_gene < df.peaks.i$posPeak[i.set_upDown[j]] & df.peaks.i$posPeak[i.set_upDown[j]] < l.df.gff.i[[1]]$up_1kb)
        
        if(length(i.set)  > 0){
          
          
          if(df.peaks.i$gene[i.set_upDown[j]] == "no"){
            df.peaks.i$promoter_1kb[i.set_upDown[j]] <- "yes"
          }
          
          df.peaks.i$gene.ID[i.set_upDown[j]] <- l.df.gff.i[[1]]$gene.ID[i.set[1]]
          
        }else{
          
          # in 5 KB upstream 
          i.set <- which(l.df.gff.i[[1]]$up_5kb < df.peaks.i$posPeak[i.set_upDown[j]] & df.peaks.i$posPeak[i.set_upDown[j]] < l.df.gff.i[[1]]$start_gene
                         | l.df.gff.i[[1]]$start_gene < df.peaks.i$posPeak[i.set_upDown[j]] & df.peaks.i$posPeak[i.set_upDown[j]] < l.df.gff.i[[1]]$up_5kb)
          
          if(length(i.set)  > 0){
            
            if(df.peaks.i$gene[i.set_upDown[j]] == "no"){
              df.peaks.i$promoter_5kb[i.set_upDown[j]] <- "yes"  
            }
            
            df.peaks.i$gene.ID[i.set_upDown[j]] <- l.df.gff.i[[1]]$gene.ID[i.set[1]]
            
          }
          
        }
        
        
        # now also duoble partition assignments 
        i.set <- which(l.df.gff.i[[1]]$end_gene < df.peaks.i$posPeak[i.set_upDown[j]] & df.peaks.i$posPeak[i.set_upDown[j]] < l.df.gff.i[[1]]$down_1kb
                       | l.df.gff.i[[1]]$down_1kb < df.peaks.i$posPeak[i.set_upDown[j]] & df.peaks.i$posPeak[i.set_upDown[j]] < l.df.gff.i[[1]]$end_gene)
        
        if(length(i.set)  > 0){
          
          if(df.peaks.i$gene[i.set_upDown[j]] == "no"){
            df.peaks.i$post_gene_1kb[i.set_upDown[j]] <- "yes"
          }
          
          df.peaks.i$gene.ID[i.set_upDown[j]] <- l.df.gff.i[[1]]$gene.ID[i.set[1]]
          
        }
        # }else{
        #   
        #   i.set <- which(l.df.gff.i[[1]]$end_gene < postTotal.significant.i$position[i.set_upDown[j]] & postTotal.significant.i$position[i.set_upDown[j]] < l.df.gff.i[[1]]$down_5kb
        #                  | l.df.gff.i[[1]]$down_5kb < postTotal.significant.i$position[i.set_upDown[j]] & postTotal.significant.i$position[i.set_upDown[j]] < l.df.gff.i[[1]]$end_gene)
        #   
        #   if(length(i.set) > 0){
        #     
        #     if(postTotal.significant.i$gene[i.set_upDown[j]] == "no"){
        #       postTotal.significant.i$post_gene_5kb[i.set_upDown[j]] <- "yes"
        #     }
        #     
        #     postTotal.significant.i$gene.ID[i.set_upDown[j]] <- l.df.gff.i[[1]]$gene.ID[i.set[1]]
        #     
        #   }
        # }
      }
      close(pb)
      
      df.peaks.i
    }
    
    stopCluster(cl)
    print(Sys.time()-strt)
    
    
    df.peak_gene_partitioning <- c()
    for(i in 1:n.chromosomes){
      df.peak_gene_partitioning <- rbind(df.peak_gene_partitioning, l.peak_gene_partitioning_per_chromosome[[i]])
    }
    
    # if not found anywehere near or in gene
    i.set <- apply(df.peak_gene_partitioning, 1, function(m) {all(m[v.partitions] == "no")})
    df.peak_gene_partitioning$non_genic[i.set] <- "yes"
    
    
    
    saveRDS(df.peak_gene_partitioning, paste("tmp/df.peak_gene_partitioning_", timeStamp, ".rds", sep = ""))
    
    
    # write.csv(df.bQTL_gene_partitioning, paste("tmp/df.bQTL_gene_partitioning_", timeStamp, ".csv", sep = ""))
    # write.csv(df.partitions, paste("tmp/df.partitions_", timeStamp, ".csv", sep = ""))
    
  }else{
    
    df.peak_gene_partitioning <- readRDS(paste("tmp/df.peak_gene_partitioning_", timeStamp, ".rds", sep = ""))  
  
    # non duplicate removal
    v.distribution <- apply(df.peak_gene_partitioning[,v.partitions], 2, table)["yes",]  
   
  }
  

}



df.ChipSeq.filtered = read.table("D:/junkDNA.ai/Projects/HASCHSEQ/HaschSeqBasePrototype/output/gene_partitioning/S2.txt", header = T)




postTotal.significant <- df.ChipSeq.filtered
postTotal.significant <- l.selection[[1]]
postTotal.significant <- df.snp_positions
postTotal.significant <- df.peaks

df.gene_conversion.AGPv3_to_AGPv4 <- df.geneID_conversion
df.gene_function <- df.gene_function

# ASBs
pos_peak = "position" 
chr_Nr = "contig"

# Peaks
pos_peak = "posPeak" 
chr_Nr = "seqnames"

# df.ChipSeq.filtered
pos_peak = "pos" 
chr_Nr = "chr"



df.ChipSeq_gene_partitioning = add_gene_partition(df.ChipSeq.filtered, 
                                                 pos_peak = "pos",
                                                 chr_Nr = "chr", 
                                                 v.genePartitions = c("gene", "five_prime_UTR",  "CDS", "three_prime_UTR", "exon"),
                                                 df.gene_annotation,
                                                 n.cpus=1)

 
# v.genePartitions = c("promoter_5kb", "promoter_1kb","gene", "five_prime_UTR", "intron", "CDS", "three_prime_UTR", "exon", "post_gene_1kb")
table(df.ChipSeq_gene_partitioning$non_genic)
test = subset(df.ChipSeq_gene_partitioning, df.ChipSeq_gene_partitioning$non_genic == "no")
length(unique(test$gene.ID))

test2 = add_features_to_gene_partitioning(test)
length(unique(test2$Arabidopsis_ortholog))

test2 = test2[,c("gene.ID", "Arabidopsis_ortholog")]
test3 = subset(test2, !is.na(test2$Arabidopsis_ortholog))

write.table(test3, "S4.txt", sep ="\t", row.names = F, quote = F)

writeRDS()
# Table S2: List of significant ZmBZR1 binding sites in meristem enriched tissue of B73 using a p-value cut-off of xxx
write.table(df.ChipSeq_gene_partitioning, paste(folder_output, "/gene_partitioning/S2.txt", sep = ""),  row.names = FALSE, quote = FALSE, sep ="\t")

####


df.peaks_gene_partitioning = add_gene_partition(df.peaks, 
                                                  pos_peak = "posPeak",
                                                  chr_Nr = "seqnames", 
                                                  v.genePartitions = c("gene", "five_prime_UTR",  "CDS", "three_prime_UTR", "exon"),
                                                  df.gene_annotation)



write.table(df.peaks_gene_partitioning, paste(folder_output, "/S4.txt", sep = ""),  row.names = FALSE, quote = FALSE, sep ="\t")



df.ASBs_gene_partitioning = add_gene_partition(df.ASBs, 
                                                pos_peak = "position",
                                                chr_Nr = "contig", 
                                                v.genePartitions = c("gene", "five_prime_UTR",  "CDS", "three_prime_UTR", "exon"),
                                                df.gene_annotation)






### double check with 


length(intersect(df.ChipSeq_gene_partitioning$gene.ID, as.character(df.bQTL_RNAseq$gene.ID[df.bQTL_RNAseq$mode == "up"])))
length(intersect(df.ChipSeq_gene_partitioning$gene.ID, as.character(df.bQTL_RNAseq$gene.ID[df.bQTL_RNAseq$mode == "down"])))




# ath chip seq - dark
df.At_BZR1_ChipSeqTargets <- read.table("data/ArabidopsisScriptsAndDatasets/At_BZR1_ChipSeqTargets.txt", header = TRUE, sep ="\t", quote = "", stringsAsFactors = FALSE)
names(df.At_BZR1_ChipSeqTargets) <- "locus"

# ath chIP chip - light
df.At_BZR1_ChIPchipTargets <- read.table("data/ArabidopsisScriptsAndDatasets/At_BZR1_ChIPchip.txt", header = TRUE, sep ="\t", quote = "", stringsAsFactors = FALSE)

length(intersect(df.bQTL_gene_partitioning$Arabidopsis_ortholog, df.At_BZR1_ChipSeqTargets$locus))
length(intersect(df.bQTL_gene_partitioning$Arabidopsis_ortholog, df.At_BZR1_ChIPchipTargets$locus))

df.bQTL_RNAseq


df.bQTL_gene_partitioning = df.ChipSeq_gene_partitioning















df.bQTL_gene_partitioning = addNearestGene(postTotal.significant, 
                                           pos_peak, chr_Nr,  
                                           df.gene_annotation = NULL, 
                                           df.gene_annotation.AGPv3 = NULL, 
                                           df.gene_function,
                                           n.cpus = 4)






# df.snp_positions - needs to have identical syntax
addNearestGene <- function(postTotal.significant, 
                           pos_peak, chr_Nr,  
                           df.gene_annotation = NULL, 
                           df.gene_annotation.AGPv3 = NULL, 
                           df.gene_function = NULL,
                           n.cpus = 4){
  
  idx.pos <- which( colnames(postTotal.significant)== pos_peak)
  idx.chr <- which( colnames(postTotal.significant) == chr_Nr)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  message("identifying gene partitions - positions of ASBs with ")
  
  df.bQTL_gene_partitioning <- c()
  
  strt<-Sys.time() 
  cl<-makeCluster(min(n.chromosomes, n.cpus))
  registerDoParallel(cl)
  
  l.bQTL_gene_partitioning_per_chromosome <-  foreach(i = 1:n.chromosomes, .packages=c("seqinr", "VariantAnnotation", "Biostrings")) %dopar% { 
    #for(i in 1:n.chromosomes){
    
    message(paste("processing chromosome", i))
    
    postTotal.significant.i <- subset(postTotal.significant, postTotal.significant[,idx.chr] == i)
    
    postTotal.significant.i["promoter_5kb"] <- "no"
    postTotal.significant.i["promoter_1kb"] <- "no"
    postTotal.significant.i["gene"] <- "no"
    postTotal.significant.i["five_prime_UTR"] <- "no"
    postTotal.significant.i["exon"] <- "no"
    postTotal.significant.i["intron"] <- "no"
    postTotal.significant.i["three_prime_UTR"] <- "no"
    postTotal.significant.i["post_gene_1kb"] <- "no"
    
    postTotal.significant.i["gene.ID"] <- NA
    
    df.gff.i <- subset(df.gene_annotation, df.gene_annotation$chr == i)
    
    l.df.gff.i <- vector(mode = "list", length = length(v.genePartitions))
    l.df.gff.i[[1]] <- subset(df.gff.i, df.gff.i$partition %in% v.genePartitions[1])
    l.df.gff.i[[2]] <- subset(df.gff.i, df.gff.i$partition %in% v.genePartitions[2])
    l.df.gff.i[[3]] <- subset(df.gff.i, df.gff.i$partition %in% v.genePartitions[3])
    l.df.gff.i[[4]] <- subset(df.gff.i, df.gff.i$partition %in% v.genePartitions[4])
    l.df.gff.i[[5]] <- subset(df.gff.i, df.gff.i$partition %in% v.genePartitions[5])
    
    # preselect all in gene 
    pb <- txtProgressBar(min = 0, max = nrow(postTotal.significant.i), style = 3)
    
    for(j in 1:nrow(postTotal.significant.i)){
      
      setTxtProgressBar(pb, j)
      
      # in gene 
      i.set <- which(l.df.gff.i[[1]]$pos.start < postTotal.significant.i[j,idx.pos] & postTotal.significant.i[j,idx.pos] < l.df.gff.i[[1]]$pos.stop)
      
      if(length(i.set) > 0){
        
        postTotal.significant.i$gene[j] <- "yes"
        postTotal.significant.i$gene.ID[j] <- l.df.gff.i[[1]]$gene.ID[i.set[1]]
        
        # in 5 prime utr
        i.set <- which(l.df.gff.i[[2]]$pos.start < postTotal.significant.i[j,idx.pos] & postTotal.significant.i[j,idx.pos] < l.df.gff.i[[2]]$pos.stop)
        if(length(i.set) > 0){
          postTotal.significant.i$five_prime_UTR[j] <- "yes"
        }else{
          # in 3 prime utr
          i.set <- which(l.df.gff.i[[4]]$pos.start < postTotal.significant.i[j,idx.pos] & postTotal.significant.i[j,idx.pos] < l.df.gff.i[[4]]$pos.stop)
          if(length(i.set) > 0){
            postTotal.significant.i$three_prime_UTR[j] <- "yes"
          }else{
            # in CDS (exon)
            i.set <- which(l.df.gff.i[[5]]$pos.start < postTotal.significant.i[j,idx.pos] & postTotal.significant.i[j,idx.pos] < l.df.gff.i[[5]]$pos.stop)
            if(length(i.set) > 0){
              postTotal.significant.i$exon[j] <- "yes"
            }else{
              # i.set <- which(l.df.gff.i[[3]]$pos.start < postTotal.significant.i$position[j] & postTotal.significant.i$position[j] < l.df.gff.i[[3]]$pos.stop)
              #if(length(i.set) > 0){  
              postTotal.significant.i$intron[j] <- "yes"
              #}
            }
          } 
        }
      }
    }
    close(pb)
    
    # identify bQTL in 5kb and 1kb distance 
    l.df.gff.i[[1]]["up_5kb"] <- 0
    l.df.gff.i[[1]]["up_1kb"] <- 0
    l.df.gff.i[[1]]["start_gene"] <- 0
    l.df.gff.i[[1]]["end_gene"] <- 0
    l.df.gff.i[[1]]["down_1kb"] <- 0
    
    for(k in 1:nrow(l.df.gff.i[[1]])){
      if(l.df.gff.i[[1]]$strand[k] == "+"){
        l.df.gff.i[[1]]$up_5kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k] - 5000)
        l.df.gff.i[[1]]$up_1kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k] - 1000)
        l.df.gff.i[[1]]$start_gene[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k])
        l.df.gff.i[[1]]$end_gene[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k])
        l.df.gff.i[[1]]$down_1kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k] + 1000)
      }else if(l.df.gff.i[[1]]$strand[k] == "-"){
        l.df.gff.i[[1]]$up_1kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k] + 1000)
        l.df.gff.i[[1]]$up_5kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k] + 5000)
        l.df.gff.i[[1]]$start_gene[k] <- as.numeric(l.df.gff.i[[1]]$pos.stop[k])
        l.df.gff.i[[1]]$end_gene[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k])
        l.df.gff.i[[1]]$down_1kb[k] <- as.numeric(l.df.gff.i[[1]]$pos.start[k] - 1000)
      }
    }
    
    
    
    i.set_upDown <- which(postTotal.significant.i$gene == "no")
    
    pb <- txtProgressBar(min = 0, max = length(i.set_upDown), style = 3)
    
    for(j in 1:length(i.set_upDown)){
      
      setTxtProgressBar(pb, j)
      
      # in 1 KB upstream 
      i.set <- which(l.df.gff.i[[1]]$up_1kb < postTotal.significant.i[i.set_upDown[j],idx.pos] & postTotal.significant.i[i.set_upDown[j],idx.pos] < l.df.gff.i[[1]]$start_gene 
                     | l.df.gff.i[[1]]$start_gene < postTotal.significant.i[i.set_upDown[j],idx.pos] & postTotal.significant.i[i.set_upDown[j],idx.pos] < l.df.gff.i[[1]]$up_1kb)
      
      if(length(i.set)  > 0){
        postTotal.significant.i$promoter_1kb[i.set_upDown[j]] <- "yes"
        postTotal.significant.i$gene.ID[i.set_upDown[j]] <- l.df.gff.i[[1]]$gene.ID[i.set[1]]
      }else{
        # in 5 KB upstream 
        i.set <- which(l.df.gff.i[[1]]$up_5kb < postTotal.significant.i[i.set_upDown[j],idx.pos] & postTotal.significant.i[i.set_upDown[j],idx.pos] < l.df.gff.i[[1]]$start_gene
                       | l.df.gff.i[[1]]$start_gene < postTotal.significant.i[i.set_upDown[j],idx.pos] & postTotal.significant.i[i.set_upDown[j],idx.pos] < l.df.gff.i[[1]]$up_5kb)
        if(length(i.set)  > 0){
          postTotal.significant.i$promoter_5kb[i.set_upDown[j]] <- "yes"
          postTotal.significant.i$gene.ID[i.set_upDown[j]] <- l.df.gff.i[[1]]$gene.ID[i.set[1]]
        }
      }
      
      
      # now also duoble partition assignments 
      i.set <- which(l.df.gff.i[[1]]$end_gene < postTotal.significant.i[i.set_upDown[j],idx.pos] & postTotal.significant.i[i.set_upDown[j],idx.pos] < l.df.gff.i[[1]]$down_1kb
                     | l.df.gff.i[[1]]$down_1kb < postTotal.significant.i[i.set_upDown[j],idx.pos] & postTotal.significant.i[i.set_upDown[j],idx.pos] < l.df.gff.i[[1]]$end_gene)
      
      if(length(i.set)  > 0){
        
        postTotal.significant.i$post_gene_1kb[i.set_upDown[j]] <- "yes"
        postTotal.significant.i$gene.ID[i.set_upDown[j]] <- l.df.gff.i[[1]]$gene.ID[i.set[1]]
        
      }else{
        i.set <- which(l.df.gff.i[[1]]$end_gene < postTotal.significant.i[i.set_upDown[j],idx.pos] & postTotal.significant.i[i.set_upDown[j],idx.pos] < l.df.gff.i[[1]]$down_5kb
                       | l.df.gff.i[[1]]$down_5kb < postTotal.significant.i[i.set_upDown[j],idx.pos] & postTotal.significant.i[i.set_upDown[j],idx.pos] < l.df.gff.i[[1]]$end_gene)
        if(length(i.set) > 0){
          postTotal.significant.i$post_gene_5kb[i.set_upDown[j]] <- "yes"
          postTotal.significant.i$gene.ID[i.set_upDown[j]] <- l.df.gff.i[[1]]$gene.ID[i.set[1]]
        }
      }
    }
    close(pb)
    postTotal.significant.i
  }
  stopCluster(cl)
  print(Sys.time()-strt)
  
  
  df.bQTL_gene_partitioning <- c()
  for(i in 1:n.chromosomes){
    df.bQTL_gene_partitioning <- rbind(df.bQTL_gene_partitioning, l.bQTL_gene_partitioning_per_chromosome[[i]])
  }
  
  return(df.bQTL_gene_partitioning)
  
}


  # 
  # i.set <- apply(df.bQTL_gene_partitioning, 1, function(m) {all(m[v.partitions] == "no")})
  # df.bQTL_gene_partitioning["nonGenic"] <- "no"
  # df.bQTL_gene_partitioning$nonGenic[i.set] <- "yes"
  # for(j in 1:length(v.partitions)){
  #   print(v.partitions[j])
  #   tmp <- subset(df.bQTL_gene_partitioning, df.bQTL_gene_partitioning[,v.partitions[j]] == "yes")
  #   print(mean(abs(tmp$POSTfreq - 0.5)))
  # }
  
  # numbers - summary
  df.partitions[1,s+1] <- table(df.bQTL_gene_partitioning$promoter_5kb)[2]
  df.partitions[2,s+1] <- table(df.bQTL_gene_partitioning$promoter_1kb)[2]
  df.partitions[3,s+1] <- table(df.bQTL_gene_partitioning$gene)[2]
  df.partitions[4,s+1] <- table(df.bQTL_gene_partitioning$five_prime_UTR)[2]
  df.partitions[5,s+1] <- table(df.bQTL_gene_partitioning$exon)[2]
  df.partitions[6,s+1] <- table(df.bQTL_gene_partitioning$intron)[2]
  df.partitions[7,s+1] <- table(df.bQTL_gene_partitioning$three_prime_UTR)[2]
  df.partitions[8,s+1] <- table(df.bQTL_gene_partitioning$post_gene_1kb)[2]
  df.partitions[9,s+1] <- table(df.bQTL_gene_partitioning$post_gene_5kb)[2]
  
  df.partitions[10,s+1] <- nrow(postTotal.significant) - sum(df.partitions[1:9,s+1])
  df.partitions[11,s+1] <- nrow(postTotal.significant)
  
  # saveRDS(df.partitions, "df.partitions.rds")
  # test <- apply(df.bQTL_gene_partitioning, 1, function(m) {length(which(m[c(17:18,20:24)] == "yes"))})
  
  # i.set <- apply( l.bQTL_gene_partitioning[[2]],1, function(m) {all(m[4:12] == "no")})
  
  print(df.partitions)
  
  length(unique(df.bQTL_gene_partitioning$gene.ID))
  
  # add the AGPv3 gene annotation
  if(!is.null(df.gene_conversion.AGPv3_to_AGPv4)){
    df.bQTL_gene_partitioning["gene.ID.AGPv3"] <- NA
    #idx.geneID.AGPv3 <- which(df.bQTL_gene_partitioning$gene.ID %in% df.gene_conversion.AGPv3_to_AGPv4$gene.ID.AGPv4)
    for(i in 1:nrow(df.bQTL_gene_partitioning)){
      idx.AGPv3 <- which(df.gene_conversion.AGPv3_to_AGPv4$gene.ID.AGPv4 == df.bQTL_gene_partitioning$gene.ID[i])
      if(length(idx.AGPv3) > 0){
        df.bQTL_gene_partitioning$gene.ID.AGPv3[i] <-  df.gene_conversion.AGPv3_to_AGPv4$gene.ID.AGPv3[idx.AGPv3[1]]
      }
    }
  }
  
  
  # gene function annotation based on AGPv4
  if(!is.null(df.gene_function)){
    df.bQTL_gene_partitioning["gene.function"] <- NA
    #idx.geneID_to_function <- which(df.bQTL_gene_partitioning$gene.ID %in% df.gene_function$gene.ID)
    for(i in 1:nrow(df.bQTL_gene_partitioning)){
      idx.geneID_to_function <- which(df.gene_function$gene.ID == df.bQTL_gene_partitioning$gene.ID[i])
      if(length(idx.geneID_to_function) > 0){
        df.bQTL_gene_partitioning$gene.function[i] <-  df.gene_function$gene.function[idx.geneID_to_function[1]]
      }
    }
  }
  
  
  
  # add arabidopsis ortholog
  if(!is.null(df.gene_orthologs)){
    df.bQTL_gene_partitioning["Arabidopsis_ortholog"] <- NA
    for(i in 1:nrow(df.bQTL_gene_partitioning)){
      idx.AGPv3 <- which(df.gene_orthologs$locusName == df.bQTL_gene_partitioning$gene.ID.AGPv3[i])
      if(length(idx.AGPv3) > 0){
        df.bQTL_gene_partitioning$Arabidopsis_ortholog[i] <-  df.gene_orthologs$Best.hit.arabi.name[idx.AGPv3[1]]
      }
    }
  }
  
  
  write.csv(df.bQTL_gene_partitioning, "output/df.bQTL_gene_partitioning_peaks.csv", row.names = FALSE)
  
  
  # gene function enrichment tests 
  # BARPLOTS mit STERNCHEN (if p < 0.05) - identify biological processes
  
  
  df.snp_genes_with_functions <- unique(df.bQTL_gene_partitioning[,c("gene.ID", "gene.function")])
  df.snp_genes_with_functions <- subset(df.snp_genes_with_functions, !is.na(df.snp_genes_with_functions$gene.ID))
  
  tb.snp_genes_with_functions  <- table(df.genes_with_functions$gene.function)
  v.snp_genes_with_functions  <- unique(df.genes_with_functions$gene.function)
  v.snp_genes_with_functions <- v.snp_genes_with_functions[!is.na(v.snp_genes_with_functions)]
  
  tb.gene_function <- table(df.gene_function$gene.function)
  
  
  # df.GO.annot.set <- subset(df.GO.annot, df.GO.annot$acc. %in% df.bQTL_gene_partitioning_with_gwas$gene.ID)
  # 
  # df.tb.annot <- as.data.frame(table(df.GO.annot$Term[df.GO.annot$Term != ""]), stringsAsFactors = FALSE)
  # df.tb.annot.set <- as.data.frame(table(df.GO.annot.set$Term[df.GO.annot.set$Term != ""]), stringsAsFactors = FALSE)
  # df.tb.annot.set <- merge(df.tb.annot.set, df.tb.annot, by = "Var1")
  # 
  ## filter go - evidence codes, biological process
  df.gene_function_enrichment <- data.frame(gene_function = v.snp_genes_with_functions, 
                                            p.val = rep(1, length(v.snp_genes_with_functions)),
                                            n.genes = numeric(length(v.snp_genes_with_functions)), 
                                            foldchange = numeric(length(v.snp_genes_with_functions)), 
                                            stringsAsFactors = FALSE)
  
  for(j in 1:length(v.snp_genes_with_functions)){  
    
    ### global enrichment test - gene basis
    hitInSample <- tb.snp_genes_with_functions[v.snp_genes_with_functions[j]]
    sampleSize <- sum(tb.snp_genes_with_functions)
    
    hitInPop <- tb.gene_function[v.snp_genes_with_functions[j]] #sum(tb.rate_limiting_domains$Freq)
    popSize <- sum(tb.gene_function)
    
    failInPop <- popSize - hitInPop #(nrow(df.global.domains) - hitInPop)
    p.val <- print(phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE))
    
    #fisher.test(matrix(c(hitInSample-1, hitInPop, failInPop, sampleSize), 2, 2), alternative='less');
    foldChange <- (hitInSample / sampleSize) / (hitInPop / popSize)
    #   mat.count <- matrix(c(n.inset,n.genes.sset - n.inset,  n.genomewide ,n.genes - n.genomewide), ncol = 2, byrow = FALSE)
    #   p.val <- fisher.test(mat.count)$p.value
    df.gene_function_enrichment$p.val[j] <- p.val
    df.gene_function_enrichment$n.genes[j] = hitInSample
    df.gene_function_enrichment$foldchange[j] = foldChange
    
  }
  
  # integrate multiple hypothesis testing
  df.gene_function_enrichment <- df.gene_function_enrichment[order(df.gene_function_enrichment$p.val),]
  df.gene_function_enrichment <- subset(df.BP_enrichment, df.BP_enrichment$p.val <= 0.05)
  
  
  head(df.bQTL_gene_partitioning)
  
  
  
  
  
  
  
}

write.csv(df.bQTL_gene_partitioning, "output/df.bQTL_gene_partitioning_ChipSeq.csv", row.names = FALSE)

test = read.csv("output/df.bQTL_gene_partitioning_ChipSeq.csv")

df.bQTL_RNAseq

length(unique(test$gene.ID))


length(intersect())










postTotal.significant.i <- subset(postTotal.significant, postTotal.significant[,idx.chr] == i)

postTotal.significant.i["promoter_5kb"] <- "no"
postTotal.significant.i["promoter_1kb"] <- "no"
postTotal.significant.i["gene"] <- "no"
postTotal.significant.i["five_prime_UTR"] <- "no"
postTotal.significant.i["exon"] <- "no"
postTotal.significant.i["intron"] <- "no"
postTotal.significant.i["three_prime_UTR"] <- "no"
postTotal.significant.i["post_gene_1kb"] <- "no"

postTotal.significant.i["gene.ID"] <- NA

df.gff.i <- subset(df.gene_annotation, df.gene_annotation$chr == i)

l.df.gff.i <- vector(mode = "list", length = length(v.genePartitions))
l.df.gff.i[[1]] <- subset(df.gff.i, df.gff.i$partition %in% v.genePartitions[1])
l.df.gff.i[[2]] <- subset(df.gff.i, df.gff.i$partition %in% v.genePartitions[2])
l.df.gff.i[[3]] <- subset(df.gff.i, df.gff.i$partition %in% v.genePartitions[3])
l.df.gff.i[[4]] <- subset(df.gff.i, df.gff.i$partition %in% v.genePartitions[4])
l.df.gff.i[[5]] <- subset(df.gff.i, df.gff.i$partition %in% v.genePartitions[5])

# preselect all in gene 
pb <- txtProgressBar(min = 0, max = nrow(postTotal.significant.i), style = 3)

for(j in 1:nrow(postTotal.significant.i)){
  
  setTxtProgressBar(pb, j)
  
  # in gene 
  i.set <- which(l.df.gff.i[[1]]$pos.start < postTotal.significant.i[j,idx.pos] & postTotal.significant.i[j,idx.pos] < l.df.gff.i[[1]]$pos.stop)
  






  
  





