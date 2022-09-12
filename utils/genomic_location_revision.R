add_genomic_location <- function(df.bQTLs, 
                                 chr = "B73-chr",
                                 pos = "B73-pos",
                                 df.gene_annotation,
                                 n.cpus = 1){
  
  v.chromosomes <- unique(df.bQTLs[,chr])
  n.chromosomes <- length(v.chromosomes)

  strt<-Sys.time() 
  cl<-makeCluster(min(n.chromosomes, n.cpus))
  registerDoParallel(cl)
  
  l.bQTL_genomic_location_per_chromosome <- foreach(i = 1:n.chromosomes, .packages=c("seqinr", "VariantAnnotation", "Biostrings")) %dopar% { 
    
    chr.i <- v.chromosomes[i]
    
    df.bQTLs.i <- subset(df.bQTLs, df.bQTLs[,chr] == chr.i)
    
    df.bQTLs.i["promoter_5kb"] <- "no"
    df.bQTLs.i["promoter_1kb"] <- "no"
    df.bQTLs.i["gene"] <- "no"
    df.bQTLs.i["five_prime_UTR"] <- "no"
    df.bQTLs.i["exon"] <- "no"
    df.bQTLs.i["intron"] <- "no"
    df.bQTLs.i["three_prime_UTR"] <- "no"
    df.bQTLs.i["post_gene_1kb"] <- "no"
    df.bQTLs.i["non_genic"] <- "no"
    
    df.bQTLs.i["gene.ID"] <- NA
    df.bQTLs.i["distance_to_gene"] <- NA
    df.bQTLs.i["strand"] <- NA
    
    # case of double annotation - only affects downstream    
    df.bQTLs.i["gene.ID.additional_downstream"] <- NA
    df.bQTLs.i["distance_to_gene.additional_downstream"] <- NA
    df.bQTLs.i["strand.additional_downstream"] <- NA
    
    df.annot.chr <- subset(df.gene_annotation, df.gene_annotation$chr == chr.i)
    
    df.gene <- subset(df.annot.chr, df.annot.chr$partition %in% "gene")
    df.5utr <- subset(df.annot.chr, df.annot.chr$partition %in% "five_prime_UTR")
    df.3utr <- subset(df.annot.chr, df.annot.chr$partition %in% "three_prime_UTR")
    df.exon <- subset(df.annot.chr, df.annot.chr$partition %in% "exon")
    
  
    # non strand sensitive part - in gene
    pb <- txtProgressBar(min = 0, max = nrow(df.bQTLs.i), style = 3)
    for(j in 1:nrow(df.bQTLs.i)){
      setTxtProgressBar(pb, j)

      # in gene
      i.gene <- which(df.gene$pos.start < df.bQTLs.i[j,pos] & df.bQTLs.i[j,pos] < df.gene$pos.stop)
      
      if(length(i.gene) > 0){
        
        df.bQTLs.i$gene[j] <- "yes"
        df.bQTLs.i$gene.ID[j] <- df.gene$gene.ID[i.gene[1]]
        df.bQTLs.i$strand[j] <- df.gene$strand[i.gene[1]]
        df.bQTLs.i$distance_to_gene[j] <- 0
        
        # in 5 prime utr
        i.set <- which(df.5utr$pos.start < df.bQTLs.i[j,pos] & df.bQTLs.i[j,pos] < df.5utr$pos.stop)
        if(length(i.set) > 0){
          df.bQTLs.i$five_prime_UTR[j] <- "yes"
        }else{
          # in 3 prime utr
          i.set <- which(df.3utr$pos.start < df.bQTLs.i[j,pos] & df.bQTLs.i[j,pos] < df.3utr$pos.stop)
          if(length(i.set) > 0){
            df.bQTLs.i$three_prime_UTR[j] <- "yes"
          }else{
            # in CDS (exon)
            i.set <- which(df.exon$pos.start < df.bQTLs.i[j,pos] & df.bQTLs.i[j,pos] < df.exon$pos.stop)
            if(length(i.set) > 0){
              df.bQTLs.i$exon[j] <- "yes"
            }else{
              df.bQTLs.i$intron[j] <- "yes"
            }
          } 
        }
      }
    }
    close(pb)
    
    # strand sensitive part #
    
    # identify bQTL in 5kb and 1kb distance 
    df.gene["up_5kb"] <- 0
    df.gene["up_1kb"] <- 0
    df.gene["start_gene"] <- 0
    df.gene["end_gene"] <- 0
    df.gene["down_1kb"] <- 0
    
    for(k in 1:nrow(df.gene)){
      
      # strand sensitive 
      if(df.gene$strand[k] == "+"){  # promoter before (-) the gene
        
        df.gene$up_5kb[k] <- as.numeric(df.gene$pos.start[k] - 5000)
        df.gene$up_1kb[k] <- as.numeric(df.gene$pos.start[k] - 1000)
        df.gene$start_gene[k] <- as.numeric(df.gene$pos.start[k])
        
        df.gene$end_gene[k] <- as.numeric(df.gene$pos.stop[k])
        df.gene$down_1kb[k] <- as.numeric(df.gene$pos.stop[k] + 1000)
        
      }else if(df.gene$strand[k] == "-"){  # promoter after (+) the gene - perspective from positive strand 
        
        df.gene$up_1kb[k] <- as.numeric(df.gene$pos.stop[k] + 1000)
        df.gene$up_5kb[k] <- as.numeric(df.gene$pos.stop[k] + 5000)
        df.gene$start_gene[k] <- as.numeric(df.gene$pos.stop[k])
        
        df.gene$end_gene[k] <- as.numeric(df.gene$pos.start[k])
        df.gene$down_1kb[k] <- as.numeric(df.gene$pos.start[k] - 1000)

      }
    }
    
    
    # search only for genes?
    i.env <- which(df.bQTLs.i$gene == "no")
    
    for(j in i.env){

      # in 1 KB upstream - both strands 
      i.set <- which(df.gene$up_1kb < df.bQTLs.i[j, pos] & df.bQTLs.i[j, pos] < df.gene$start_gene 
                     | df.gene$start_gene < df.bQTLs.i[j, pos] & df.bQTLs.i[j, pos] < df.gene$up_1kb)
      
      if(length(i.set) > 0){
        
        df.bQTLs.i$promoter_1kb[j] <- "yes"
        df.bQTLs.i$gene.ID[j] <- df.gene$gene.ID[i.set[1]]
        df.bQTLs.i$strand[j] <- df.gene$strand[i.set[1]]
        
        if(df.bQTLs.i$strand[j] == "+"){
          df.bQTLs.i$distance_to_gene[j] <- df.gene$pos.start[i.set[1]] - df.bQTLs.i[j, pos]
        }else{
          df.bQTLs.i$distance_to_gene[j] <- df.bQTLs.i[j, pos] - df.gene$pos.stop[i.set[1]]
        }
        
      }else{
        
        # in 5 KB upstream - consider both strands
        i.set <- which(df.gene$up_5kb < df.bQTLs.i[j, pos] & df.bQTLs.i[j,pos] < df.gene$start_gene
                       | df.gene$start_gene < df.bQTLs.i[j, pos] & df.bQTLs.i[j, pos] < df.gene$up_5kb)
        
        if(length(i.set)  > 0){

          df.bQTLs.i$promoter_5kb[j] <- "yes"  
          df.bQTLs.i$gene.ID[j] <- df.gene$gene.ID[i.set[1]]
          df.bQTLs.i$strand[j] <- df.gene$strand[i.set[1]]
          
          if(df.bQTLs.i$strand[j] == "+"){
            df.bQTLs.i$distance_to_gene[j] <- df.gene$pos.start[i.set[1]] - df.bQTLs.i[j, pos]
          }else{
            df.bQTLs.i$distance_to_gene[j] <- df.bQTLs.i[j, pos] - df.gene$pos.stop[i.set[1]]
          }
          
        }
        
      }
      
      # now also double partition assignments (additional downstream 1kb)
      qtl.exist <- any(df.bQTLs.i[j, c("promoter_5kb", "promoter_1kb", "gene", "five_prime_UTR", "exon", "intron", "three_prime_UTR")] == "yes")

      if(!qtl.exist){
        col.gene_id <- "gene.ID"
        col.distance <- "distance_to_gene"
        col.strand <- "strand"
      }else{
        col.gene_id <- "gene.ID.additional_downstream"
        col.distance <- "distance_to_gene.additional_downstream"
        col.strand <- "strand.additional_downstream"
      }

      i.set <- which(df.gene$end_gene < df.bQTLs.i[j, pos] & df.bQTLs.i[j, pos] < df.gene$down_1kb
                     | df.gene$down_1kb < df.bQTLs.i[j, pos] & df.bQTLs.i[j, pos] < df.gene$end_gene)
      
      if(length(i.set)  > 0){
          
        df.bQTLs.i$post_gene_1kb[j] <- "yes"
          
        df.bQTLs.i[j, col.gene_id] <- df.gene$gene.ID[i.set[1]]
        df.bQTLs.i[j, col.strand] <- df.gene$strand[i.set[1]]
        
        if(df.bQTLs.i[j, col.strand] == "-"){ # inverse because post gene
          df.bQTLs.i[j, col.distance] <- df.gene$pos.start[i.set[1]] - df.bQTLs.i[j, pos]
        }else{
          df.bQTLs.i[j, col.distance] <- df.bQTLs.i[j, pos] - df.gene$pos.stop[i.set[1]]
        }
        
      }
      
    }

    df.bQTLs.i
    
  }
  
  stopCluster(cl)
  print(Sys.time()-strt)
  
  
  df.bQTL_genomic_location <- c()
  for(i in 1:n.chromosomes){
    df.bQTL_genomic_location <- rbind(df.bQTL_genomic_location, l.bQTL_genomic_location_per_chromosome[[i]])
  }
  
  # if not found anywhere near or in gene
  i.set <- apply(df.bQTL_genomic_location, 1, function(m) {all(m[c("promoter_5kb", "promoter_1kb", "gene", "five_prime_UTR", "exon", "intron", "three_prime_UTR", "post_gene_1kb", "non_genic")] == "no")})
  i.set <- which(i.set == TRUE)
  df.bQTL_genomic_location$non_genic[i.set] <- "yes"
  
  return(df.bQTL_genomic_location)
}



add_genomic_location_bQTLs <- function(df.bQTLs, 
                                       df.gene_annotation,
                                       n.cpus = 1,
                                       b.low_affinity = F,
                                       b.B73_only = F){
  
  
  if(b.B73_only){
    print("no separation by postfrequency - using only B73 coordinates for bQTLs and gene annotation")
    df.bQTL_genomic_location <- add_genomic_location(df.bQTLs, 
                                                     chr = "B73-chr",
                                                     pos = "B73-pos",
                                                     df.gene_annotation,
                                                     n.cpus = n.cpus)
    
  }else{
    # added >= 0.5 postfrequency for genomic location to include 0.5 allelic bias SNPs (will not affect ASBs anyways)
    df.bQTLs.B73 <- subset(df.bQTLs, df.bQTLs$POSTfreq >= 0.5)
    df.bQTLs.Mo17 <- subset(df.bQTLs, df.bQTLs$POSTfreq < 0.5)
    
    if(!b.low_affinity){
      df.bQTLs.affinity <- df.bQTLs.B73
    }else{
      df.bQTLs.affinity <- df.bQTLs.Mo17
    }
    
    
    df.bQTL_genomic_location <- add_genomic_location(df.bQTLs.affinity,
                                                     chr = "B73-chr",
                                                     pos = "B73-pos",
                                                     df.gene_annotation,
                                                     n.cpus = n.cpus)
    
    if(!b.low_affinity){
      df.bQTLs.affinity <- df.bQTLs.Mo17
    }else{
      df.bQTLs.affinity <- df.bQTLs.B73
    }
    
    df.bQTL_genomic_location <- rbind(df.bQTL_genomic_location, 
                                      add_genomic_location(df.bQTLs.affinity,
                                                           chr = "Mo17-chr",
                                                           pos = "Mo17-pos",
                                                           df.gene_annotation,
                                                           n.cpus = n.cpus))
  }
  return(df.bQTL_genomic_location)
}


# 
# get_partitions <- function(df.bQTL_genomic_location,
#                            v.genePartitions = c("gene", "five_prime_UTR",  "CDS", "three_prime_UTR", "exon")){
#   
#   for(j in 1:length(v.partitions)){
#     print(v.partitions[j])
#     tmp <- subset(df.bQTL_genomic_location, df.bQTL_genomic_location[,v.partitions[j]] == "yes")
#     print(mean(abs(tmp$POSTfreq - 0.5)))
#   }
#   
#   # numbers - summary
#   df.partitions[1,s+1] <- table(df.bQTL_genomic_location$promoter_5kb)[2]
#   df.partitions[2,s+1] <- table(df.bQTL_genomic_location$promoter_1kb)[2]
#   df.partitions[3,s+1] <- table(df.bQTL_genomic_location$gene)[2]
#   df.partitions[4,s+1] <- table(df.bQTL_genomic_location$five_prime_UTR)[2]
#   df.partitions[5,s+1] <- table(df.bQTL_genomic_location$exon)[2]
#   df.partitions[6,s+1] <- table(df.bQTL_genomic_location$intron)[2]
#   df.partitions[7,s+1] <- table(df.bQTL_genomic_location$three_prime_UTR)[2]
#   df.partitions[8,s+1] <- table(df.bQTL_genomic_location$post_gene_1kb)[2]
#   # df.partitions[9,s+1] <- table(df.bQTL_gene_partitioning$post_gene_5kb)[2]
#   
#   # adapted - no post 5 kb
#   df.partitions[9,s+1] <- nrow(df.bQTL_genomic_location) - sum(df.partitions[1:8,s+1])
#   df.partitions[10,s+1] <- nrow(df.bQTL_genomic_location)
#   
#   df.partitions["percentage_significant"] <- df.partitions[,2] / df.partitions[10,2]
#   df.partitions["percentage_non_significant"] <- df.partitions[,3] / df.partitions[10,3]
#   
#   return(df.partitions)
# }
# 

