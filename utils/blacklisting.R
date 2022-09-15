
read_invert_rle <- function(gr.reads, start, end, rm.zero = F){
  
  # remove 0 scores (no reads) and inverse RLE
  df.reads <- as.data.frame(gr.reads)
  
  if(rm.zero)
    df.reads <- subset(df.reads, df.reads$score > 0)
  
  idx <- 1
  v.pos <- numeric(sum(df.reads$width))
  for(j in 1:nrow(df.reads)){
    df.reads$start[j] <- max(start, df.reads$start[j])
    df.reads$end[j] <- min(df.reads$end[j], end)
    df.reads$width[j] <- df.reads$end[j] - df.reads$start[j] + 1
  }
  
  idx <- 1
  v.pos <- numeric(sum(df.reads$width))
  for(j in 1:nrow(df.reads)){
    count <- idx + df.reads$width[j]
    v.pos[idx:(count - 1)] <- df.reads$start[j]:df.reads$end[j]
    idx <- count
  }
  
  v.scores <- inverse.rle(list(values = df.reads$score, lengths = df.reads$width))
  df.reads.irle <- data.frame(pos = v.pos, score = v.scores)
  
  return(df.reads.irle)
}


define_split <- function(v.snps.pos,
                         chr,
                         gr.reads, 
                         th.bp_ranges){
  
  p.start <- min(v.snps.pos) - th.bp_ranges
  p.end <- max(v.snps.pos) + th.bp_ranges
  
  q=GRanges(seqnames=chr,
            ranges=IRanges(start = p.start, end = p.end),
            strand="*")
  
  return(subsetByOverlaps(gr.reads, q))
  
}


estimate_avg_region_postfreq <- function(pos,
                                         chr,
                                         gr.reads, 
                                         th.bp_ranges){
                  
    start <- pos - th.bp_ranges
    end <- pos + th.bp_ranges
    
    q=GRanges(seqnames=chr,
              ranges=IRanges(start = start, end = end),
              strand="*")
    
    return(mean(score(subsetByOverlaps(gr.reads, q))))

}


partial_estimate <- function(df.snps.partial, 
                             gr.reads.B73.partial,
                             gr.reads.Mo17.partial,
                             th.bp_ranges,
                             n.snps.split){

  v.snp.pos.B73 <- df.snps.partial$`B73-pos`
  v.snps.pos.Mo17 <- df.snps.partial$`Mo17-pos`
  
  chr.B73 <- df.snps.partial$`B73-chr`[1]
  chr.Mo17 <- df.snps.partial$`Mo17-chr`[1]
  
  v.snps.pos <- v.snp.pos.B73[1:min(length(v.snp.pos.B73), n.snps.split)]
  gr.reads.B73.split <- define_split(v.snps.pos, chr.B73, gr.reads.B73.partial, th.bp_ranges)
  
  v.snps.pos <- v.snps.pos.Mo17[1:min(length(v.snps.pos.Mo17), n.snps.split)]
  gr.reads.Mo17.split <- define_split(v.snps.pos, chr.Mo17, gr.reads.Mo17.partial, th.bp_ranges)
  
  message("estimate region based postfrequency ")
  
  v.region.post_frequency <- numeric(length = nrow(df.snps.partial))
  names(v.region.post_frequency) <- df.snps.partial$id
  pb <- txtProgressBar(min = 0, max = length(v.region.post_frequency), style = 3)
  for(j in 1:length(v.snp.pos.B73)){
    setTxtProgressBar(pb, j)
    
    if(j %% n.snps.split == 0){

      gr.reads.B73.split <- define_split(v.snp.pos.B73[j:min(length(v.snp.pos.B73), (j + n.snps.split))], 
                                         chr.B73, gr.reads.B73.partial, th.bp_ranges)

      gr.reads.Mo17.split <- define_split(v.snps.pos.Mo17[j:min(length(v.snps.pos.Mo17), (j + n.snps.split))], 
                                          chr.Mo17, gr.reads.Mo17.partial, th.bp_ranges)
      
    }
    
    mu.B73 <- estimate_avg_region_postfreq(v.snp.pos.B73[j],
                                           chr.B73,
                                           gr.reads.B73.split, 
                                           th.bp_ranges)
    
    mu.Mo17 <- estimate_avg_region_postfreq(v.snps.pos.Mo17[j],
                                            chr.Mo17,
                                            gr.reads.Mo17.split, 
                                            th.bp_ranges)
    
    v.region.post_frequency[j] <- mu.B73 / (mu.B73 + mu.Mo17)
    
  }  
  close(pb)
  
  return(v.region.post_frequency)
}


blacklisting_background_distribution <- function(df.snps, 
                                                 path.raw_input_reads,
                                                 n.chromosomes = 10,
                                                 n.snps.split = 100,
                                                 th.bp_ranges = 500,
                                                 n.cpus = 2){
                
  # blacklist regions / snps of systematic bias
  # use B73/Mo17 coordinates from complete SNP data for correspondence in hybrid genome
  # use all snps to create background distribution (to check position of regions surrounding ASBs)
  gr.reads <- import.bw(path.raw_input_reads)
  
  # constraint for speed - significant distribution is enough
  
  # to speed up blacklisting use samples from B73 / Mo17 chr 1 - 10 only
  df.snps <- subset(df.snps, df.snps$`Mo17-chr` %in% paste("Mo17-chr", 1:10, sep = ""))
  df.snps["id"] <- seq(1, nrow(df.snps))
  df.snps["B73.norm_reads"] <- 0
  df.snps["Mo17.norm_reads"] <- 0

  # subset selection chr == chr
  df.snps.selection <- c()
  for(i in 1:n.chromosomes){
    
    print(i)
    
    gc()

    chr.B73 <- paste("B73-chr", i, sep = "")
    chr.Mo17 <- paste("Mo17-chr", i, sep = "")
    ## for computational efficiency only - subet to chrB73 equals chr Mo17 (should be enough background samples still)
    df.snps.i <- subset(df.snps, df.snps$`B73-chr` == chr.B73 & df.snps$`Mo17-chr` == chr.Mo17)
  
    v.snp.pos.B73 <- df.snps.i$`B73-pos`
    v.snp.pos.Mo17 <- df.snps.i$`Mo17-pos`
    
    print("B73")
    df.reads.irle <- read_invert_rle(chr.B73, v.snp.pos.B73 , gr.reads)
    gc()

    df.reads.irle <- subset(df.reads.irle, df.reads.irle$pos %in% v.snp.pos.B73)
    rownames(df.snps.i) <- as.character(v.snp.pos.B73)
    df.snps.i[as.character(df.reads.irle$pos), "B73.norm_reads"] <- df.reads.irle$score
    df.snps.i <- df.snps.i[as.character(v.snp.pos.B73), ]
    gc()
    
    print("Mo17")
    df.reads.irle <- read_invert_rle(chr.Mo17, v.snp.pos.Mo17, gr.reads)
    gc()

    df.reads.irle <- subset(df.reads.irle, df.reads.irle$pos %in% v.snp.pos.Mo17)
    rownames(df.snps.i) <- as.character(v.snp.pos.Mo17)
    df.snps.i[as.character(df.reads.irle$pos), "Mo17.norm_reads"] <- df.reads.irle$score
    df.snps.i <- df.snps.i[as.character(v.snp.pos.Mo17), ]
    gc()
    
    
    print("bind")
    df.snps.selection <- rbind(df.snps.selection, df.snps.i)
    
  } 
    
  # saveRDS(df.snps.selection, "df.snps.selection.rds")
    
  
  df.snps.B73 <- subset(df.snps.B73, df.snps.B73)
  
  for(i in 1:n.chromosomes){
    
    print(i)
    
    gc()
    
    chr.B73 <- paste("B73-chr", i, sep = "")
    chr.Mo17 <- paste("Mo17-chr", i, sep = "")
    ## for computational efficiency only - subet to chrB73 equals chr Mo17 (should be enough background samples still)
    df.snps.i <- subset(df.snps, df.snps$`B73-chr` == chr.B73 & df.snps$`Mo17-chr` == chr.Mo17)
    
    df.snps.B73.i <- subset(df.snps.B73, df.snps.B73$`B73-chr` == chr.B73)
    df.snps.B73.i <- df.snps.B73.i[order(df.snps.B73.i$`B73-pos`),]
    
    
    v.diff.pos <- diff(df.snps.B73.i$`B73-pos`)
    
    
    
    
    chr.B73 <- paste("B73-chr", i, sep = "")
    chr.Mo17 <- paste("Mo17-chr", i, sep = "")
    ## for computational efficiency only - subet to chrB73 equals chr Mo17 (should be enough background samples still)
    df.snps.i <- subset(df.snps, df.snps$`B73-chr` == chr.B73 & df.snps$`Mo17-chr` == chr.Mo17)
    
    
    
  }
  
  
  plot(df.snps.i$B73.norm_reads[1:1000], df.snps.i$Mo17.norm_reads[1:1000])
  
  x <- df.snps.i$B73.norm_reads[1:10000]
  y <- df.snps.i$Mo17.norm_reads[1:10000]
  
  pf <- df.snps.i$B73.norm_reads / (df.snps.i$B73.norm_reads + df.snps.i$Mo17.norm_reads)
  
  plot(df.snps.i$`B73-pos`, pf) # flag genomic ranges?
  
  library(MASS)
  den3d <- kde2d(x, y)
  
  # the new part:
  library(plotly)
  plot_ly(x=den3d$x, y=den3d$y, z=den3d$z) %>% add_surface()
  
  # v.diff.pos <- diff(df.snps.B73.i$`B73-pos`)
  
  
    ## Mo17 
    v.snps.chrs.Mo17 <- df.snps.B73.i$`Mo17-chr`
    v.snps.pos.Mo17 <- df.snps.B73.i$`Mo17-pos`
    
    v.chrs.Mo17 <- unique(df.snps.B73.i$`Mo17-chr`)
    l.chr.Mo17 <- vector(mode = "list", length = length(v.chrs.Mo17))
    names(l.chr.Mo17) <- v.chrs.Mo17
    message("prepare Mo17 chromosomal ranges")
    pb <- txtProgressBar(min = 0, max = length(v.chrs.Mo17), style = 3)
    for(j in 1:length(v.chrs.Mo17)){
      setTxtProgressBar(pb, j)
      df.snps.j <- subset(df.snps.B73.i, df.snps.B73.i$`Mo17-chr` == v.chrs.Mo17[j])
      q=GRanges(seqnames= v.chrs.Mo17[j],
                ranges=IRanges(start = min(df.snps.j$`Mo17-pos`), end = max(df.snps.j$`Mo17-pos`)),
                strand="*")
      
      l.chr.Mo17[[j]] <- subsetByOverlaps(gr.reads, q)
    }
    close(pb)
    gc()
    

    v.diff.pos <- diff(df.snps.B73.i$`B73-pos`)
    
    
    # divide into windows 
    # define windows 
    n.genomic_window <- 100000
    
    n.genomic_window # add until, set breaks... 
    sum(v.diff.pos[1:900]) # shift by ... 
    
    q=GRanges(seqnames=chr.B73,
              ranges=IRanges(start = df.snps.B73.i$`B73-pos`[1], end = df.snps.B73.i$`B73-pos`[900]),
              strand="*")
    
    test <-  score(subsetByOverlaps(gr.reads.B73, q))
    
    
    pb <- txtProgressBar(min = 0, max = nrow(df.snps.B73.i), style = 3)
    for(j in 1:nrow(df.snps.B73.i)){
      setTxtProgressBar(pb, j)
      
      pos <- df.snps.B73.i$`B73-pos`[j]
      
      df.snps.B73.i
      
      q=GRanges(seqnames=chr.B73,
                ranges=IRanges(start = pos, end = pos),
                strand="*")
      
      subsetByOverlaps(gr.reads.B73, q)
      
      
    }
    close(pb)
    
  }
  
  
  
bQTL_scatterplot <- function(df.bQTLs, n.chromosomes = 10){
  
  # TODO: general plot 
  # lines
  # per chromosome 
  
  
  
  ### scatter plot of the postfrequencies - should be updated to > , < 0.5
  v.offset <- numeric(n.chromosomes)
  v.start <- numeric(n.chromosomes)
  
  for(i in 1:n.chromosomes){
    chr <- paste("B73-chr", i, sep = "")
    df.bQTLs.i <- subset(df.bQTLs, df.bQTLs$`B73-chr` == chr)
    v.offset[i] <-  max(df.bQTLs.i$`B73-pos`)
    v.start[i] <- sum(v.offset[1:i])
  }
  
  v.start <- c(0, v.start)
  
  # directionality plot 
  df.scatterplot <- subset(df.bQTLs, df.bQTLs$`B73-chr` == "B73-chr1")
  
  for(i in 1:(n.chromosomes - 1)){
    chr <- paste("B73-chr", (i+1), sep = "")
    df.bQTLs.i <- subset(df.bQTLs, df.bQTLs$`B73-chr` == chr)
    df.bQTLs.i$`B73-pos` <- df.bQTLs.i$`B73-pos` + sum(v.offset[1:i])
    df.scatterplot <- rbind(df.scatterplot, df.bQTLs.i)
  }
  
  p6 <- ggplot(df.scatterplot, aes(x = `B73-pos`, y = POSTfreq, fill = POSTfreq)) +
    geom_point(shape = 21, size = 1,  alpha = 0.5, stroke = 0.0) + theme_bw() 
  
  p6 <- p6 + scale_fill_continuous(low = "blue", high = "green2")
  p6 <- p6 + geom_vline(xintercept = v.start, col='black', lwd=0.5, linetype="dashed")
  
  p6
}




 
# 
#   
#   message("Processing B73 biased SNPs")
#   df.snps.B73 <- subset(df.snps, df.snps$POSTfreq > 0.5)
#   for(i in 1:n.chromosomes){
#     
#     chr.B73 <- paste("B73-chr", i, sep = "")
#     chr.Mo17 <- paste("Mo17-chr", i, sep = "")
#     
#     gr.reads.B73 <- gr.reads[seqnames(gr.reads) == chr.B73]
#     gr.reads.Mo17 <- gr.reads[seqnames(gr.reads) == chr.Mo17]
#     
#     message(chr.B73)
#     
#     df.snps.B73.i <- subset(df.snps.B73, df.snps.B73$`B73-chr` == chr.B73 & 
#                                          df.snps.B73$`Mo17-chr` == chr.Mo17)
#     df.snps.B73.i <- df.snps.B73.i[order(df.snps.B73.i$`B73-pos`),]
#     
#     
#     i.sets <- round(nrow(df.snps.B73.i) / 5000, 0)
#   
#     strt<-Sys.time() 
#     cl<-makeCluster(min(n.chromosomes, n.cpus))
#     registerDoParallel(cl)
#     
#     l.region.post_frequency <- foreach(i = 1:n.chromosomes) %dopar% {     
#   
#     for(j in 1:5000){
#       
#         print(j)
#       
#         df.snps.partial <- df.snps.B73.i[1:i.sets,]
#         
#         gr.reads.B73.partial <- define_split(df.snps.partial$`B73-pos`, chr.B73, gr.reads.B73, th.bp_ranges)
#         gr.reads.Mo17.partial <- define_split(df.snps.partial$`Mo17-pos`, chr.Mo17, gr.reads.Mo17, th.bp_ranges)
#         
#         v.region.post_frequency <- partial_estimate(df.snps.partial, 
#                                                      gr.reads.B73.partial,
#                                                      gr.reads.Mo17.partial,
#                                                      th.bp_ranges,
#                                                      n.snps.split = 50)
#             
#     }
#     
#     stopCluster(cl)
#     print(Sys.time()-strt)
#                         
#     gc()                            
#     
#     
#     v.snp.pos.B73 <- df.snps.B73.i$`B73-pos`
#     v.snps.pos.Mo17 <- df.snps.B73.i$`Mo17-pos`
#     
#     
#     v.snps.pos <- v.snp.pos.B73[1:min(length(v.snp.pos.B73), n.snps.split)]
#     gr.reads.B73.split <- define_split(v.snps.pos, chr.B73, gr.reads.B73, th.bp_ranges)
# 
#     v.snps.pos <- v.snps.pos.Mo17[1:min(length(v.snps.pos.Mo17), n.snps.split)]
#     gr.reads.Mo17.split <- define_split(v.snps.pos, chr.Mo17, gr.reads.Mo17, th.bp_ranges)
# 
#     message("estimate region based postfrequency ")
#     
#     v.region.post_frequency <- numeric(length = length(v.snp.pos.B73))
#     pb <- txtProgressBar(min = 0, max = length(v.region.post_frequency), style = 3)
#     for(j in 1:length(v.snp.pos.B73)){
#       setTxtProgressBar(pb, j)
#       
#       print(j)
#       
#       if(j %% n.snps.split == 0){
#         
#         v.snps.pos <- v.snp.pos.B73[j:min(length(v.snp.pos.B73), (j + n.snps.split))]
#         gr.reads.B73.split <- define_split(v.snps.pos, chr.B73, gr.reads.B73, th.bp_ranges)
#         
#         v.snps.pos <- v.snps.pos.Mo17[j:min(length(v.snps.pos.Mo17), (j + n.snps.split))]
#         gr.reads.Mo17.split <- define_split(v.snps.pos, chr.Mo17, gr.reads.Mo17, th.bp_ranges)
# 
#       }
#       
#       mu.B73 <- estimate_avg_region_postfreq(v.snp.pos.B73[j],
#                                              chr.B73,
#                                              gr.reads.B73.split, 
#                                              th.bp_ranges)
#       
#       mu.Mo17 <- estimate_avg_region_postfreq(v.snps.pos.Mo17[j],
#                                               chr.Mo17,
#                                               gr.reads.Mo17.split, 
#                                               th.bp_ranges)
#       
#       v.region.post_frequency[j] <- mu.B73 / (mu.B73 + mu.Mo17)
#       
#     }  
#     close(pb)
#     
#     
#     
#   }
#   
#   
#   message("Processing Mo17 biased SNPs")
#   df.snps.Mo17 <- subset(df.snps, df.snps$POSTfreq < 0.5)
#   for(i in 1:n.chromosomes){
#     
#     chr.B73 <- paste("B73-chr", i, sep = "")
#     chr.Mo17 <- paste("Mo17-chr", i, sep = "")
#     
#     gr.reads.B73 <- gr.reads[seqnames(gr.reads) == chr.B73]
#     gr.reads.Mo17 <- gr.reads[seqnames(gr.reads) == chr.Mo17]
#     
#     message(chr.Mo17)
#     
#     df.snps.Mo17.i <- subset(df.snps.Mo17, df.snps.Mo17$`B73-chr` == chr.B73 & 
#                              df.snps.Mo17$`Mo17-chr` == chr.Mo17)
#     
#     df.snps.Mo17.i <- df.snps.Mo17.i[order(df.snps.Mo17.i$`B73-pos`),]
#     v.snps.pos.Mo17 <- df.snps.B73.i$`Mo17-pos`
#     v.snps.pos <- v.snps.pos.Mo17[1:min(length(v.snps.pos.Mo17), n.snps.split)]
#     gr.reads.Mo17.split <- define_split(v.snps.pos, chr.Mo17, gr.reads.Mo17, th.bp_ranges)
#     
#     
#     df.snps.Mo17.i <- df.snps.Mo17.i[order(df.snps.Mo17.i$`B73-pos`),]
#     v.snp.pos.B73 <- df.snps.Mo17.i$`B73-pos`
#     v.snps.pos <- v.snp.pos.B73[1:min(length(v.snp.pos.B73), n.snps.split)]
#     gr.reads.B73.split <- define_split(v.snps.pos, chr.B73, gr.reads.B73, th.bp_ranges)
#     
#     
#     message("estimate region based postfrequency ")
#     
#     v.region.post_frequency <- numeric(length = length(v.snp.pos.B73))
#     pb <- txtProgressBar(min = 0, max = length(v.region.post_frequency), style = 3)
#     for(j in 1:length(v.snp.pos.B73)){
#       setTxtProgressBar(pb, j)
#       
#       if(j %% n.snps.split == 0){
#         
#         v.snps.pos <- v.snp.pos.B73[j:min(length(v.snp.pos.B73), (j + n.snps.split))]
#         gr.reads.B73.split <- define_split(v.snps.pos, chr.B73, gr.reads.B73, th.bp_ranges)
#         
#         ## Mo17 
#         
#         v.snps.pos <- v.snps.pos.Mo17[j:min(length(v.snps.pos.Mo17), (j + n.snps.split))]
#         gr.reads.Mo17.split <- define_split(v.snps.pos, chr.Mo17, gr.reads.Mo17, th.bp_ranges)
#         
#       }
#       
#       mu.B73 <- estimate_avg_region_postfreq(v.snp.pos.B73[j],
#                                              chr.B73,
#                                              gr.reads.B73.split, 
#                                              th.bp_ranges)
#       
#       
#       mu.Mo17 <- estimate_avg_region_postfreq(v.snps.pos.Mo17[j],
#                                               chr.Mo17,
#                                               gr.reads.Mo17.split, 
#                                               th.bp_ranges)
#       
#       v.region.post_frequency[j] <- mu.B73 / (mu.B73 + mu.Mo17)
#       
#     }  
#     close(pb)
#     
#   }
#   
#   # SNP ORDER 
#   
#   # PLOT CORRELATION 
#   
# }
# 
  
  
extract_avg_postfrequency <- function(parent.chr,
                                      parent.pos, 
                                      df.bQTLs,
                                      path.raw_input_reads, 
                                      n.snps.split, 
                                      th.bp_ranges){
  
  
  v.chromosomes <- unique(df.bQTLs[,parent.chr])
  
  print(table(df.bQTLs[,parent.chr]))
  
  v.region.post_frequency <- c()
  for(i in 1:length(v.chromosomes)){
    gc()
    
    chr = v.chromosomes[i]
    message("processing: ", chr)
    
    df.bQTLs.i <- subset(df.bQTLs, df.bQTLs[,parent.chr] == chr)
    df.bQTLs.i <- df.bQTLs.i[order(df.bQTLs.i[,parent.pos]),]
    v.bQTLs.pos <- df.bQTLs.i[,parent.pos]
    # gr.reads.chr <- gr.reads[seqnames(gr.reads) == chr]
    gr.chr <- GRanges(seqnames = chr, ranges = IRanges(start = min(v.bQTLs.pos), end = max(v.bQTLs.pos)))
    gr.reads.chr <- import.bw(path.raw_input_reads, which = gr.chr)
    
    v.seq <- v.bQTLs.pos[1:min(length(v.bQTLs.pos), n.snps.split)]
    
    gr.reads.split <- define_split(v.seq,
                                   chr,
                                   gr.reads.chr, 
                                   th.bp_ranges)
    
    message("estimate region based postfrequency ")
    v.pf.chr <- numeric(length = length(v.bQTLs.pos))
    names(v.pf.chr) <- df.bQTLs.i$id
    pb <- txtProgressBar(min = 0, max = length(v.pf.chr), style = 3)
    for(j in 1:length(v.bQTLs.pos)){
      setTxtProgressBar(pb, j)
      
      if(j %% n.snps.split == 0){
        j.end <- min(length(v.bQTLs.pos), (j + n.snps.split))
        v.seq <- v.bQTLs.pos[j:j.end]
        gr.reads.split <- define_split(v.seq,
                                       chr,
                                       gr.reads.chr, 
                                       th.bp_ranges)
      }  
      
      start <- v.bQTLs.pos[j] - th.bp_ranges
      end <- v.bQTLs.pos[j] + th.bp_ranges
      
      q=GRanges(seqnames=chr,
                ranges=IRanges(start = start, end = end),
                strand="*")
      
      gr.read.snp <- subsetByOverlaps(gr.reads.split, q)
      df.reads.snp <- read_invert_rle(gr.read.snp, start, end)
      v.pf.chr[j] <- mean(df.reads.snp$score)
      #v.pf.chr[j] <- mean(score(subsetByOverlaps(gr.reads.split, q)))
    }  
    close(pb)
    
    v.region.post_frequency <- c(v.region.post_frequency, v.pf.chr)
  }
  
  rm(gr.reads.chr)
  
  v.region.post_frequency
  
}

  
estimate_genomic_postfrequency <- function(df.bQTLs, 
                                           path.raw_input_reads,
                                            n.snps.split = 50,
                                            th.bp_ranges = 500){

  gc()
  
  
  df.bQTLs["id"] <- seq(1, nrow(df.bQTLs))

  message("Processing B73")
  
  parent.chr = "B73-chr"
  parent.pos = "B73-pos"
  
  v.region.post_frequency <- extract_avg_postfrequency(parent.chr,
                                                       parent.pos, 
                                                       df.bQTLs,
                                                       path.raw_input_reads, 
                                                       n.snps.split,
                                                       th.bp_ranges)
                  
  df.bQTLs["B73.avg.norm_reads"] <- 0
  df.bQTLs$B73.avg.norm_reads[as.numeric(names(v.region.post_frequency))] <- v.region.post_frequency
  gc()
  
  message("Processing Mo17")
  
  parent.chr = "Mo17-chr"
  parent.pos = "Mo17-pos"
  
  v.region.post_frequency <- extract_avg_postfrequency(parent.chr,
                                                       parent.pos, 
                                                       df.bQTLs,
                                                       path.raw_input_reads, 
                                                       n.snps.split,
                                                       th.bp_ranges)
  
  df.bQTLs["Mo17.avg.norm_reads"] <- 0
  df.bQTLs$Mo17.avg.norm_reads[as.numeric(names(v.region.post_frequency))] <- v.region.post_frequency
  gc()
  
  df.bQTLs["pf.region"] <- df.bQTLs$B73.avg.norm_reads / (df.bQTLs$B73.avg.norm_reads + df.bQTLs$Mo17.avg.norm_reads)
  
  plot(df.bQTLs$POSTfreq, df.bQTLs$pf.region)
  
 
  gc()
  
  return(df.bQTLs)
}


identify_systematic_bias <- function(df.bQTLs, 
                                     df.bgSNPs.candidates,
                                     path.raw_input_reads,
                                     n.snps.split = 50,
                                     th.bp_ranges = 1000,
                                     n.bg.multiplier = 1,
                                     bQTL_only = F){

  message("estimate post frequency on bQTLs in genomic region from input read data")

  df.bQTLs.bl <- estimate_genomic_postfrequency(df.bQTLs, 
                                                path.raw_input_reads,
                                                n.snps.split,
                                                th.bp_ranges)
  
  if(bQTL_only)
    return(list(df.bQTLs=df.bQTLs.bl))
  
  message("sample background post frequency distribution from read data input")
  
  source("utils/qtl_background_sampling.R")
  
  df.bgSNPs <- create_background_QTLs(df.bQTLs,
                                      df.bgSNPs.candidates, 
                                      n.bg.multiplier, 
                                      v.partitions = c(),
                                      seed = 1234)
  
  df.bgSNPs.bl <- estimate_genomic_postfrequency(df.bgSNPs, 
                                                 path.raw_input_reads,
                                                  n.snps.split,
                                                  th.bp_ranges)
            
  
  
  gc()
  
  return(list(df.bQTLs=df.bQTLs.bl, df.bgSNPs=df.bgSNPs.bl))
}


filter_systematic_bias <- function(df.bQTLs,
                                   v.bg.genomic_postfrequency,
                                   th.prob = 0.05,
                                   plot.trendline = T){
  
  ### apply filter - given probability cutoffs on the genomic postfrequency 
  
  th.lower <- quantile(v.bg.genomic_postfrequency, th.prob)
  th.upper  <- quantile(v.bg.genomic_postfrequency, 1 - th.prob)

  library(ggplot2)
  p <- ggplot(df.bQTLs,aes(POSTfreq, pf.region)) + 
    geom_point(alpha= 0.3) + theme_bw() + 
    annotate("rect", xmin=0, xmax=1,ymin=th.upper,ymax=Inf,fill="green",alpha=0.2) + 
    annotate("rect", xmin=0, xmax=1,ymin=0,ymax=th.lower,fill="blue",alpha=0.2)
    
  if(plot.trendline)
    p <- p + geom_smooth(method="lm", linetype = "dashed", col = "red")
    
  plot(p)
  
  plot(allelic_bias_per_chromosome(df.bQTLs, n.chromosomes = 10))

  df.bQTLs["outlier"] <- T
  idx <- which(th.lower <= df.bQTLs$pf.region & df.bQTLs$pf.region <= th.upper)
  df.bQTLs$outlier[idx] <- F

  plot(allelic_bias_per_chromosome_input(df.bQTLs, n.chromosomes = 10))
  
  df.bQTLs <- subset(df.bQTLs, df.bQTLs$outlier == F)
  
  plot(allelic_bias_per_chromosome_input(df.bQTLs, n.chromosomes = 10))
  
  gc()


  return(df.bQTLs)
}



allelic_bias_per_chromosome_input <- function(df.bQTLs, 
                                              n.chromosomes = 10){
  
  # if flag_outliers
  
  ### scatter plot of the postfrequencies - should be updated to > , < 0.5
  v.offset <- numeric(n.chromosomes)
  v.start <- numeric(n.chromosomes)
  
  for(i in 1:n.chromosomes){
    chr <- paste("B73-chr", i, sep = "")
    df.bQTLs.i <- subset(df.bQTLs, df.bQTLs$`B73-chr` == chr)
    v.offset[i] <-  max(df.bQTLs.i$`B73-pos`)
    v.start[i] <- sum(v.offset[1:i])
  }
  
  v.start <- c(0, v.start)
  
  # directionality plot 
  df.scatterplot <- subset(df.bQTLs, df.bQTLs$`B73-chr` == "B73-chr1")
  
  for(i in 1:(n.chromosomes - 1)){
    chr <- paste("B73-chr", (i+1), sep = "")
    df.bQTLs.i <- subset(df.bQTLs, df.bQTLs$`B73-chr` == chr)
    df.bQTLs.i$`B73-pos` <- df.bQTLs.i$`B73-pos` + sum(v.offset[1:i])
    df.scatterplot <- rbind(df.scatterplot, df.bQTLs.i)
  }
  
  df.bQTLs.outliers <- subset(df.scatterplot, df.scatterplot$outlier == T)
  #df.bQTLs.outliers <- gapminder %>% 
  #  filter(gdpPercap>=59000)
  
  p6 <- ggplot(df.scatterplot, aes(x = `B73-pos`, y = POSTfreq, fill = POSTfreq)) +
    geom_point(shape = 21, size = 1,  alpha = 0.5, stroke = 0.0) + theme_bw() # shape = 21
  
  p6 <- p6 + scale_fill_continuous(low = "blue", high = "green2")
  p6 <- p6 + geom_vline(xintercept = v.start, col='black', lwd=0.5, linetype="dashed")
  
  p6 <- p6 + geom_point(data = df.bQTLs.outliers, aes(x=`B73-pos`,y=POSTfreq), size = 0.7, color="red", alpha = 0.5)

  p6

}




allelic_bias_per_chromosome <- function(df.bQTLs, 
                                        n.chromosomes = 10){

  # if flag_outliers
  
  ### scatter plot of the postfrequencies - should be updated to > , < 0.5
  v.offset <- numeric(n.chromosomes)
  v.start <- numeric(n.chromosomes)
  
  for(i in 1:n.chromosomes){
    chr <- paste("B73-chr", i, sep = "")
    df.bQTLs.i <- subset(df.bQTLs, df.bQTLs$`B73-chr` == chr)
    v.offset[i] <-  max(df.bQTLs.i$`B73-pos`)
    v.start[i] <- sum(v.offset[1:i])
  }
  
  v.start <- c(0, v.start)
  
  # directionality plot 
  df.scatterplot <- subset(df.bQTLs, df.bQTLs$`B73-chr` == "B73-chr1")
  
  for(i in 1:(n.chromosomes - 1)){
    chr <- paste("B73-chr", (i+1), sep = "")
    df.bQTLs.i <- subset(df.bQTLs, df.bQTLs$`B73-chr` == chr)
    df.bQTLs.i$`B73-pos` <- df.bQTLs.i$`B73-pos` + sum(v.offset[1:i])
    df.scatterplot <- rbind(df.scatterplot, df.bQTLs.i)
  }
  
  p6 <- ggplot(df.scatterplot, aes(x = `B73-pos`, y = POSTfreq, fill = POSTfreq)) +
    geom_point(shape = 21, size = 1,  alpha = 0.5, stroke = 0.0) + theme_bw() 
  
  p6 <- p6 + scale_fill_continuous(low = "blue", high = "green2")
  p6 <- p6 + geom_vline(xintercept = v.start, col='black', lwd=0.5, linetype="dashed")
  
  p6
}


bQTL_scatterplot <- function(df.bQTLs, n.chromosomes = 10){
  
  ### scatter plot of the postfrequencies - should be updated to > , < 0.5
  v.offset <- numeric(n.chromosomes)
  v.start <- numeric(n.chromosomes)
  
  for(i in 1:n.chromosomes){
    chr <- paste("B73-chr", i, sep = "")
    df.bQTLs.i <- subset(df.bQTLs, df.bQTLs$`B73-chr` == chr)
    v.offset[i] <-  max(df.bQTLs.i$`B73-pos`)
    v.start[i] <- sum(v.offset[1:i])
  }
  
  v.start <- c(0, v.start)
  
  # directionality plot 
  df.scatterplot <- subset(df.bQTLs, df.bQTLs$`B73-chr` == "B73-chr1")
  
  for(i in 1:(n.chromosomes - 1)){
    chr <- paste("B73-chr", (i+1), sep = "")
    df.bQTLs.i <- subset(df.bQTLs, df.bQTLs$`B73-chr` == chr)
    df.bQTLs.i$`B73-pos` <- df.bQTLs.i$`B73-pos` + sum(v.offset[1:i])
    df.scatterplot <- rbind(df.scatterplot, df.bQTLs.i)
  }
  
  p6 <- ggplot(df.scatterplot, aes(x = `B73-pos`, y = POSTfreq, fill = POSTfreq)) +
    geom_point(shape = 21, size = 1,  alpha = 0.5, stroke = 0.0) + theme_bw() 
  
  p6 <- p6 + scale_fill_continuous(low = "blue", high = "green2")
  p6 <- p6 + geom_vline(xintercept = v.start, col='black', lwd=0.5, linetype="dashed")
  
  p6
}
