
gene_function_enrichment <- function(df.bQTL_gene_partitioning, th.pval = 0.05){

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
  df.gene_function_enrichment <- subset(df.BP_enrichment, df.BP_enrichment$p.val <= th.pval)
  
  return(df.gene_function_enrichment)
}