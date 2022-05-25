add_gene_features <- function(df.bQTL_gene_partitioning,
                              df.gene_function = NULL,
                              df.gene_orthologs = NULL,
                              df.gene_conversion.AGPv3_to_AGPv4 = NULL){
  
  
  # add the AGPv3 gene annotation
  if(exists("df.gene_conversion.AGPv3_to_AGPv4") && !is.null(df.gene_conversion.AGPv3_to_AGPv4)){
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
  if(exists("df.gene_function") && !is.null(df.gene_function)){
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
  if(exists("df.gene_orthologs") && !is.null(df.gene_orthologs)){
    df.bQTL_gene_partitioning["Arabidopsis_ortholog"] <- NA
    for(i in 1:nrow(df.bQTL_gene_partitioning)){
      idx.AGPv3 <- which(df.gene_orthologs$locusName == df.bQTL_gene_partitioning$gene.ID.AGPv3[i])
      if(length(idx.AGPv3) > 0){
        df.bQTL_gene_partitioning$Arabidopsis_ortholog[i] <-  df.gene_orthologs$Best.hit.arabi.name[idx.AGPv3[1]]
      }
    }
  }
  
  
  return(df.bQTL_gene_partitioning)
}