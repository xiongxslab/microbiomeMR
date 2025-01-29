##########Functional annotation GO
enrich.go <- enrichGO(gene = matched_gene,  #list fot enrichment    
                      OrgDb = 'org.Hs.eg.db',  #gene database
                      keyType = 'SYMBOL',  
                      ont = 'ALL',  #GO Ontology， BP、MF、CC
                      pAdjustMethod = "fdr",  #
                      universe = Background_1,
                      pvalueCutoff = 1,  # p value（1 for outputting all）
                      qvalueCutoff = 1,  #q value（1 for outputting all）
                      readable = FALSE)
write.table(enrich.go, file = print(paste0("/data/slurm/wanghc/microbiome_QTL/KEGG_and_GO_analysis/specific_regulate/Tissue-specific/GO_newbackground_OR/",microbiome_name,"-",category_name,".go.txt")), quote = FALSE, sep = "\t", row.names = FALSE)
enrich.go.cutoff <- enrich.go %>% as.data.frame() %>% filter(pvalue< 0.05)

calculate_overlap <- function(genes1, genes2) {
  # Split the geneID strings into vectors of individual gene IDs
  genes1 <- unlist(strsplit(genes1, "/"))
  genes2 <- unlist(strsplit(genes2, "/"))
  
  # Calculate the number of overlapping genes
  overlap <- length(intersect(genes1, genes2))
  
  # Determine the smaller number of genes between the two sets
  min_genes <- min(length(genes1), length(genes2))
  
  # Calculate the overlap percentage
  overlap_percentage <- (overlap / min_genes) * 100
  
  return(overlap_percentage)
}

rows_to_keep <- rep(TRUE, nrow(enrich.go.cutoff))

for (i in 2:nrow(enrich.go.cutoff)) {
  for (j in 1:(i - 1)) {
    if (rows_to_keep[j]) {
      overlap_percentage <- calculate_overlap(enrich.go.cutoff$geneID[j], enrich.go.cutoff$geneID[i])
      if (overlap_percentage > 50) {
        if (enrich.go.cutoff$pvalue[i] < enrich.go.cutoff$pvalue[j]) {
          rows_to_keep[j] <- FALSE
        } else {
          rows_to_keep[i] <- FALSE
        }
        break
      }
    }
  }
}

filtered_df <-enrich.go.cutoff[rows_to_keep, ]
write.table(filtered_df, file = print(paste0("/data/slurm/wanghc/microbiome_QTL/KEGG_and_GO_analysis/specific_regulate/Tissue-specific/GO_newbackground_OR/",microbiome_name,"-",category_name,".go.cutoff.txt")), quote = FALSE, sep = "\t", row.names = FALSE)









