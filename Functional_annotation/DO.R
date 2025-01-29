output="/data/slurm/huyh/DO_result/s_and_i"
enrichgene <- fread(file="/data/slurm/huyh/DO/data/Tissue-specific/s_and_i_1.txt", sep = "\t", header = F)

BACKGROUND <- fread(file="/data/slurm/wanghc/microbiome_QTL/KEGG_and_GO_analysis/specific_regulate/Tissue-specific/Background.txt", sep = "\t", header = F)
background1=bitr(BACKGROUND$V2,fromType='SYMBOL',toType='ENTREZID',OrgDb=org.Hs.eg.db)
background2=background1$ENTREZID

geneid <- bitr(enrichgene$V2,fromType='SYMBOL',toType='ENTREZID',OrgDb=org.Hs.eg.db)# convert the geneid to the ENTREZID format
gene <- geneid$ENTREZID

result=enrichDO(gene=gene,
                ont="DO",
                pAdjustMethod = "fdr",
                universe=background2,
                readable=TRUE,
                minGSSize =5,
                pvalueCutoff = 1,  
                qvalueCutoff = 1
                )

res <- setReadable(result, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
res <- data.frame(res)
# enrich.go <- enrich.go %>% as.data.frame()
write.table(res, file =  print(paste0(output,".DO.txt")), quote = FALSE, sep = "\t", row.names = FALSE)

res.cutoff <- res %>% as.data.frame() %>% filter(pvalue < 0.05) 
write.table(res.cutoff, file =  print(paste0(output,".DO.cutoff.txt")), quote = FALSE, sep = "\t", row.names = FALSE)

cnetplot(result,categorySize="pvalue",node_label = "all",
         showCategory=nrow(res.cutoff))
