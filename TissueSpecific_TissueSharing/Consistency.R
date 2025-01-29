########consistency of effect sizes
 
         MR_data_Transverse <- fread(file = paste0("/data/slurm/wanghc/microbiome_QTL/SMR/Transverse2microbiome/Transverse-microbiome-smr/Transverse-microbiome.1.smr"), sep = "\t", header = FALSE, 
                                     col.names = c("exposure",       "outcome", "method",  "nsnp",   "b_Transverse" ,     "se_Transverse",     "pval_Transverse",   "adjusted_pval_Transverse",   "gene_name_Transverse","CHR_Transverse","START_Transverse","END_Transverse"))
         
         MR_data_Ileum <- fread(file = paste0("/data/slurm/wanghc/microbiome_QTL/SMR/Ileum2microbiome/Ileum-microbiome-smr/Ileum-microbiome.1.smr"), sep = "\t", header = FALSE, 
                                col.names = c("exposure",       "outcome", "method",  "nsnp",   "b_Ileum" ,     "se_Ileum",     "pval_Ileum",   "adjusted_pval_Ileum",   "gene_name_Ileum","CHR_Ileum","START_Ileum","END_Ileum"))
         
         
        
         # Merge the two datasets by "exposure" and "outcome"
         common_rows <- merge(MR_data_Transverse, MR_data_Ileum, by = c("exposure", "outcome"))
         
         # Check sign consistency of "b" values
         common_rows[, consistent_sign := sign(b_Transverse) == sign(b_Ileum)]
         
         # Calculate the percentage of rows with consistent signs
         total_rows <- nrow(common_rows)
         consistent_rows <- sum(common_rows$consistent_sign, na.rm = TRUE)
         percentage_consistent <- (consistent_rows / total_rows) * 100
         
         
         
         
         
         #################################### comparison of effect sizes between transverse and sigmoid
         
           MR_data_Transverse <- fread(file = paste0("/data/slurm/wanghc/microbiome_QTL/SMR/Transverse2microbiome/Transverse-microbiome-smr/Transverse-microbiome.1.smr"), sep = "\t", header = FALSE, 
                                     col.names = c("exposure",       "outcome", "method",  "nsnp",   "b" ,     "se",     "pval",   " adjusted_pval",   "gene_name","CHR","START","END"))
         
         MR_data_Sigmoid <- fread(file = paste0("/data/slurm/wanghc/microbiome_QTL/SMR/Sigmoid2microbiome/Sigmoid-microbiome-smr/Sigmoid-microbiome.1.smr"), sep = "\t", header = FALSE, 
                                                   col.names = c("exposure", "outcome", "method", "nsnp", "b", "se", "pval", "adjusted_pval", "gene_name","CHR","START","END"))
                                  
         
        
         common_rows <- merge(MR_data_Transverse, MR_data_Sigmoid, by = c("exposure", "outcome"))
         
         b_transverse <- common_rows$`b.x`  # b values from Transverse
         b_Sigmoid <- common_rows$`b.y`# b values from Sigmoid 
         
         df <- data.frame(b_transverse, b_Sigmoid)
        
         regression_model <- lm(b_Sigmoid ~ b_transverse, data = df)
         
         model_summary <- summary(regression_model)
         


         #################################### comparison of effect sizes between ileum and sigmoid

 MR_data_Ileum <- fread(file = paste0("/data/slurm/wanghc/microbiome_QTL/SMR/Ileum2microbiome/Ileum-microbiome-smr/Ileum-microbiome.1.smr"), sep = "\t", header = FALSE, 
                                     col.names = c("exposure",       "outcome", "method",  "nsnp",   "b" ,     "se",     "pval",   " adjusted_pval",   "gene_name","CHR","START","END"))
         
         MR_data_Sigmoid <- fread(file = paste0("/data/slurm/wanghc/microbiome_QTL/SMR/Sigmoid2microbiome/Sigmoid-microbiome-smr/Sigmoid-microbiome.1.smr"), sep = "\t", header = FALSE, 
                                                   col.names = c("exposure", "outcome", "method", "nsnp", "b", "se", "pval", "adjusted_pval", "gene_name","CHR","START","END"))
                                  
         
        
         common_rows <- merge(MR_data_Ileum, MR_data_Sigmoid, by = c("exposure", "outcome"))
         
         b_Ileum <- common_rows$`b.x`  # b values from Transverse
         b_Sigmoid <- common_rows$`b.y`# b values from Sigmoid 
         
         df <- data.frame(b_transverse, b_Sigmoid)
        
         regression_model <- lm(b_Sigmoid ~ b_Ileum, data = df)
         
         model_summary <- summary(regression_model)
         
         
