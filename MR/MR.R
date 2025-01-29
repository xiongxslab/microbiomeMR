
library(data.table)
library(dplyr)
library(TwoSampleMR)
library(ieugwasr)
library(plinkbinr)
library(foreach)
library(doParallel)





#####generate the gene-to-microbiome pairs in the sigmoid colon  

tissues <- c("Transverse", "Sigmoid", "ileum")
samples <-  fread(file = "/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/gut_microbiome_genus/genus.txt",sep = "\t", header = FALSE)


 for (a in tissues) {
for (h in 1:nrow(samples)) {
  
  for (i in 1:22) {
    bacteria <- samples$V1[h]
    
    # Read the first file into a data frame
    
    Phe1_dat <- paste0("/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/Colon-microbiome/eQTL-microbiome/Colon_",a,".eQTL__expo_new/chr", i, ".2.txt")
    Phe2_dat <- paste0("/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/gut_microbiome_genus/",bacteria,"/",bacteria,"_chr", i, ".1.txt")
    Phe1_dat_test <- fread(file =  Phe1_dat, sep = "\t", header = FALSE, 
                           col.names = c("Phenotype", "SNP", "beta", "pval","se","maf","effect_allele", "other_allele")) %>% mutate(samplesize = 318)
    # Read the second file into a data frame
    
    
    Phe2_dat_test<- fread(file = Phe2_dat, sep = "\t", header = FALSE, 
                          col.names = c("Phenotype", "SNP", "beta", "effect_allele", "other_allele", "pval", "se", "samplesize")) 
    
    
    merged_data <- merge(Phe1_dat_test, Phe2_dat_test, by.x = "SNP", by.y = "SNP", all = FALSE)
    
  
    final_data <- merged_data[, c("Phenotype.x", "Phenotype.y")]
    
    data <-  final_data %>% group_by(Phenotype.x, Phenotype.y) %>% summarise_all(list(~paste(., collapse = ","))) %>%
      ungroup() 
    
    
    
    # Write the combined phenotype vector to a new file
    write.table(data, file=print(paste0("/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/Sigmoid-microbiome_genus/",a,"_eQTL-",bacteria,"/",a,"_eQTL-",bacteria,".chr",i,".pair")), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)}
}}

# MR analysis in three tissues

maf <- fread(file="/data/slurm/licy/QTL_integration/1.Database/Ref_panel/1kg_v3/EUR.frq", sep = "\t", header = TRUE) %>% dplyr::select(SNP, MAF) %>% dplyr::rename(eaf = MAF)
maf <- data.frame(maf)
for (a in tissues) {
for (h in 1:nrow(samples)) {
  
  
  for (j in 1:22) {
  
  bacteria <- samples$V1[h]
  Phe1.dat <- paste0("/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/Colon-microbiome/eQTL-microbiome/Colon_",a,".eQTL__expo_new/chr", j, ".2.txt")
  Phe2.dat <- paste0("/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/gut_microbiome_genus/",bacteria,"/",bacteria,"_chr", j, ".1.txt")
  mr.pairs <- paste0("/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/Sigmoid-microbiome_genus/",a,"_eQTL-",bacteria,"/",a,"_eQTL-",bacteria,".chr",j,".pair")
  

  Phe1_dat_test <- fread(file = Phe1.dat, sep = "\t", header = FALSE, 
                         col.names = c("Phenotype", "SNP", "beta", "pval","se","eaf","effect_allele", "other_allele"))%>%
  mutate(samplesize = 318)  
  Phe1_dat_test <- data.frame(Phe1_dat_test)
  
  Phe2_dat_test<- fread(file = Phe2.dat, sep = "\t", header = FALSE, 
                        col.names = c("Phenotype", "SNP", "beta", "effect_allele", "other_allele", "pval", "se", "samplesize"))
 
 Phe2_dat_test <- data.frame(Phe2_dat_test)
  
  mr_pairs <- fread(file= mr.pairs, 
                    sep = "\t", header = FALSE, col.names = c("eQTL", "bac"))
  mr_pairs <- data.frame(mr_pairs)
  
  b_file <- paste0("/data/slurm/wanghc/b_file/chr",j)
  perform_mr <- function(i){
    Phe2_out_dat_test<- Phe2_dat_test %>% filter(Phenotype == mr_pairs[i,'bac'])
    Phe1_exp_dat_test <- Phe1_dat_test %>% filter(Phenotype == mr_pairs[i,'eQTL']) %>% filter(SNP %in% Phe2_out_dat_test$SNP)
    
    check <- 1
    if(length(Phe1_exp_dat_test$SNP) > 1){
      check <<- try(
        Phe1_exp_dat_test<- ld_clump_local(Phe1_exp_dat_test %>% dplyr::rename(rsid = SNP), clump_kb = 100, clump_r2 = 0.01, clump_p = 1,
                                           bfile = b_file, plink_bin = get_plink_exe()) %>% dplyr::rename(SNP = rsid),
        silent = FALSE
      )
    }
    if("try-error" %in% class(check)){
      next
    }
    

    Phe2_out_dat_test <- Phe2_out_dat_test %>% filter(SNP %in% Phe1_exp_dat_test$SNP)
    Phe1_exp_dat_test <- format_data(Phe1_exp_dat_test, type = "exposure", phenotype_col = "Phenotype", snp_col = "SNP", beta_col = "beta", se_col = "se",
                                     effect_allele_col = "effect_allele", other_allele_col = "other_allele",eaf_col = "eaf", pval_col = "pval", samplesize_col = "samplesize")
    Phe2_out_dat_test <- format_data(Phe2_out_dat_test, type = "outcome", phenotype_col = "Phenotype", snp_col = "SNP", beta_col = "beta", se_col = "se",
                                     effect_allele_col = "effect_allele", other_allele_col = "other_allele",pval_col = "pval", samplesize_col = "samplesize")
    mr_dat <- harmonise_data(exposure_dat = Phe1_exp_dat_test, outcome_dat = Phe2_out_dat_test) %>% filter(mr_keep == TRUE)
    if(length(mr_dat$SNP) >= 1){
      if(length(mr_dat$SNP) == 1){
        mr_res <- mr(mr_dat, method_list = "mr_wald_ratio")
      }else if(length(mr_dat$SNP) > 1){
        mr_res <- mr(mr_dat, method_list = "mr_ivw")
      }
      mr_dat_file <- paste0( "/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/"a,"-microbiome_genus/",a,"_eQTL-",bacteria,"/mr.dat/",mr_pairs[i,'eQTL'], "-", mr_pairs[i,'bac'], ".mr.dat")
      write.table(mr_dat, file = mr_dat_file, quote = FALSE, sep = "\t", row.names = FALSE)
      mr_res <- mr_res%>% dplyr::select(-id.exposure, -id.outcome) %>% dplyr::select(exposure, everything())
      return(mr_res)
    }
  }
  cl <- makeCluster(8)
  registerDoParallel(cl)
  clusterEvalQ(cl, .libPaths("/data/slurm/wanghc/R"))
  res <- foreach(i = 1:nrow(mr_pairs), .combine = rbind, .packages = c("dplyr", "TwoSampleMR", "ieugwasr", "plinkbinr"), .errorhandling = "pass") %dopar% {
    perform_mr(i)
  }
  write.table(res, file = print(paste0("/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/Sigmoid-microbiome_genus/",a,"_eQTL-",bacteria,"/",a,"_eQTL-",bacteria,".chr",j,".MR")), quote = FALSE, sep = "\t", row.names = FALSE)
  stopImplicitCluster()
  stopCluster(cl)
    }}}


    ####### FDR 

    tissues <- c("Transverse", "Sigmoid", "ileum")
    for (a in tissues) {
    for (h in samples) {
  


MR.result <- paste0("/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/Colon-microbiome/eQTL-microbiome/",a,"_eQTL-Microbiome/", a,"_eQTL-", h,"/final.test")



MR_Result <- fread(file=MR.result, sep = "\t", header = FALSE, col.names = c('exposure',        'outcome',      'method',       'nsnp', 'b',    'se',   'pval'))

MR_Result$adjusted_pval <- p.adjust(MR_Result$pval, method="fdr")

significant_results <- MR_Result[MR_Result$adjusted_pval < 0.05, ]

write.table(significant_results, file = print(paste0("/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/Colon-microbiome/eQTL-microbiome/",a,"_eQTL-Microbiome/", a,"_eQTL-", h,a,
                                                   "/_eQTL-",h,".result.MR")), quote = FALSE, sep = "\t", row.names = FALSE)}}

