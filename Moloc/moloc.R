######Moloc for the micobiome, methylation and gene expression
library(data.table)
library(dplyr)
library(moloc)
library(foreach)
library(doParallel)


moloc <- function(i){

  maf_file <- paste0("/dell/wanghc/b_file/EUR/",moloc_pairs[i, 5],".frq")

  maf <- fread(maf_file, sep = "\t", header = FALSE) %>%  dplyr::select(SNP = V2,A1 = V3,A2 = V4, MAF = V5)

  maf <- data.frame(maf)



  gtex <- fread(maf_file, sep = "\t", header = FALSE) %>%  dplyr::select(SNP = V2,A1 = V3,A2 = V4)
  gtex <- data.frame(gtex)

  Phe1_dat <- fread(print(paste0("/dell/wanghc/microbiome_QTL/MR_gut_microbiome/gut_microbiome_genus/",samples,"/",samples,"_", moloc_pairs[i, 5], ".1.txt.gz")), sep = "\t", header = FALSE, col.names = c("Phenotype", "SNP", "BETA", "effe
ct_allele", "other_allele", "PVAL", "SE","N")) %>%
    inner_join(y = maf, by = "SNP") %>% mutate(BETA = ifelse(effect_allele == A1, BETA, -BETA))
  Phe2_file <- paste0("/dell/wanghc/microbiome_QTL/MR_gut_microbiome/Colon-microbiome/mQTL-microbiome/mQTL_out/",moloc_pairs[i, 5],".2.txt")
  Phe2_dat <- fread(Phe2_file, sep = "\t", header = FALSE, col.names = c("Phenotype", "SNP", "MAF", "PVAL","BETA","SE","N","effect_allele", "other_allele"))

  gwas_sum_file <- paste0("/dell/wanghc/microbiome_QTL/MR_microbiome_gut/Micorbiome_eQTL/Micorbiome_Transverse/Transverse/",moloc_pairs[i, 5],".1.txt.gz")

  gwas_dat<- fread(file = gwas_sum_file, sep = "\t", header = T ) %>% select (Phenotype, SNP=rsid,  BETA=beta,  PVAL=pval, SE=se,MAF= maf,N=samplesize,  effect_allele,       other_allele) %>% filter( Phenotype == moloc_pairs[[i, 6]])



  Phe1_moloc_dat <- Phe1_dat %>% filter(Phenotype == moloc_pairs[i, exposure])
  Phe2_moloc_dat <- Phe2_dat %>% filter(Phenotype == moloc_pairs[i, outcome])
  if (nrow(Phe1_moloc_dat) == 0) {
    return(NULL)
  }

  if (nrow(Phe2_moloc_dat) == 0) {
    # You can return NULL or an empty data frame, or simply use return() to skip this iteration
  return(NULL)
  }

  com_SNP <- Reduce(intersect, list(Phe1_moloc_dat$SNP, Phe2_moloc_dat$SNP, gwas_dat$SNP))

  if(length(com_SNP) > 1 && sort(gwas_dat[SNP %in% com_SNP,]$PVAL)[1] < 1e-4){
test <- gwas_dat[SNP %in% com_SNP,]
    duplicate_rows <- test%>%
      group_by(SNP) %>%
      filter(n() > 1)
    test_1 <- Phe2_moloc_dat[SNP %in% com_SNP,]
    duplicate_rows_1 <- test_1%>%
      group_by(SNP) %>%
      filter(n() > 1)


    new_file <- duplicate_rows %>%
      inner_join(y = maf, by = "SNP") %>% filter((effect_allele == A1 & other_allele == A2) | (other_allele == A1 & effect_allele == A2)) %>%
      dplyr::select(Phenotype,      SNP,    BETA,   PVAL,   SE,     MAF.x,  N,      effect_allele,other_allele)  %>% rename(MAF = MAF.x)

    new_file_1 <- duplicate_rows_1 %>%
      inner_join(y = maf, by = "SNP") %>% filter((effect_allele == A1 & other_allele == A2) | (other_allele == A1 & effect_allele == A2)) %>%
      dplyr::select(Phenotype,      SNP,    BETA,   PVAL,   SE,     MAF.x,  N,      effect_allele,other_allele)  %>% rename(MAF = MAF.x)

    if (nrow(new_file) > 0&& nrow(new_file_1) > 0) {

      original_data_cleaned <- anti_join(test, duplicate_rows, by = "SNP")
      original_data_cleaned_1 <- anti_join(test_1, duplicate_rows_1, by = "SNP")
      final_data <- bind_rows(original_data_cleaned, new_file)
      final_data_1<- bind_rows(original_data_cleaned_1, new_file_1)
      listData <- list(Phe1_moloc_dat[SNP %in% com_SNP,], final_data_1, final_data)
      moloc_res <- moloc_test(listData, prior_var = "default", save.SNP.info = FALSE)
      c_res <- c(t(moloc_res$best_snp))
      names(c_res) <- c("pp_a", "best_SNP_a", "pp_b", "best_SNP_b", "pp_ab", "best_SNP_ab", "pp_c", "best_SNP_c", "pp_ac", "best_SNP_ac", "pp_bc", "best_SNP_bc", "pp_abc", "best_SNP_abc")
      moloc_combine_res0 <- data.frame(exposure = moloc_pairs[i, exposure], outcome = moloc_pairs[i, outcome],gene= moloc_pairs[i,gene],chr=moloc_pairs[i,chr], nsnps = moloc_res$nsnps, t(c_res))

      return(moloc_combine_res0)} else if (nrow(new_file) == 0 && nrow(new_file_1) > 0)  {filtered_data <- test %>%
        group_by(SNP) %>%
        filter(PVAL == min(PVAL)) %>%
        ungroup()
      original_data_cleaned_1 <- anti_join(test_1, duplicate_rows_1, by = "SNP")
      final_data_1<- bind_rows(original_data_cleaned_1, new_file_1)




      listData <- list(Phe1_moloc_dat[SNP %in% com_SNP,], final_data_1,filtered_data )
      moloc_res <- moloc_test(listData, prior_var = "default", save.SNP.info = FALSE)
      c_res <- c(t(moloc_res$best_snp))
      names(c_res) <- c("pp_a", "best_SNP_a", "pp_b", "best_SNP_b", "pp_ab", "best_SNP_ab", "pp_c", "best_SNP_c", "pp_ac", "best_SNP_ac", "pp_bc", "best_SNP_bc", "pp_abc", "best_SNP_abc")
      moloc_combine_res0 <- data.frame(exposure = moloc_pairs[i, exposure], outcome = moloc_pairs[i, outcome], gene= moloc_pairs[i,gene],chr=moloc_pairs[i,chr],nsnps = moloc_res$nsnps, t(c_res))
 return(moloc_combine_res0)

      } else if (nrow(new_file) > 0 && nrow(new_file_1) == 0) {filtered_data_1 <- test_1 %>%
        group_by(SNP) %>%
        filter(PVAL == min(PVAL)) %>%
        ungroup()
      original_data_cleaned <- anti_join(test, duplicate_rows, by = "SNP")
      final_data<- bind_rows(original_data_cleaned, new_file)




      listData <- list(Phe1_moloc_dat[SNP %in% com_SNP,], filtered_data_1,final_data )
      moloc_res <- moloc_test(listData, prior_var = "default", save.SNP.info = FALSE)
      c_res <- c(t(moloc_res$best_snp))
      names(c_res) <- c("pp_a", "best_SNP_a", "pp_b", "best_SNP_b", "pp_ab", "best_SNP_ab", "pp_c", "best_SNP_c", "pp_ac", "best_SNP_ac", "pp_bc", "best_SNP_bc", "pp_abc", "best_SNP_abc")
      moloc_combine_res0 <- data.frame(exposure = moloc_pairs[i, exposure], outcome = moloc_pairs[i, outcome], gene= moloc_pairs[i,gene],chr=moloc_pairs[i,chr],nsnps = moloc_res$nsnps, t(c_res))



      } else {filtered_data_1 <- test_1 %>%
        group_by(SNP) %>%
        filter(PVAL == min(PVAL)) %>%
        ungroup()
      filtered_data <- test %>%
        group_by(SNP) %>%
        filter(PVAL == min(PVAL)) %>%
        ungroup()

      listData <- list(Phe1_moloc_dat[SNP %in% com_SNP,], filtered_data_1,filtered_data )
      moloc_res <- moloc_test(listData, prior_var = "default", save.SNP.info = FALSE)
      c_res <- c(t(moloc_res$best_snp))
      names(c_res) <- c("pp_a", "best_SNP_a", "pp_b", "best_SNP_b", "pp_ab", "best_SNP_ab", "pp_c", "best_SNP_c", "pp_ac", "best_SNP_ac", "pp_bc", "best_SNP_bc", "pp_abc", "best_SNP_abc")
      moloc_combine_res0 <- data.frame(exposure = moloc_pairs[i, exposure], outcome = moloc_pairs[i, outcome], gene= moloc_pairs[i,gene],chr=moloc_pairs[i,chr],nsnps = moloc_res$nsnps, t(c_res))

      return(moloc_combine_res0)

      }

     }
  else {return(NULL)
  }
}






args <- commandArgs(T)
samples <- args[1]

chr_file <- paste0("/dell/wanghc/microbiome_QTL/MR_gut_microbiome/Colon-microbiome/mQTL-microbiome/moloc_mbQTL_mQTL_pair/","moloc-",samples,"-Transverse_mQTL.3.pair")

if (!file.exists(chr_file)) {
  return(NULL) # Skip this iteration
}

moloc_pair <- fread(file=chr_file, sep = "\t", header = F)
if(nrow(moloc_pair)==0){
  return(NULL)
}


moloc_pairs <- fread(chr_file, sep = "\t", header = FALSE) %>% dplyr::select(c(1,2,8,15,16,17)) %>% dplyr::rename(exposure = V1, outcome = V2, fdr = V8, PP4 = V15,chr = V16, gene = V17)



out_file <- paste0("/data/slurm/wanghc/Microbiome_eQTL/Moloc_microbiome_mQTL_eQTL/",
                   "/moloc-",samples,"-mQTL-eQTL/","moloc-",samples,"-mQTL-eQTL",".moloc.pp")


cl <- makeCluster(8)
registerDoParallel(cl)
clusterEvalQ(cl, .libPaths("/dell/wanghc/R_1/4.2"))
res <- foreach(i = 1:nrow(moloc_pairs), .combine = rbind, .packages = c("data.table", "dplyr", "moloc"), .errorhandling = "pass") %dopar% {
  moloc_result <- moloc(i)
  if (is.null(moloc_result)) {
    return(NULL)
  }
  moloc_result
}

write.table(res, file = out_file, quote = FALSE, sep = "\t", row.names = FALSE)
stopImplicitCluster()
stopCluster(cl)



######Moloc for gene expression, microbiome and diseases
for (h in samples){
  for (j in 1:22) {
    
    moloc <- function(i){
      Phe1_moloc_dat <- Phe1_dat %>% filter(Phenotype == moloc_pairs[i, exposure])
      Phe2_moloc_dat <- Phe2_dat %>% filter(Phenotype == moloc_pairs[i, outcome])
      if (nrow(Phe1_moloc_dat) == 0) {
       
        return(NULL)
      }
      
      if (nrow(Phe2_moloc_dat) == 0) {
        
        return(NULL)
      }
      
      com_SNP <- Reduce(intersect, list(Phe1_moloc_dat$SNP, Phe2_moloc_dat$SNP, gwas_dat$SNP))
      
      if(length(com_SNP) > 1 && sort(gwas_dat[SNP %in% com_SNP,]$PVAL)[1] < 1e-4){
        listData <- list(Phe1_moloc_dat[SNP %in% com_SNP,], Phe2_moloc_dat[SNP %in% com_SNP,], gwas_dat[SNP %in% com_SNP,])
        moloc_res <- moloc_test(listData, prior_var = "default", save.SNP.info = FALSE)
        c_res <- c(t(moloc_res$best_snp))
        names(c_res) <- c("pp_a", "best_SNP_a", "pp_b", "best_SNP_b", "pp_ab", "best_SNP_ab", "pp_c", "best_SNP_c", "pp_ac", "best_SNP_ac", "pp_bc", "best_SNP_bc", "pp_abc", "best_SNP_abc")
        moloc_combine_res0 <- data.frame(exposure = moloc_pairs[i, exposure], outcome = moloc_pairs[i, outcome], nsnps = moloc_res$nsnps, t(c_res))
       
        return(moloc_combine_res0)
      }
      else {return(NULL)  # Stop further execution for this iteration
      }
    }
    
    
    chr_file <- paste0("/data/slurm/wanghc/microbiome_QTL/moloc_gut_microbiome/moloc_Ileum_microbiome/moloc_Ileum-eQTL_microbiome/moloc_Ileum-eQTL_microbiome_pair/moloc_Ileum-eQTL-",h,"-pair",
                       "/moloc_Ileum-eQTL-",h,".pair")
   
    Phe1_file <- paste0("/data/slurm/wanghc/microbiome_QTL/MR_microbiome_gut/Micorbiome_eQTL/Micorbiome_Ileum/Ileum/chr", j, ".3.txt")
    Phe2_file <- paste0("/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/Colon-microbiome/","family.",h,"/","family.",h,"_chr", j, ".1.txt")
    maf_file <- paste0("/data/slurm/licy/QTL_integration/1.Database/Ref_panel/1kg_v3/EUR.frq")
    gwas_sum_file <- paste0("/data/slurm/wanghc/microbiome_QTL/moloc_gut_microbiome/moloc_Ileum_microbiome/moloc_Ileum-eQTL_microbiome/phecode-555-IBD/","phecode-555-both_sexes.16.tsv")
    
    out_file <- paste0("/data/slurm/wanghc/microbiome_QTL/moloc_gut_microbiome/moloc_Ileum_microbiome/moloc_Ileum-eQTL_microbiome/", "phecode-555-IBD/","moloc_Ileum-eQTL-",h,"-IBD/",
                       "moloc_Ileum-eQTL-",h,"-IBD",".chr",j,".moloc.pp")
    
    moloc_pairs <- fread(chr_file, sep = "\t", header = FALSE) %>% dplyr::select(c(1,2,8,15)) %>% dplyr::rename(exposure = V1, outcome = V2, fdr = V8, PP4 = V15)
    maf <- fread(maf_file, sep = "\t", header = TRUE) %>% dplyr::select(SNP, A1, A2, MAF)
    maf <- data.frame(maf)
   
    Phe1_dat <- fread(Phe1_file, sep = "\t", header = FALSE, col.names = c("Phenotype",	"SNP",	"BETA",	"PVAL",	"SE",	"MAF",	"N",	"effect_allele",	"other_allele")) 
    
    Phe2_dat <- fread(Phe2_file, sep = "\t", header = FALSE, col.names = c("Phenotype", "SNP", "BETA", "effect_allele", "other_allele", "PVAL", "SE","N")) %>%
      inner_join(y = maf, by = "SNP") %>% mutate(BETA = ifelse(effect_allele == A1, BETA, -BETA))
    
    gtex <-  fread(maf_file, sep = "\t", header = TRUE) %>% dplyr::select(SNP, A1, A2)
    gtex <- data.frame(gtex)
   
    gwas_dat <- fread(gwas_sum_file, sep = "\t", header = TRUE, col.names = c("SNP", "A1_gwas", "A2_gwas", "BETA", "SE", "MAF", "PVAL", "N" ,"Ncases"))  %>%
      inner_join(y = gtex, by = "SNP") %>%filter((A1_gwas == A1 & A2_gwas == A2) | (A2_gwas == A1 & A1_gwas == A2))
    
    gwas_dat <- gwas_dat %>% mutate(MAF = ifelse(A1_gwas == A1, MAF, 1 - MAF)) %>% select(SNP, A1, A2, MAF, PVAL, N ,BETA, SE, Ncases)
    
    Phe1_dat <- Phe1_dat[complete.cases(Phe1_dat), ]
    Phe2_dat <- Phe2_dat[complete.cases(Phe2_dat), ]
    gwas_dat <- gwas_dat[complete.cases(gwas_dat), ]
    
    
    cl <- makeCluster(8)
    registerDoParallel(cl)
    clusterEvalQ(cl, .libPaths("/data/slurm/wanghc/R"))
    res <- foreach(i = 1:nrow(moloc_pairs), .combine = rbind, .packages = c("data.table", "dplyr", "moloc"), .errorhandling = "pass") %dopar% {
      moloc_result <- moloc(i)
      if (is.null(moloc_result)) {
        return(NULL)
      }
      moloc_result
    }
    
    write.table(res, file = out_file, quote = FALSE, sep = "\t", row.names = FALSE)
   
    stopImplicitCluster()
    stopCluster(cl)}}
