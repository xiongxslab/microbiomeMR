######MR analysis 
.libPaths("/data/slurm/wanghc/R_1/4.2")
library(data.table)
library(dplyr)
library(TwoSampleMR)
library(ieugwasr)
library(plinkbinr)
library(foreach)
library(doParallel)

#####generate pairs (take transverse colon as example)

args <- commandArgs(T)
samples <- args[1]

Chr <- args[2]
 Phe1.dat <- paste0("/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/Colon-microbiome/eQTL-microbiome/Colon_Transerve.eQTL__expo/",Chr, ".1.txt")
  Phe2.dat <- paste0("/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/gut_microbiome_genus/",samples,"/",samples,"_", Chr, ".1.txt.gz")
    Phe1_dat_test <- fread(file =  Phe1.dat, sep = "\t", header = FALSE, 
                           col.names = c("Phenotype", "SNP", "beta", "pval","se","eaf","effect_allele", "other_allele","samplesize"))
    # Read the second file into a data frame
    
    
    Phe2_dat_test<- fread(file = Phe2.dat, sep = "\t", header = FALSE, 
                          col.names = c("Phenotype", "SNP", "beta", "effect_allele", "other_allele", "pval", "se", "samplesize")) 
    
    
    merged_data <- merge(Phe1_dat_test, Phe2_dat_test, by.x = "SNP", by.y = "SNP", all = FALSE)
    
  
    final_data <- merged_data[, c("Phenotype.x", "Phenotype.y")]
    
    data <-  final_data %>% group_by(Phenotype.x, Phenotype.y) %>% summarise_all(list(~paste(., collapse = ","))) %>%
      ungroup() 
    
  
    write.table(data, file=print(paste0("/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/Transverse-microbiome_genus_stringent/","Transverse_eQTL-",samples,"/Transverse
_eQTL-",samples,".",Chr,".pair")), quote = 
                  FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


######## MR analysis (take sigmoid tissue as an example)
args <- commandArgs(T)
samples <- args[1]

Chr <- args[2]

    Phe1.dat <- paste0("/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/Colon-microbiome/eQTL-microbiome/Colon_Sigmoid.eQTL__expo/",Chr, ".2.txt")
    Phe2.dat <- paste0("/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/gut_microbiome_genus/",samples,"/",samples,"_", Chr, ".1.txt.gz")
    mr.pairs <- paste0("/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/Sigmoid-microbiome_genus_stringent/Sigmoid_eQTL-",samples,"/Sigmoid_eQTL-",samples,".",Chr,".pair")
 maf <- fread(file=print(paste0("/data/slurm/wanghc/b_file/EUR/",Chr,".frq")), sep = "\t", header = FALSE) %>%  dplyr::select(SNP = V2, MAF = V5)%>% dplyr::rename(eaf = 
MAF)
    
Phe1_dat_test <- fread(file = Phe1.dat, sep = "\t", header = FALSE,
                           col.names = c("Phenotype", "SNP", "beta", "pval","se","eaf","effect_allele", "other_allele"))%>% mutate(samplesize = 318)
Phe2_dat_test<- fread(file = Phe2.dat, sep = "\t", header = FALSE,
                          col.names = c("Phenotype", "SNP", "beta", "effect_allele", "other_allele", "pval", "se", "samplesize"))
mr_pairs <- fread(file= mr.pairs,
                      sep = "\t", header = FALSE, col.names = c("eQTL", "bac"))
    mr_pairs <- data.frame(mr_pairs)
b_file <- paste0("/data/slurm/wanghc/b_file/",Chr)
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
 }else if(length(mr_dat$SNP) > 1){
          mr_res <- mr(mr_dat, method_list = "mr_ivw")
        }
 mr_res <- mr_res%>% dplyr::select(-id.exposure, -id.outcome) %>% dplyr::select(exposure, everything())
        return(mr_res)
      }
    }
    cl <- makeCluster(8)
    registerDoParallel(cl)
    clusterEvalQ(cl, .libPaths("/data/slurm/wanghc/R_1/4.2/"))
    res<- foreach(i = 1:nrow(mr_pairs), .combine = rbind, .packages = c("dplyr", "TwoSampleMR", "ieugwasr", "plinkbinr","data.table"), .errorhandling = "pass") %dopar% {
      perform_mr(i)
    }
   write.table(res, file = print(paste0("/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/Sigmoid-microbiome_genus_stringent/","Sigmoid_eQTL-",samples,"/Sigmoid_eQTL-",samples,".",Chr,".MR")), quote = FALSE, sep = "\t", row.names = FALSE)
    stopImplicitCluster()
    stopCluster(cl)


###### sensitivity test 
####Pleiotropy test
args <- commandArgs(T)
samples <- args[1]

Chr <- args[2]
Phe1.dat <- paste0("/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/Colon-microbiome/eQTL-microbiome/Colon_Transerve.eQTL__expo/",Chr, ".1.txt")
Phe2.dat <- paste0("/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/gut_microbiome_genus/",samples,"/",samples,"_", Chr, ".1.txt.gz")
mr.pairs <- paste0("/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/Transverse-microbiome_genus_stringent/Transverse_eQTL-",samples,"/Transverse_eQTL-",samples,".",Chr,".pair")

maf <- fread(file=print(paste0("/data/slurm/wanghc/b_file/EUR/",Chr,".frq")), sep = "\t", header = FALSE) %>%  dplyr::select(SNP = V2, MAF = V5)%>% dplyr::rename(eaf = MAF)maf <- data.frame(maf)
Phe1_dat_test <- fread(file = Phe1.dat, sep = "\t", header = FALSE,
                       col.names = c("Phenotype", "SNP", "beta", "pval","se","eaf","effect_allele", "other_allele","samplesize"))
Phe2_dat_test<- fread(file = Phe2.dat, sep = "\t", header = FALSE,
                      col.names = c("Phenotype", "SNP", "beta", "effect_allele", "other_allele", "pval", "se", "samplesize"))
mr_pairs <- fread(file= mr.pairs,
                  sep = "\t", header = FALSE, col.names = c("eQTL", "bac"))
b_file <- paste0("/data/slurm/wanghc/b_file/",Chr)
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
                                   effect_allele_col = "effect_allele", other_allele_col = "other_allele",eaf_col = "eaf", pval_col = "pval", samplesize_col = "samplesize")  Phe2_out_dat_test <- format_data(Phe2_out_dat_test, type = "outcome", phenotype_col = "Phenotype", snp_col = "SNP", beta_col = "beta", se_col = "se",
                                   effect_allele_col = "effect_allele", other_allele_col = "other_allele",pval_col = "pval", samplesize_col = "samplesize")
  mr_dat <- harmonise_data(exposure_dat = Phe1_exp_dat_test, outcome_dat = Phe2_out_dat_test) %>% filter(mr_keep == TRUE)
  if(length(mr_dat$SNP) >= 3){
    pleio <- mr_pleiotropy_test(mr_dat)
  }else  {return(NULL)}
}
cl <- makeCluster(8)
registerDoParallel(cl)
clusterEvalQ(cl, .libPaths("/data/slurm/wanghc/R_1/4.2/"))
res<- foreach(i = 1:nrow(mr_pairs), .combine = rbind, .packages = c("dplyr", "TwoSampleMR", "ieugwasr", "plinkbinr","data.table"), .errorhandling = "pass") %dopar% {
  perform_mr(i)
}
write.table(res, file = print(paste0("/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/Transverse-microbiome_genus_stringent/","Transverse_eQTL-",samples,"/Transverse_eQTL-",samples,".",Chr,".pleio")), quote = FALSE, sep = "\t", row.names = FALSE)
stopImplicitCluster()
stopCluster(cl)
#####leave one out test
args <- commandArgs(T)
samples <- args[1]
Chr <- args[2]
Phe1.dat <- paste0("/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/Colon-microbiome/eQTL-microbiome/Colon_Transerve.eQTL__expo/",Chr, ".1.txt")
Phe2.dat <- paste0("/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/gut_microbiome_genus/",samples,"/",samples,"_", Chr, ".1.txt.gz")
mr.pairs <- paste0("/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/Transverse-microbiome_genus_stringent/Transverse_eQTL-",samples,"/Transverse_eQTL-",samples,".",Chr,".pair")
maf <- fread(file=print(paste0("/data/slurm/wanghc/b_file/EUR/",Chr,".frq")), sep = "\t", header = FALSE) %>%  dplyr::select(SNP = V2, MAF = V5)%>% dplyr::rename(eaf = MAF)

Phe1_dat_test <- fread(file = Phe1.dat, sep = "\t", header = FALSE,
                       col.names = c("Phenotype", "SNP", "beta", "pval","se","eaf","effect_allele", "other_allele","samplesize"))

Phe2_dat_test<- fread(file = Phe2.dat, sep = "\t", header = FALSE,
                      col.names = c("Phenotype", "SNP", "beta", "effect_allele", "other_allele", "pval", "se", "samplesize"))
mr_pairs <- fread(file= mr.pairs,
                  sep = "\t", header = FALSE, col.names = c("eQTL", "bac"))
b_file <- paste0("/data/slurm/wanghc/b_file/",Chr)
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
                                   effect_allele_col = "effect_allele", other_allele_col = "other_allele",eaf_col = "eaf", pval_col = "pval", samplesize_col = "samplesize")  Phe2_out_dat_test <- format_data(Phe2_out_dat_test, type = "outcome", phenotype_col = "Phenotype", snp_col = "SNP", beta_col = "beta", se_col = "se",
                                   effect_allele_col = "effect_allele", other_allele_col = "other_allele",pval_col = "pval", samplesize_col = "samplesize")
  mr_dat <- harmonise_data(exposure_dat = Phe1_exp_dat_test, outcome_dat = Phe2_out_dat_test) %>% filter(mr_keep == TRUE)
  if(length(mr_dat$SNP) >= 4){
    loo <- mr_leaveoneout(mr_dat)
     }else {return(NULL)}
}
cl <- makeCluster(8)
registerDoParallel(cl)
clusterEvalQ(cl, .libPaths("/data/slurm/wanghc/R_1/4.2/"))
res<- foreach(i = 1:nrow(mr_pairs), .combine = rbind, .packages = c("dplyr", "TwoSampleMR", "ieugwasr", "plinkbinr","data.table"), .errorhandling = "pass") %dopar% {
  perform_mr(i)
}
 write.table(res, file = print(paste0("//data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/Transverse-microbiome_genus_stringent/","Transverse_eQTL-",samples,"/Transverse_eQTL-",samples,".",Chr,".rmloo")), quote = FALSE, sep = "\t", row.names = FALSE,col.names = T)
stopImplicitCluster()
stopCluster(cl)
#####Pass the sensitivity test and HEIDI test and study-wise FDR

MR_1 <-  fread(file = "/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/study_wise_BH_MR/Transverse_eQTL.all.MR" , sep = "\t", header = T)%>% mutate(Tissue="Transverse")
MR_2 <-  fread(file = "/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/study_wise_BH_MR/Sigmoid_eQTL.all.MR" , sep = "\t", header = T)%>% mutate(Tissue="Sigmoid")

MR_3 <-  fread(file = "/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/study_wise_BH_MR/Ileum_eQTL.all.MR" , sep = "\t", header = T)%>% mutate(Tissue="Ileum")

MR <- bind_rows(MR_1, MR_2, MR_3)

Pleio_1 <- fread(file = "/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/study_wise_BH_MR/Transverse_eQTL.all.pleio" , sep = "\t", header = T) %>% mutate(Tissue="Transverse")

Pleio_2 <- fread(file = "/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/study_wise_BH_MR/Sigmoid_eQTL.all.pleio" , sep = "\t", header = T) %>% mutate(Tissue="Sigmoid")
Pleio_3 <- fread(file = "/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/study_wise_BH_MR/Ileum_eQTL.all.pleio" , sep = "\t", header = T) %>% mutate(Tissue="Ileum")

Pleio <- bind_rows(Pleio_1, Pleio_2, Pleio_3)



Pleio_filtered <- Pleio %>%
  filter(pval < 0.05) %>%
  select(exposure, outcome) %>%
  distinct()

MR_clean <- MR %>%
  anti_join(Pleio_filtered, by = c("exposure", "outcome"))


loo_1 <- fread(file = "/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/study_wise_BH_MR/Transverse_eQTL.all.rmloo" , sep = "\t", header = T) %>% mutate(Tissue="Transverse")

loo_2 <- fread(file = "/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/study_wise_BH_MR/Sigmoid_eQTL.all.rmloo" , sep = "\t", header = T) %>% mutate(Tissue="Sigmoid")
loo_3 <- fread(file = "/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/study_wise_BH_MR/Ileum_eQTL.all.rmloo" , sep = "\t", header = T) %>% mutate(Tissue="Ileum")

res <- loo %>%
  group_by(exposure, outcome) %>%
  filter(
    any(SNP == "All" & p < 0.05) &

      sum(SNP != "All" & p > 0.05) == 1
  ) %>%
  ungroup()%>% filter(SNP == "All") %>%
  select(exposure, outcome) %>%
  distinct()

MR_clean <- MR_clean %>%
  anti_join(res, by = c("exposure", "outcome"))

MR_clean$adjusted_pval <- p.adjust(MR_clean$pval, method="fdr")

significant_results <- MR_clean[MR_clean$adjusted_pval < 0.05, ]

HEIDI <- fread(file = "/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/study_wise_BH_MR/Tissue.all.smr" , sep = "\t", header = T)  %>% rename(exposure = Gene)
  
HEIDI_filtered <- HEIDI %>%
      filter( p_HEIDI< 0.01) %>%
      select(exposure, outcome) %>%
      distinct()
MR_final <- MR_clean %>% anti_join(HEIDI_filtered, by = c("exposure", "outcome"))

MR_transverse <- MR_final %>% filter(Tissue == "Transverse")
write.table(MR_transverse, file = "/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/study_wise_BH_MR/study_wise_BH/Transverse_eQTL.filtered.smr.gene.loo.MR", quote = FALSE, sep = "\t", row.names = FALSE)
MR_sigmoid <- MR_final %>% filter(Tissue == "Sigmoid")
write.table(MR_sigmoid, file = "/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/study_wise_BH_MR/study_wise_BH/Sigmoid_eQTL.filtered.smr.gene.loo.MR", quote = FALSE, sep = "\t", row.names = FALSE)
MR_ileum <- MR_final %>% filter(Tissue == "Ileum")
write.table(MR_ileum, file = "/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/study_wise_BH_MR/study_wise_BH/Ileum_eQTL.filtered.smr.gene.loo.MR", quote = FALSE, sep = "\t", row.names = FALSE)


