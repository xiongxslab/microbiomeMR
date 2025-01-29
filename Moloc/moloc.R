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
