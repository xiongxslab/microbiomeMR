samples <-  c("Bifidobacteriaceae.id.433","Enterobacteriaceae.id.3469", "Oxalobacteraceae.id.2966", "Rikenellaceae.id.967","Verrucomicrobiaceae.id.4036",
              "Christensenellaceae.id.1866", "Erysipelotrichaceae.id.2149", "Pasteurellaceae.id.3689","Ruminococcaceae.id.2050","Victivallaceae.id.2255",
              "Clostridiaceae1.id.1869","Peptococcaceae.id.2024","Streptococcaceae.id.1850","ClostridialesvadinBB60group.id.11286","Peptostreptococcaceae.id.2042",
              "Coriobacteriaceae.id.811","Lachnospiraceae.id.1987",
              "Porphyromonadaceae.id.943","Defluviitaleaceae.id.1924", 
              "Lactobacillaceae.id.1836","Prevotellaceae.id.960","BacteroidalesS24.7group.id.11173",
              "Desulfovibrionaceae.id.3169","Methanobacteriaceae.id.121","Rhodospirillaceae.id.2717","Veillonellaceae.id.2172","Acidaminococcaceae.id.2166","Actinomycetaceae.id.421",
              "Alcaligenaceae.id.2875","Bacteroidaceae.id.917""Desulfovibrionaceae.id.3169", "Methanobacteriaceae.id.121", "Rhodospirillaceae.id.2717" ,"Veillonellaceae.id.2172")

tissues <- c("Transverse", "Sigmoid", "ileum")
for (a in tissues) {
for (h in samples) {
  for (j in 1:22) {
    
    Phe1_dat <- paste0("/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/Colon-microbiome/","family.",h,"/","family.",h,"_chr", j, ".1.txt")
    
    Phe2_dat <- paste0("/data/slurm/wanghc/microbiome_QTL/MR_microbiome_gut/Micorbiome_eQTL/Micorbiome_",a,"/",a,"/chr", j, ".3.txt")
    
   
    chr_file <- paste0("/data/slurm/wanghc/microbiome_QTL/MR_microbiome_gut/Micorbiome_eQTL/Micorbiome_",a,"/",h,"-",a,"/",h,"-",a,".chr",j,".MR")
    
    maf <- fread(file="/data/slurm/licy/QTL_integration/1.Database/Ref_panel/1kg_v3/EUR.frq", sep = "\t", header = TRUE) %>% select(SNP, MAF)
    maf <- data.frame(maf)
    
    
    
    Phe1_dat_test <- fread(file=Phe1_dat, sep = "\t", header = FALSE, col.names = c("Phenotype", "SNP", "beta", "effect_allele", "other_allele", "pval", "se","samplesize")) %>% 
      inner_join(y = maf, by = "SNP") %>%
      mutate(varbeta = se*se)
    Phe1_dat_test <- data.frame(Phe1_dat_test)
    
    
   
    Phe2_dat_test <- fread(file=Phe2_dat, sep = "\t", header = FALSE, col.names = c("Phenotype",	"SNP",	"beta",	"pval",	"se",	"MAF",	"samplesize",	"effect_allele",	"other_allele")) %>%
     
      mutate(varbeta = se*se)
    Phe2_dat_test <- data.frame(Phe2_dat_test)
   
    is.readable.file <- function(file) {
      file.exists(file) && file.info(file)$size > 1 && length(readLines(file, n = 1)) > 0
    }
    
    
    
    if (is.readable.file(chr_file)) {
      
      coloc_pairs <- fread(file=chr_file, sep = "\t", header = TRUE) %>% select(exposure, outcome)
      coloc_pairs <- data.frame(coloc_pairs)
    } else {
      NULL
    }
    
    
    coloc_abf <- function(i){
      Phe1_coloc_dat <- Phe1_dat_test %>% filter(Phenotype == coloc_pairs[i,'exposure'])
      Phe2_coloc_dat <- Phe2_dat_test %>% filter(Phenotype == coloc_pairs[i,'outcome'])
      coloc_dat <- inner_join(Phe1_coloc_dat, Phe2_coloc_dat, by = "SNP", suffix = c("_bac", "_eQTL"))
      
      
      if (nrow(coloc_dat) == 0) {
       
        return(NULL)
      }
      
      coloc_dat <- coloc_dat[
        (coloc_dat$effect_allele_bac == coloc_dat$effect_allele_eQTL | 
           coloc_dat$effect_allele_bac == coloc_dat$other_allele_eQTL) & 
          (coloc_dat$other_allele_bac == coloc_dat$effect_allele_eQTL | 
             coloc_dat$other_allele_bac == coloc_dat$other_allele_eQTL),]
      
      
      if (nrow(coloc_dat) == 0) {
       
        return(NULL)
      }
      
      
      
      
      coloc_res <- coloc.abf(dataset1 = list(snp = coloc_dat$SNP, beta = coloc_dat$beta_bac, varbeta = coloc_dat$varbeta_bac,
                                             MAF = coloc_dat$MAF_bac, N =coloc_dat$samplesize_bac , type = "quant"),
                             dataset2 = list(snp = coloc_dat$SNP, beta = coloc_dat$beta_eQTL, varbeta = coloc_dat$varbeta_eQTL,
                                             MAF = coloc_dat$MAF_eQTL, N = coloc_dat$samplesize_eQTL, type = "quant"))
      
      
      coloc_results_file <- paste0("/data/slurm/wanghc/microbiome_QTL/Coloc_microbiome_gut/Coloc_Microbiome-",a,"_eQTL/","Coloc_",h,"-",a,"_eQTL", "/", coloc_pairs[i,'exposure'], "-", coloc_pairs[i,'outcome'], ".coloc.results")
      coloc_results <- coloc_res$results %>% rename(SNP = snp) %>% arrange(desc(SNP.PP.H4))
      
      write.table(coloc_results, file = coloc_results_file, quote = FALSE, sep = "\t", row.names = FALSE)
      coloc_combine_res <- data.frame(exposure = coloc_pairs[i,'exposure'], outcome = coloc_pairs[i,'outcome'],
                                      lead_SNP = coloc_results[1,"SNP"], t(as.data.frame(coloc_res$summary)))
      return(coloc_combine_res)
    }
    
    
    cl <- makeCluster(8)
    registerDoParallel(cl)
    clusterEvalQ(cl, .libPaths("/data/slurm/wanghc/R"))
    
    res <- foreach(i = 1:nrow(coloc_pairs), .combine = rbind, .packages = c("tidyverse", "coloc"), .errorhandling = "pass") %dopar% {
      coloc_abf(i)
    }
    write.table(res, file =  print(paste0("/data/slurm/wanghc/microbiome_QTL/Coloc_microbiome_gut/Coloc_Microbiome-",a,"_eQTL/","Coloc_", h,"-",a,"_eQTL", "/","Coloc_", h,"-Sigmoid_eQTL",".chr",j,".coloc")), quote = FALSE, sep = "\t", row.names = FALSE)
    
    stopImplicitCluster()
    stopCluster(cl)}}}
