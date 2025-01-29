####find associated eQTL based on MR result  tissue type: transverse sigmoid ileum 
Transverse_expo_LDSC <- fread(file = paste0("/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/Colon-microbiome/eQTL-microbiome/Colon_Transerve.eQTL__expo_new/Transverse_expo_LDSC.txt"), sep = "\t", header = FALSE, 
                         col.names = c("Phenotype", "SNP", "beta", "pval","se","eaf","samplesize","effect_allele", "other_allele"))


MR_data <- fread(file = paste0("/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/Transverse-microbiome_genus/final_MR/Transverse_eQTL-genus.microbiota.result.sorted.MR"), sep = "\t", header = FALSE, 
                         col.names = c("exposure",       "outcome", "method",  "nsnp",   " b" ,     " se ",     "pval",   " adjusted_pval",   "gene_name"))

MR_signif_data <- MR_data %>%
  select(exposure) %>% 
  distinct()



filtered_data <- Transverse_expo_LDSC %>%
  filter(Phenotype %in% MR_signif_data$exposure)
nrow(MR_signif_data)

write.table(filtered_data,paste0("/data/slurm/wanghc/microbiome_QTL/LDscore/Transverse_microbiome_expo_bed/Transverse_expo_LDSC.txt"),quote=F,sep="\t",col.names=T,row.names=F)



######convert files to bed files

vcf = read.table('/sugon/nanjh/scAPA_project/01genotype_data/liftover/10WGSfromBroad_match_samples_merged_chrall_qc_hg38_dedup.vcf.gz',head=F)
vcf = fread(file = '/sugon/nanjh/scAPA_project/01genotype_data/liftover/10WGSfromBroad_match_samples_merged_chrall_qc_hg38_dedup.vcf.gz',sep = "\t", header = TRUE)
vcf = vcf[,c(1:5)]
colnames(vcf) = c('chr','pos','rsid','Ref','Alt')



####ct.list = c('Ileum','Sigmoid','Transverse')

for(i in 1:22){
  
  qtl = Phe1_dat_test <- fread(file = paste0("/data/slurm/wanghc/microbiome_QTL/LDscore/Transverse_microbiome_expo_bed/Transverse_expo_LDSC.txt"), sep = "\t", header = FALSE, 
                               col.names = c("Phenotype", "SNP", "beta", "pval","se","eaf","samplesize","effect_allele", "other_allele")) 
  qtl = merge(qtl,vcf,by.x='SNP',by.y='rsid')
  qtl$off = qtl$pos - 1
  qtl$chr = paste0('chr',qtl$chr)
  qtl = qtl[,c('chr','off','pos')]
  qtl = qtl[!duplicated(qtl),]  # de-duplicate
  write.table(qtl,paste0("/data/slurm/wanghc/microbiome_QTL/LDscore/Transverse_microbiome_expo_bed/","Transverse.signif.bed"),quote=F,sep="\t",col.names=F,row.names=F)
}



######split the files 

split_data <- qtl%>%
  group_split(chr)
for (i in seq_along(split_data)) {
  # Extract the chromosome value from the first row of the current subset
  chr_value <- unique(split_data[[i]]$chr)
  
  # Construct the filename based on the chromosome value
  filename <- paste0("/data/slurm/wanghc/microbiome_QTL/LDscore/Transverse_microbiome_expo_bed/Transverse_microbiome_expo_bed/Transverse.", chr_value, ".signif.bed")
  
  # Save the current subset to a file
  write.table(split_data[[i]], filename, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}



#######exend the file 

#ct.list = c('Ast','Exc','Inh','Oli','Mic','Opc')

exts = c(1000,10000) # try extend for 1k and 10k

for(i in 1:22){
  qtl = fread(paste0("/data/slurm/wanghc/microbiome_QTL/LDscore/Ileum_microbiome_expo_bed/Ileum_microbiome_expo_bed/Ileum.chr",i,".signif.hg19.bed"),sep = "\t", header = FALSE)
  for(ext in exts){
    qtl$left = qtl$V3 - ext
    qtl$right = qtl$V3 + ext
    write.table(qtl[,c('V1','left','right')],paste0("/data/slurm/wanghc/microbiome_QTL/LDscore/Ileum_microbiome_expo_bed/Ileum_microbiome_expo_bed/Ileum.chr",i,'.signif.hg19.ext',ext,'.bed'),quote=F,sep="\t",col.names=F,row.names=F)
  }
}


######## convert the file to hg19 version
/data/slurm/software/ucsc/liftOver $ct.apaQTL.signif.bed /data/slurm/public/Reference/human/hg38ToHg19.over.chain.gz $ct.apaQTL.signif.hg19.bed $ct.unmap


####### produce annotation files 

for s in $(seq 1 22)
do
/data/slurm/wanghc/ldsc/make_annot.py --bed-file /data/slurm/wanghc/microbiome_QTL/LDscore/Ileum_expo_bed/Ileum.chr$s.signif.hg19.ext1000.bed --bimfile /data/slurm/xiongxs/PublicData/LDSC/1000G_EUR_Phase3_plink/1000G.EUR.QC.$s.bim  --annot-file /data/slurm/wanghc/microbiome_QTL/LDscore/annotation/Ileum.chr$s.signif.hg19.ext1000.annot.gz   1>Ileum.$s.out 2>Ileum.$s.err & 
  done


########ldscore 

for a in Transverse Ileum Sigmoid
do


for s  in $(seq 1 22)
do
python    /data/slurm/wanghc/ldsc/ldsc.py   --l2 --bfile /data/slurm/xiongxs/PublicData/LDSC/1000G_EUR_Phase3_plink/1000G.EUR.QC.$s  --ld-wind-cm 1 --annot annotation/$a.chr$s.signif.hg19.ext1000.annot.gz --thin-annot --out    annotation/$a.signif.hg19.ext1000.chr$s  --print-snps /data/slurm/xiongxs/PublicData/LDSC/listHM3.txt 1>$a.$s.out 2>$a.$s.err & 
  done
wait
done

########process the outcome files

samples <-  fread(file = "/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/gut_microbiome_genus/genus.txt",sep = "\t", header = FALSE)
for (h in 1:nrow(samples)) {
  
  
  for (j in 1:22) {
    bacteria <- samples$V1[h]
Phe2.dat <- paste0("/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/gut_microbiome_genus/",bacteria,"/",bacteria,"_chr", j, ".1.txt")

Phe2_dat_test<- fread(file = Phe2.dat, sep = "\t", header = FALSE, 
                      col.names = c("Phenotype", "SNP", "beta", "effect_allele", "other_allele", "pval", "se", "samplesize"))
Phe2_dat_test$Z <- Phe2_dat_test$beta / abs(Phe2_dat_test$beta) * abs(qnorm(Phe2_dat_test$pval/2))

selected_data <- Phe2_dat_test %>%
  select(SNP, A1 = effect_allele, A2 = other_allele, Z, N = samplesize)


write.table(selected_data, file = print(paste0("/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/gut_microbiome_genus/",bacteria,"/",bacteria,"_chr", j, ".sumstats")), quote = FALSE, sep = "\t", row.names = FALSE)}}


atlas<-  fread(file = "/data/slurm/wanghc/microbiome_QTL/moloc_gut_microbiome/moloc_Ileum_microbiome_genus/moloc_Ileum-eQTL_microbiome_genus_diseases/atlas.txt",sep = "\t", header = FALSE)
gwas_meta_file <-paste0 ("/data/slurm/wanghc/microbiome_QTL/moloc_gut_microbiome/diseases/ukb.info.table.new.txt")
gwas_meta <- fread(gwas_meta_file, sep = "\t", header = TRUE, col.names = c("ID", "Description", "category", "pops", "num_pops", "Ncases","n_controls_EUR","filename","N"))


str(atlas)
atlas <- as.data.frame(atlas)
for (h in 1:nrow(atlas)) {
  
  
 # for (j in 1:22) {
    atlas_data <- atlas$V1[h]
    
    
    gwas_sum_file <- paste0("/data/slurm/wanghc/microbiome_QTL/moloc_gut_microbiome/diseases/ieu_atlas/",atlas_data,"/",atlas_data,".sumstats.gz")
    
    
    gwas_meta <- fread(gwas_meta_file, sep = "\t", header = FALSE, col.names = c("ID", "Ncases", "Nctrl", "N", "Type", "Description"))
    
    gwas_dat <- fread(gwas_sum_file, sep = "\t", header = TRUE, col.names = c("SNP", "A1", "A2", "BETA", "SE", "MAF", "PVAL"))
    
    gwas_dat$Z <- gwas_dat$BETA / abs(gwas_dat$BETA) * abs(qnorm(gwas_dat$PVAL/2))
    
   gwas_dat <- gwas_dat %>% mutate(N = gwas_meta[ID == atlas_data, N]) %>%dplyr::select(SNP, A1, A2, Z, N)

    
 
    
   # selected_data <- Phe2_dat_test %>%
    #  select(SNP, A1 = effect_allele, A2 = other_allele, Z, N = samplesize)
    
    
    write.table(gwas_dat, file = print(paste0("/data/slurm/wanghc/microbiome_QTL/LDscore/diseases/",atlas_data,".sumstats")), quote = FALSE, sep = "\t", row.names = FALSE)}


gwas_meta_file <-paste0 ("/data/slurm/wanghc/microbiome_QTL/moloc_gut_microbiome/diseases/ukb.info.table.new.txt")
gwas_meta <- fread(gwas_meta_file, sep = "\t", header = TRUE, col.names = c("ID", "Description", "category", "pops", "num_pops", "Ncases","n_controls_EUR","filename","N"))


ukb<-  fread(file = "/data/slurm/wanghc/microbiome_QTL/moloc_gut_microbiome/moloc_Ileum_microbiome_genus/moloc_Ileum-eQTL_microbiome_genus_diseases/ukb.txt",sep = "\t", header = FALSE)

ld <- function(i){ ukb_data <- ukb$V1[41]


gwas_sum_file <- paste0("/data/slurm/wanghc/microbiome_QTL/moloc_gut_microbiome/diseases/ukb/",ukb_data,"/",ukb_data,".tsv.gz")


gwas_dat <- fread(gwas_sum_file, sep = "\t", header = TRUE, col.names = c("SNP", "BETA","SE","PVAL","MAF","A1", "A2"))

gwas_dat <- fread(gwas_sum_file, sep = "\t", header = TRUE, col.names = c(
  "SNP","A1", "A2",    "BETA",    "SE",     " MAF",     "PVAL",   "N",      "Ncases"
))
gwas_dat$A1 <- toupper(gwas_dat$A1)
gwas_dat$A2 <- toupper(gwas_dat$A2)
gwas_dat$Z <- gwas_dat$BETA / abs(gwas_dat$BETA) * abs(qnorm(gwas_dat$PVAL/2))

gwas_dat <- gwas_dat%>% mutate(N = gwas_meta[ID == ukb_data, N])  %>%dplyr::select(SNP, A1, A2, Z, N)




# selected_data <- Phe2_dat_test %>%
#  select(SNP, A1 = effect_allele, A2 = other_allele, Z, N = samplesize)


write.table(gwas_dat, file = print(paste0("/data/slurm/wanghc/microbiome_QTL/LDscore/diseases/",ukb_data,".sumstats")), quote = FALSE, sep = "\t", row.names = FALSE)}


gwas_meta_file <-paste0 ("/data/slurm/licy/QTL_integration/1.Database/GWAS/GWAS_atlas/atlas_sp.list")


gwas_meta <- fread(gwas_meta_file, sep = "\t", header = FALSE, col.names = c("ID", "Ncases", "Nctrl", "N", "Type", "Description"))


ld<- function(i){ atlas_data <- atlas$V1[i]

gwas_sum_file <- paste0("/data/slurm/wanghc/microbiome_QTL/moloc_gut_microbiome/diseases/ieu_atlas/",atlas_data,"/",atlas_data,".sumstats.gz")



gwas_dat <- fread(gwas_sum_file, sep = "\t", header = TRUE, col.names = c("SNP", "A1", "A2", "BETA", "SE", "MAF", "PVAL"))

gwas_dat$A1 <- toupper(gwas_dat$A1)
gwas_dat$A2 <- toupper(gwas_dat$A2)

#gwas_dat<- gwas_dat[!is.na(gwas_dat$PVAL), ]
gwas_dat$PVAL <- as.numeric(as.character(gwas_dat$PVAL))
gwas_dat$Z <- gwas_dat$BETA / abs(gwas_dat$BETA) * abs(qnorm(gwas_dat$PVAL/2))

gwas_dat <- gwas_dat%>% mutate(N = gwas_meta[ID == atlas_data , N])  %>%dplyr::select(SNP, A1, A2, Z, N)




# selected_data <- Phe2_dat_test %>%
#  select(SNP, A1 = effect_allele, A2 = other_allele, Z, N = samplesize)


write.table(gwas_dat, file = print(paste0("/data/slurm/wanghc/microbiome_QTL/LDscore/diseases/",atlas_data,".sumstats")), quote = FALSE, sep = "\t", row.names = FALSE)}



ieu<-  fread(file = "/data/slurm/wanghc/microbiome_QTL/moloc_gut_microbiome/moloc_Ileum_microbiome_genus/moloc_Ileum-eQTL_microbiome_genus_diseases/ieu.txt",sep = "\t", header = FALSE)



ld <- function(i){ ieu_data <- ieu$V1[i]

gwas_sum_file <- paste0("/data/slurm/wanghc/microbiome_QTL/moloc_gut_microbiome/diseases/ieu_atlas/",ieu_data,"/",ieu_data,".sumstats.gz")



gwas_dat <- fread(gwas_sum_file, sep = "\t", header = TRUE, col.names = c("SNP", "A1", "A2", "BETA", "SE", "MAF", "PVAL", "N"))

gwas_dat$A1 <- toupper(gwas_dat$A1)
gwas_dat$A2 <- toupper(gwas_dat$A2)

#gwas_dat<- gwas_dat[!is.na(gwas_dat$PVAL), ]
gwas_dat$PVAL <- as.numeric(as.character(gwas_dat$PVAL))
gwas_dat$Z <- gwas_dat$BETA / abs(gwas_dat$BETA) * abs(qnorm(gwas_dat$PVAL/2))

gwas_dat <- gwas_dat %>%dplyr::select(SNP, A1, A2, Z, N)




# selected_data <- Phe2_dat_test %>%
#  select(SNP, A1 = effect_allele, A2 = other_allele, Z, N = samplesize)


write.table(gwas_dat, file = print(paste0("/data/slurm/wanghc/microbiome_QTL/LDscore/diseases/",ieu_data,".sumstats")), quote = FALSE, sep = "\t", row.names = FALSE)}


cl <- makeCluster(8)
registerDoParallel(cl)
clusterEvalQ(cl, .libPaths("/data/slurm/wanghc/R"))
res <- foreach(i = 1:nrow(ieu), .packages = c("data.table", "dplyr"), .errorhandling = "pass") %dopar% {
ld(i)
  
}



stopImplicitCluster()
stopCluster(cl)

#######merge the file

samples <-  fread(file = "/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/gut_microbiome_genus/genus.txt",sep = "\t", header = FALSE)
for (h in 1:nrow(samples)) {
data_list <- list()
for (i in 1:22) {
  
  bacteria <- samples$V1[h]
microbiome_qtl <- fread(file = print(paste0("/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/gut_microbiome_genus/",bacteria,"/",bacteria,"_chr", i, ".sumstats")) , sep = "\t", header = TRUE) 

data_list[[i]] <- microbiome_qtl
}

merged_data <- rbindlist(data_list)

write.table(merged_data, file = print(paste0("/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/gut_microbiome_genus/",bacteria,"/",bacteria, ".sumstats")), quote = FALSE, sep = "\t", row.names = FALSE)

}

#######concatenate the results

cts.list = c('Transverse','Ileum','Sigmoid')

samples <-  fread(file = "/data/slurm/wanghc/microbiome_QTL/LDscore/diseases.txt",sep = "\t", header = FALSE)


samples_data = as.character(samples$V1)

df.merge = data.frame()
for(ct in cts.list){
  for(trait in samples_data){
    print(ct)
    print(trait)
    df = read.table(paste0("/data/slurm/wanghc/microbiome_QTL/LDscore/LDSC_out_MR_diseases_ext1k/",ct,"/",ct,".",trait,".baseline_v2.2.coeff.results"),head=T)
    df = df[1,]
    df$diseases_or_traits = trait
    df$tissue = ct
    df.merge = rbind(df.merge,df)
  }
}

saveRDS(df.merge,'/data/slurm/wanghc/microbiome_QTL/LDscore/LDSC_out_MR_diseases_ext1k/MR_diseases.LDSC_output_ext1k.rds')
write.table(df.merge, file = print(paste0("/data/slurm/wanghc/microbiome_QTL/LDscore/LDSC_out_MR_diseases_ext1k/MR_diseases.LDSC_output_ext1k.txt")), quote = FALSE, sep = "\t", row.names = FALSE)





LDSC_out_ext10k <- fread(file = print(paste0("/data/slurm/wanghc/microbiome_QTL/LDscore/LDSC_out_MR_diseases_ext10k/MR_diseases.LDSC_output_ext10k.txt")) , sep = "\t", header = TRUE) 



filtered_df <- LDSC_out_ext10k %>%
  filter(Enrichment > 1 & Enrichment_p < 0.05)

LDSC_out_ext10k$adjusted_pval <- p.adjust(LDSC_out_ext10k$Enrichment_p, method="fdr")

#significant_results <- LDSC_out_ext10k[LDSC_out_ext10k$adjusted_pval < 0.1, ]
#filtered_df <- significant_results[significant_results$Enrichment>1]

filtered_df <- LDSC_out_ext10k %>%
  filter(Enrichment > 1 & adjusted_pval < 0.1)



write.table(filtered_df, file = print(paste0("/data/slurm/wanghc/microbiome_QTL/LDscore/LDSC_out_MR_diseases_ext1k/MR_diseases.LDSC_output_ext1k.filtered.txt")), quote = FALSE, sep = "\t", row.names = FALSE)

ukb_meta_file <-paste0 ("/data/slurm/wanghc/microbiome_QTL/moloc_gut_microbiome/diseases/ukb.info.table.new.txt")
ukb_meta <- fread(ukb_meta_file, sep = "\t", header = TRUE, col.names = c("ID", "Description", "category", "pops", "num_pops", "Ncases","n_controls_EUR","filename","N"))


merged_data <- filtered_df  %>%
  left_join(ukb_meta, by = c("diseases_or_traits" = "ID"))


merged_data$diseases_or_traits <- coalesce(merged_data$Description, merged_data$diseases_or_traits)

merged_data <- merged_data %>% dplyr::select(-Description,-category, -pops, -num_pops, -Ncases,-n_controls_EUR,-filename,-N)

ieu_meta_file <- paste0("/data/slurm/licy/QTL_integration/4.Moloc/GWAS.metadata_S1")
ieu_meta <- fread(ieu_meta_file, sep = "\t", header = TRUE) %>% mutate(Ncase = gsub(",", "", Ncase)) %>% mutate_at("Ncase", as.numeric)

merged_data <- merged_data %>%
  left_join(ieu_meta, by = c("diseases_or_traits" = "ID"))

merged_data$diseases_or_traits <- coalesce(merged_data$"Trait tag", merged_data$diseases_or_traits)
merged_data <- merged_data %>% dplyr::select(-Decription,-Consortium, -Population, -Category.y, -Ncase,-Ncontrol,-Nsnps,-N,-Sex,-PMID,-Year,-URL,-Remarks,-Type,-`Trait tag`)

write.table(merged_data, file = print(paste0("/data/slurm/wanghc/microbiome_QTL/LDscore/LDSC_out_MR_diseases_ext1k/MR_diseases.LDSC_output_ext1k.filtered_des.txt")), quote = FALSE, sep = "\t", row.names = FALSE)





