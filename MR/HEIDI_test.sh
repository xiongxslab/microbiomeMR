####Heidi test 
####produce ma file 
.libPaths("/data/slurm/wanghc/R")
library(data.table)
library(dplyr)
library(foreach)
library(tidyr)

args <- commandArgs(T)
microbiome <- args[1]
Chr_txt <-  fread(file = "/data/slurm/wanghc/microbiome_QTL/SMR/microbiome_esd/chr.txt",sep = "\t", header = FALSE)
all_data <- list()
for (i in 1:nrow(Chr_txt)) {
        chr<-Chr_txt$V1[i]
maf_file <- paste0("/data/slurm/licy/QTL_integration/1.Database/Ref_panel/GTEx_v8/plink_maf/",chr,".frq")
maf <- fread(maf_file, sep = "\t", header = TRUE) %>% dplyr::select(CHR, SNP, A1, A2, MAF) %>% rename(Chr = CHR, Freq = MAF)

Phe1.dat <- paste0("/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/gut_microbiome_genus/",microbiome,"/",microbiome,"-", chr, ".txt.2.smr.gz")
Phe_dat <-fread(file = Phe1.dat, sep = "\t", header = F,
      col.names = c("Phenotype", "SNP", "rsid", "Beta","effect_allele", "other_allele","p","se","samplesize")) %>%
  inner_join(y = maf, by = "SNP") %>%
  mutate(Beta = ifelse(effect_allele == A1, Beta, -Beta)) %>%
  select( SNP, A1, A2, Freq, Beta, se, p,samplesize)%>%  rename(freq = Freq, b = Beta, n = samplesize)
  Phe_esd_file <- paste0("/data/slurm/wanghc/microbiome_QTL/SMR/","microbiome_esd/genus/", microbiome,"-",chr,".ma")
  write.table(Phe_dat, file = Phe_esd_file, quote = FALSE, sep = "\t", row.names = FALSE)
  all_data[[chr]] <- Phe_dat
}

merged_data <- bind_rows(all_data)

# Save the merged data to a file
merged_file <- paste0("/data/slurm/wanghc/microbiome_QTL/SMR/microbiome_esd/", microbiome, ".ma")
write.table(merged_data, file = merged_file, quote = FALSE, sep = "\t", row.names = FALSE)

####produce epi file     
.libPaths("/data/slurm/wanghc/R") 
library(data.table)
library(dplyr)
library(foreach)
library(tidyr)
input_file <- "/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/study_wise_FDR_MR/Transverse_MR/Transverse_eQTL.filtered.pos.gene.MR"
output_file <- "/data/slurm/wanghc/microbiome_QTL/SMR/Sigmoid_eQTL_exp_besd_new/Sigmoid.epi"

# Read the input file
df <- fread(input_file, header = TRUE, sep = "\t")

# Optional: remove "chr" prefix from column 8 if needed
df[[8]] <- sub("^chr", "", df[[8]])

# Compute the mean of column 11 and 12
mean_pos <- (df[[11]] + df[[12]]) / 2

# Construct output data frame
output_df <- data.table(
  chr   = df[[8]],
  gene1 = df[[1]],
  zero  = 0,
  mean  = mean_pos,
  gene2 = df[[1]],
  strand = "+"
)%>% unique()

write.table(output_df, file = output_file, quote = FALSE, sep = "\t", row.names = FALSE,col.names = T)


df <- fread("/data/slurm/wanghc/microbiome_QTL/SMR/Transverse_eQTL_exp_besd_new/Transverse.epi",sep = "\t", header = T)
df$chr <- as.character(df$chr)
chr_list <- unique(df$chr)

output_dir <- "/data/slurm/wanghc/microbiome_QTL/SMR/Transverse_eQTL_exp_besd_new/"
dir.create(output_dir, showWarnings = FALSE)

for (chr_val in chr_list) {
  chr_df <- df[df$chr == chr_val, ]
  output_file <- file.path(output_dir, paste0("chr", chr_val, ".epi"))
  write.table(chr_df, file = output_file, sep = "\t",quote = FALSE,  row.names = FALSE,col.names = F)
}

#####produce esi file
for i in `seq 22`; do 
awk 'BEGIN {OFS="\t"} NR==FNR {a[$1]; next} $1 in a' /data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/study_wise_BH_MR/Transverse_MR/chr${i}.txt   /data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/Colon-microbiome/eQTL-microbiome/Colon_Transerve.eQTL__expo/chr${i}.txt  |awk 'BEGIN {OFS="\t"} {
    split($2, b, ":");  sub(/^chr/, "", b[1]); 
    print $0, b[1], b[2], b[3], b[4]
}'  | awk 'BEGIN {OFS="\t"} {
    print $11, $2, 0, $12, $8, $9, $6
}' >/data/slurm/wanghc/microbiome_QTL/SMR/Transverse_eQTL_exp_besd_new/chr${i}.esi

done
#####produce flist file
for i in `seq 22`; do awk -v tis=$tis -v OFS='\t' 'BEGIN{print "Chr","ProbeID","GeneticDistance","ProbeBp","Gene","Orientation","PathOfEsd"}{print $0,"/data/slurm/wanghc/microbiome_QTL/SMR/Transverse_eQTL_exp_esd_new/"$2".esd"}' chr${i}.epi > chr${i}.flist ;  done 
##### produce esd file 

.libPaths("/data/slurm/wanghc/R")
library(data.table)
library(dplyr)
library(foreach)
library(tidyr)

args <- commandArgs(T)
chr <- args[1]

smr_pairs<- fread(file = paste0("/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/study_wise_BH_MR/Transverse_MR/",chr, ".txt"), sep = "\t", header = T
                            )%>% select(exposure, outcome)%>% distinct(exposure)
#View(smr_pairs)

maf_file <- paste0("/data/slurm/licy/xQTL_MR/1.Database/Ref_panel/GTEx_v8/plink_maf/",chr,".frq")
#maf <- fread(maf_file, sep = "\t", header = TRUE) %>% dplyr::select(SNP, A1, A2, MAF)
#maf <- fread(maf_file, sep = "\t", header = T) %>%  dplyr::select(Chr = CHR, Freq = MAF)
maf <- fread(maf_file, sep = "\t", header = TRUE) %>% dplyr::select(CHR, SNP, A1, A2, MAF) %>% rename(Chr = CHR, Freq = MAF)

Phe1.dat <- paste0("/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/Colon-microbiome/eQTL-microbiome/Colon_Transerve.eQTL__expo/", chr, ".txt")
Phe_dat <-fread(file = Phe1.dat, sep = "\t", header = T,
      col.names = c("Phenotype", "SNP", "Beta", "p","se","maf","samplesize","effect_allele", "other_allele","rsid")) %>%
  #mutate(Chr = gsub("chr", "", chr)) %>% 
  separate(SNP, c(NA, "Bp", NA, NA), ":", convert = TRUE, remove = FALSE) %>% 
 inner_join(y = maf, by = "SNP") %>%
 mutate(Beta = ifelse(effect_allele == A1, Beta, -Beta)) %>%
  #       z = case_when(Beta > 0 ~ qnorm(1-p/2), Beta < 0 ~ -qnorm(1-p/2)), se = Beta/z) %>%
  select(Phenotype, Chr, SNP, Bp, A1, A2, Freq, Beta, se, p)

for(i in 1:nrow(smr_pairs)){
  Phe_esd <- Phe_dat  %>% dplyr::filter(Phenotype == as.character(smr_pairs[i,"exposure"])) %>% select(-Phenotype)
  Phe_esd_file <- paste0("/data/slurm/wanghc/microbiome_QTL/SMR/","Transverse_eQTL_exp_esd_new/", smr_pairs[i,"exposure"], ".esd")
  write.table(Phe_esd, file = Phe_esd_file, quote = FALSE, sep = "\t", row.names = FALSE)
} 

microbiome_file="/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/gut_microbiome_genus/genus.txt"
microbiome_list=$(sed -n '1,116p' "$microbiome_file")

for microbiome in $microbiome_list; do for i in `seq 22`; do smr --bfile /data/slurm/licy/QTL_integration/1.Database/Ref_panel/GTEx_v8/plink/chr${i}  --beqtl-summary /data/slurm/wanghc/microbiome_QTL/SMR/Transverse_eQTL_exp_besd/chr${i}  
--gwas-summary  /data/slurm/wanghc/microbiome_QTL/SMR/microbiome_esd/${microbiome}.ma --peqtl-smr 0.99 --peqtl-heidi 0.99 --heidi-min-m 1 --out Transverse-${microbiome}/chr${i} --smr-multi --ld-multi-snp 0.01 --thread-num 5; done; done  
