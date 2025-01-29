######Heidi test
microbiome_file="/data/slurm/wanghc/microbiome_QTL/MR_gut_microbiome/gut_microbiome_genus/genus.txt"
microbiome_list=$(sed -n '1,116p' "$microbiome_file")

for microbiome in $microbiome_list; do for i in `seq 22`; do smr --bfile /data/slurm/licy/QTL_integration/1.Database/Ref_panel/GTEx_v8/plink/chr${i}  --beqtl-summary /data/slurm/wanghc/microbiome_QTL/SMR/Ileum_eQTL_exp_besd/chr${i}  
--gwas-summary  /data/slurm/wanghc/microbiome_QTL/SMR/microbiome_esd/${microbiome}.ma --peqtl-smr 0.99 --peqtl-heidi 0.99 --heidi-min-m 1 --out Ileum-${microbiome}/chr${i} --smr-multi --ld-multi-snp 0.01 --thread-num 5; done; done  
