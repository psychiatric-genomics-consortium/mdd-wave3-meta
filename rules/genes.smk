## GeneMatrix.tsv.gz downloaded from https://figshare.com/articles/dataset/geneMatrix/13335548 on 31st January 2022

## Remove double quotes, select only protein_coding genes (N = 19878) and those with known positions on hg19
zcat geneMatrix.tsv.gz | sed 's/\"//g' | awk '($5 == "protein_coding" && $12 ~ /^chr/)' | awk -v range=3 '{print substr($12,range+1),$13,$14,$3}' > geneMatrixForFastBAT.txt

## Remove double quotes, select all genes (N = 58844) with known positions on hg19
zcat geneMatrix.tsv.gz | sed 's/\"//g' | awk '($12 ~ /^chr/)' | awk -v range=3 '{print substr($12,range+1),$13,$14,$3}' > AllgeneMatrixForFastBAT.txt


## Extract header row and apply MAF >= 0.01 and infoscore > 0.8 to summary stats

zcat daner_pgc_mdd_full_eur_hg19_v3.49.24.11.neff.gz | head -n 1 | awk '{print $2, $11}' > mdd_fastbat_neff_maf01_info80.txt
zcat daner_pgc_mdd_full_eur_hg19_v3.49.24.11.neff.gz | awk '($7 >= 0.01 && $7 <= 0.99 && $8 >= 0.8)' | awk '{print $2, $11}' >> mdd_fastbat_neff_maf01_info80.txt


## Run fastBAT within GCTA using qc'd summary stats and geneMatrix protein coding gene list (N = 19878)

./gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta_v1.94.0Beta_linux_kernel_3_x86_64_static \
--bfile all_phase3.eur \
--fastBAT mdd_fastbat_neff_maf01_info80.txt \
--fastBAT-gene-list geneMatrixForFastBAT.txt \
--out mdd_fastbat_geneMatrix


## Run fastBAT within GCTA using qc'd summary stats and all geneMatrix gene list (N = 58844)

./gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta_v1.94.0Beta_linux_kernel_3_x86_64_static \
--bfile all_phase3.eur \
--fastBAT mdd_fastbat_neff_maf01_info80.txt \
--fastBAT-gene-list AllgeneMatrixForFastBAT.txt \
--out mdd_fastbat_AllgeneMatrix
