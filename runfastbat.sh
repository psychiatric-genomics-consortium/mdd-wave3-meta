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
