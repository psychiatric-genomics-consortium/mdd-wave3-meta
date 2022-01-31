## Run fastBAT within GCTA using qc'd summary stats and geneMatrix gene list

./gcta_v1.94.0Beta_linux_kernel_4_x86_64/gcta_v1.94.0Beta_linux_kernel_4_x86_64_static \
--bfile ldref/g1000_eur \
--fastBAT mdd_fastbat_neff_maf01_info80.txt \
--fastBAT-gene-list geneMatrixForFastBAT.txt \
--out mdd_fastbat_geneMatrix
