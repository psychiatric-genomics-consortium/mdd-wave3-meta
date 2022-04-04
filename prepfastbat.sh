## Extract header row and apply MAF >= 0.01 and infoscore > 0.8 to summary stats

zcat daner_pgc_mdd_full_eur_hg19_v3.49.24.11.neff.gz | head -n 1 | awk '{print $2, $11}' > mdd_fastbat_neff_maf01_info80.txt
zcat daner_pgc_mdd_full_eur_hg19_v3.49.24.11.neff.gz | awk '($7 >= 0.01 && $7 <= 0.99 && $8 >= 0.8)' | awk '{print $2, $11}' >> mdd_fastbat_neff_maf01_info80.txt
