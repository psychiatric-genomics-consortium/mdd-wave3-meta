
### Go to finemapping results

cd results/finemapping

### Get header

gunzip -c locus_results_full_eur_hg19_v3.49.24.11.rp.chr1.80715434_80893134.gz | head -1 > locus_results_full_eur_hg19_v3.49.24.11.rp.wg.results.table

### Get top PIP in each region

for file in locus_results_full_eur_hg19_v3.49.24.11.rp.chr*.gz
do
    gunzip -c $file | head -2 | tail -1 >> locus_results_full_eur_hg19_v3.49.24.11.rp.wg.results.table
done

### Get cojo SNPs

for file in locus_results_full_eur_hg19_v3.49.24.11.rp.chr*.gz
do
    LANG=C zgrep -wf <(awk 'NR > 1 {print $3}' ../../docs/tables/meta_snps_full_eur.cojo.format.txt) $file >> locus_results_full_eur_hg19_v3.49.24.11.rp.wg.results.table
done

### Sort and make unique

sort -k1,1 -k3,3 -g locus_results_full_eur_hg19_v3.49.24.11.rp.wg.results.table | uniq > TEMP
mv TEMP locus_results_full_eur_hg19_v3.49.24.11.rp.wg.results.table
