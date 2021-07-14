. /etc/profile.d/modules.sh

module load igmm/apps/R/3.6.1

Rscript /exports/igmm/eddie/GenScotDepression/shen/Tools/PRSice/PRSice.R \
    --prsice /exports/igmm/eddie/GenScotDepression/shen/Tools/PRSice/PRSice_linux \
    --base $1 \
    --target $2 \
    --thread 5 \
    --stat OR \
    --pheno $3 \
    --pheno-col dummy_pheno \
    --clump-kb 250 \
    --clump-r2 0.25 \
    --bar-levels 0.00000005,0.0000001,0.000001,0.00001,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1 \
    --fastscore \
    --all-score \
    --out $4
