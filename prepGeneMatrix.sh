## GeneMatrix.tsv.gz downloaded from https://figshare.com/articles/dataset/geneMatrix/13335548 on 31st January 2022

## Remove double quotes, select only protein_coding genes and those with known positions on hg19
zcat geneMatrix.tsv.gz | sed 's/\"//g' | awk '($5 == "protein_coding" && $12 ~ /^chr/)' | awk -v range=3 '{print substr($12,range+1),$13,$14,$3}' > geneMatrixForFastBAT.txt