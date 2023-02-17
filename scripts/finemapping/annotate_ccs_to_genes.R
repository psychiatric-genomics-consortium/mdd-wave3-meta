### Script: Annotate locus credible causal sets with genes
### Author: JRIC
### Date: 2023-02-17

### Run as Rscript annotate_ccs_to_genes.R [Credible Causal Set Boundaries file] [Reduced Gene Matrix file] [Output file]

### Strings aren't factors

options(stringsAsFactors=F)

### Initialise args

args <- commandArgs(trailingOnly=TRUE)

### Install packages - only run once

# install.packages("data.table")
# install.packages("BiocManager")
# BiocManager::install("GenomicRanges")

### Load packages

library(data.table)
library(GenomicRanges)

### Load data

## Expected format of Boundaries file - no header, columns for Locus_ID, CHR, CCS_Start_BP, CCS_End_BP in that order

Boundaries <-fread(args[1], data.table=F, head=F)
names(Boundaries) <- c("locus", "chr", "start", "end")

## Expected format for reduced Gene Matrix file - header: gene_id [ENSEMBL gene ID], gene_name [HUGO gene name], hg19g0 [hg19 Chr], g1 [hg19 gene start], g2 [hg19 gene end]
## Original file from Pat Sullivan: https://figshare.com/articles/dataset/geneMatrix/13335548

Genes <- fread("args[2]", data.table=F)

### Construct Genomic Ranges objects

Genes_GR <- GRanges(
    seqnames = Rle(Genes$hg19g0),
    ranges = IRanges(start = Genes$g1, end = Genes$g2), 
    gene_name = Genes$gene_name, 
    gene_id = Genes$gene_id)

Boundaries_GR <- GRanges(
    seqnames = Rle(Boundaries$chr),
    ranges = IRanges(start=Boundaries$start, end=Boundaries$end), 
    locus = Boundaries$locus)
    
### Identify overlaps of Boundaries and Genes files
## Any overlap of a boundary with the region of a gene --> Locus
## The boundary is fully within the region of a gene --> High Confidence (HC) Locus

Boundaries_Genes_Overlaps_Extended <- findOverlaps(Boundaries_GR, Genes_GR, type="any")
Boundaries_Genes_Overlaps_HC_Extended <- findOverlaps(Boundaries_GR, Genes_GR, type="within")

### Extract rows with overlaps from Boundaries and Genes objects and column bind to a data frame

Boundaries_Genes_Locus <- as.data.frame(cbind(as.data.frame(Boundaries_GR[c(Boundaries_Genes_Overlaps_Extended@from), ]), as.data.frame(Genes_GR[c(Boundaries_Genes_Overlaps_Extended@to), ])))
Boundaries_Genes_HC_Locus <- as.data.frame(cbind(as.data.frame(Boundaries_GR[c(Boundaries_Genes_Overlaps_HC_Extended@from), ]), as.data.frame(Genes_GR[c(Boundaries_Genes_Overlaps_HC_Extended@to), ])))

### Add new indicator columns to Boundaries for final output

Boundaries$HC.Locus.Gene.Names <- NA
Boundaries$HC.Locus.Gene.IDs <- NA
Boundaries$Locus.Gene.Names <- NA
Boundaries$Locus.Gene.IDs <- NA

### Fill indicator columns

for(i in Boundaries$locus){
    Boundaries$HC.Locus.Gene.Names[i] <- ifelse(is.na(match(i, Boundaries_Genes_HC_Locus$locus)), NA, paste(Boundaries_Genes_HC_Locus[Boundaries_Genes_HC_Locus$locus == i, "gene_name"], collapse=","))
    Boundaries$HC.Locus.Gene.IDs[i] <- ifelse(is.na(match(i, Boundaries_Genes_HC_Locus$locus)), NA, paste(Boundaries_Genes_HC_Locus[Boundaries_Genes_HC_Locus$locus == i, "gene_id"], collapse=","))
    Boundaries$Locus.Gene.Names[i] <- ifelse(is.na(match(i, Boundaries_Genes_Locus$locus)), NA, paste(Boundaries_Genes_Locus[Boundaries_Genes_Locus$locus == i, "gene_name"], collapse=","))
    Boundaries$Locus.Gene.IDs[i] <- ifelse(is.na(match(i, Boundaries_Genes_Locus$locus)), NA, paste(Boundaries_Genes_Locus[Boundaries_Genes_Locus$locus == i, "gene_id"], collapse=","))
}

### Write output

fwrite(Boundaries, args[3], col.names=T, row.names=F, quote=F, na="NA", sep="\t")

