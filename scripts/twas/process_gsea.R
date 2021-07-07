#!/usr/bin/Rscript

library(data.table)

# Read in competitive results
gsea<-fread('results/twas/gsea/PGC_MDD3_twas_AllTissues_GSEA.competitive.txt')

# There is nothing significant after fdr correction
# The top gene set is neurotrophin binding, a previously implicated set in depression
gsea_GO<-gsea[grepl('GO.', gsea$GeneSet),]

# There is apparent enrichment of brain sets but the non-independence of observations should be accoutned for.

# Compare with MESC results
mesc<-fread('results/twas/mesc/PGC_MDD3_TWAS.MESC.sets.categories.h2med')

# Extract GO sets for comparison with gsea results
mesc_GO<-mesc[grepl('GO_',mesc$Gene_category),]

# Convert mesc GO names to match gsea
mesc_GO$GeneSet<-gsub("[[:punct:]]", ".", mesc_GO$Gene_category)

# Merge results
both<-merge(gsea_GO, mesc_GO, by='GeneSet')

cor(both$h2med_enrichment, both$T)

both_enrich<-both[both$h2med_enrichment > 0,]
both_enrich<-both_enrich[order(both_enrich$h2med_enrichment_pvalue),]
both_enrich[,c('GeneSet','h2med_enrichment_pvalue','P')]

# The results are correlated 0.14, with many of the MESC enriched results showing enrichment in TWAS-GSEA. However, the most significant TWAS-GSEA sets are not significant in MESC. This indicates TWAS-GSEA is picking up signal that MESC is not. The sets found by TWAS-GSEA seem unrelated to depression, whereas MESC pathways are all neuronal. This may be due to MESC restricting its analysis to mediating genetic effects on expression, whereas TWAS-GSEA can be influenced by pleiotropy and linkage effects. One way of improving TWAS-GSEA might be to leverage COLOC results to remove signals that are driven by linkage.

####
# Split MESC results into pathways, broad categories and tissues
####

# Broad
mesc_broad<-mesc[(mesc$Gene_category %in% c('universe','gwascatalog','mgi_essential','CEGv2_subset_universe',"wang_essential_in_culture",'exac_lof_intolerant',"samocha_missense_constrained","cassa_selected_against_lof",'clinvar_path_likelypath','fda_approved_drug_targets')),]

mesc_broad$P.CORR<-p.adjust(mesc_broad$h2med_enrichment_pvalue, method='bonferroni')
mesc_broad<-mesc_broad[order(mesc_broad$h2med_enrichment_pvalue),]

write.csv(mesc_broad, 'results/twas/mesc/PGC_MDD3_TWAS.MESC.broad_categories.csv', row.names=F, quote=F)

# Heritability
mesc_h2<-mesc[(mesc$Gene_category %in% c('h2cis_bin_1','h2cis_bin_2','h2cis_bin_3','h2cis_bin_4',"h2cis_bin_5")),]

# Tissue set don't appear to be tested here

# There nearly 1600 pathway sets in this file, although the MESC paper restricts its analysis to ~800 sets. We can restrict our analysis to these sets using the supplementary table 7.

h2cis_bin_1
Adipose_Subcutaneous
mesc$Gene_category[grepl('GO_POSITIVE_REGULATION_OF_RESPONSE_TO_STIMULUS',mesc$Gene_category )]

