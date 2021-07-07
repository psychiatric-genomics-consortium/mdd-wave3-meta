#!/usr/bin/Rscript

###
# Create a table with colocalisation results for all significant features
###
library(data.table)
twas_sign <- fread("results/twas/twas_results/PGC_MDD3_twas_AllTissues_TWSig_CLEAN.txt")

twas_sign$Location<-paste0('chr',twas_sign$CHR,':',twas_sign$P0,'-',twas_sign$P1)   

# order CHR and P0
twas_sign_ordered <- twas_sign[order(twas_sign$CHR, twas_sign$P0), ]

twas_sign_ordered <- twas_sign_ordered[, c("Location", "ID", "PANEL_clean_short", "TWAS.Z", "TWAS.P", "COLOC.PP0", "COLOC.PP1", "COLOC.PP2", "COLOC.PP3", "COLOC.PP4")]

###
#Create a couple of additional columns specifying whether the feature is colocalised or not
###

#to specify coloc pp4 > 0.8 (see gusev et al (2019) Nat Genet on epithelial ovarian cancer)
twas_sign_ordered$Colocalised <- "No"
twas_sign_ordered$Colocalised[twas_sign_ordered$COLOC.PP4 > 0.8] <- "Yes"
cat(sum(twas_sign_ordered$Colocalised == "Yes"),'features present with PP4 greater than 0.8\n')
cat(length(unique(twas_sign_ordered$ID[twas_sign_ordered$Colocalised == "Yes"])),'unique genes present with PP4 greater than 0.8\n')

#to specify coloc pp3 < 0.2 
twas_sign_ordered$Low_PP3_0.2 <- "No"
twas_sign_ordered$Low_PP3_0.2[twas_sign_ordered$COLOC.PP3 < 0.2] <- "Yes"
cat(sum(twas_sign_ordered$Low_PP3_0.2 == "Yes"),'features present a PP3 greater than 0.2\n')

###
# Clean and Save 
###
twas_sign_ordered <- twas_sign_ordered[, c("Location", "ID", "PANEL_clean_short", "TWAS.Z", "TWAS.P", "COLOC.PP0", "COLOC.PP1", "COLOC.PP2", "COLOC.PP3", "COLOC.PP4", "Low_PP3_0.2", "Colocalised")]

write.csv(twas_sign_ordered, "results/twas/twas_results/PGC_MDD3_TWAS_colocalisation.csv", row.names = F, quote=F)
