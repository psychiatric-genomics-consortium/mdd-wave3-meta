#!/usr/bin/Rscript

###################
# PLINK format 1000 genomes phase 3 data
###################
# Based on instructions from Hannah Meyer https://cran.r-project.org/web/packages/plinkQC/vignettes/Genomes1000.pdf
# With adaptations from Joni Coleman's reference resource

output_dir<-'resources/twas/1kg'
dir.create(output_dir)

# Download the 1000Genomes data provided by PLINK (PLINK 2 format)
system(paste0('wget -q -o /dev/null -O ',output_dir,'/all_phase3.pgen.zst https://www.dropbox.com/s/afvvf1e15gqzsqo/all_phase3.pgen.zst?dl=1'))
system(paste0('wget -q -o /dev/null -O ',output_dir,'/all_phase3.pvar.zst https://www.dropbox.com/s/op9osq6luy3pjg8/all_phase3.pvar.zst?dl=1'))
system(paste0('wget -q -o /dev/null -O ',output_dir,'/all_phase3.psam https://www.dropbox.com/s/yozrzsdrwqej63q/phase3_corrected.psam?dl=1'))

# Decompress the pgen file
system(paste0('plink2 --zst-decompress ',output_dir,'/all_phase3.pgen.zst > ',output_dir,'/all_phase3.pgen'))

# Convert to plink 1 format
system(paste0('plink2 --pfile ',output_dir,'/all_phase3 vzs --max-alleles 2 --make-bed --out ',output_dir,'/all_phase3'))

# Delete the plink2 format data
system(paste0('rm ',output_dir,'/all_phase3.pgen*'))
system(paste0('rm ',output_dir,'/all_phase3.pvar*'))

# Make FID == IID
fam<-read.table(paste0(output_dir,'/all_phase3.fam'), header=F)
fam$V1<-fam$V2
write.table(fam, paste0(output_dir,'/all_phase3.fam'), col.names=F, row.names=F, quote=F)

# Split the plink data by chromosome (autosome only)
for(chr in 1:22){
  system(paste0('plink --bfile ',output_dir,'/all_phase3 --chr ',chr,' --allow-extra-chr --make-bed --out ',output_dir,'/all_phase3.chr',chr))
}

# Delete genome wide data
system(paste0('rm ',output_dir,'/all_phase3.bed'))
system(paste0('rm ',output_dir,'/all_phase3.bim'))
system(paste0('rm ',output_dir,'/all_phase3.fam'))
system(paste0('rm ',output_dir,'/all_phase3.log'))

# For autosomes: Remove variants with 3+ alleles or multiple positions
for(chr in 1:22){
  system(paste0('plink --bfile ',output_dir,'/all_phase3.chr',chr,' --bmerge ',output_dir,'/all_phase3.chr',chr,' --merge-mode 6 --make-bed --out ',output_dir,'/all_phase3.chr',chr,'_noDup'))
}

missnp<-NULL
for(chr in 1:22){
  log_chr<-readLines(paste0(output_dir, '/all_phase3.chr',chr,'_noDup.log'))
  log_chr<-log_chr[grepl('Warning: Multiple positions seen for variant',log_chr)]
  log_chr<-gsub("'.","",gsub("Warning: Multiple positions seen for variant '", "", log_chr))
  print(chr)
  print(log_chr)
  missnp<-c(missnp, log_chr)
}
write.table(missnp, paste0(output_dir, '/all_phase3_noDup.missnp'), col.names=F, row.names=F, quote=F)

for(chr in 1:22){
  system(paste0('plink --bfile ',output_dir,'/all_phase3.chr',chr,' --exclude ',output_dir,'/all_phase3_noDup.missnp --make-bed --out ',output_dir,'/all_phase3.cleaned.chr',chr))
}

system(paste0('rm ',output_dir,'/*.missnp'))
system(paste0('rm ',output_dir,'/*.log'))

for(i in 1:22){
  system(paste0('mv ',output_dir,'/all_phase3.cleaned.chr',i,'.bed ', output_dir,'/all_phase3.chr',i,'.bed'))
  system(paste0('mv ',output_dir,'/all_phase3.cleaned.chr',i,'.bim ', output_dir,'/all_phase3.chr',i,'.bim'))
  system(paste0('mv ',output_dir,'/all_phase3.cleaned.chr',i,'.fam ', output_dir,'/all_phase3.chr',i,'.fam'))
}

####################
# Prepare keep files for super_populations and populations
####################

# Create a keep file listing each population super population from the reference.
system(paste0('mkdir -p ',output_dir,'/super_pop_keep_files'))

pop_data<-read.table(paste0(output_dir, '/all_phase3.psam'), header=F, stringsAsFactors=F)

for(i in unique(pop_data$V5)){
  write.table(cbind(pop_data$V1[pop_data$V5 == i],pop_data$V1[pop_data$V5 == i]), paste0(output_dir,'/super_pop_keep_files/',i,'_samples.keep'), col.names=F, row.names=F, quote=F)
}

####################
# Subset EUR individuals from the 1KG Phase 3 data
####################

# Retain only variants with a MAF > 0.001
system(paste0('mkdir ',output_dir,'/EUR'))

for(i in 1:22){
  system(paste0('plink --bfile ', output_dir,'/all_phase3.chr',i,' --keep ',output_dir,'/super_pop_keep_files/EUR_samples.keep --maf 0.001 --make-bed --out ',output_dir,'/EUR/EUR_phase3.MAF_001.chr',i))
}

####################
# Compute allele frequencies across all individuals
####################

for(chr in 1:22){
  system(paste0('plink --bfile ',output_dir,'/EUR/EUR_phase3.MAF_001.chr',chr,' --freq --out ',output_dir,'/EUR/EUR_phase3.MAF_001.chr',chr))
}

# Delete the log and nosex files
system(paste0('rm ',output_dir,'/EUR/EUR_phase3.MAF_001.chr*.log'))

