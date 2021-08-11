setwd('/exports/igmm/eddie/GenScotDepression/shen/SData/UKB/ii.PGRS/IDP_5thrls_40k_1.5k/PGCMDD3-meta-PheWAS/')
library(dplyr)
library(pbapply)

# Load data    -----------------------------------------------------------------------------------------------
fields.loose=read.csv('data/data_dictionary/fields_loose_for_process.csv',header=T,sep='\t',quote="'\"")

# Loose phenotypes   -----------------------------------------------------------------------------------------
# Function: interpolate by instances, keep interpolation record ***
extract_pheno <- function(field_dat,dat.dic,instances=0:2){
        field.id=field_dat$field.id
        master.dat=get(field_dat$dat)
            tmp.dic=filter(dat.dic,field_name==field.id)
        # get data **
        tmp.block.dat <- master.dat %>%
              select(matches(paste0('f.',field.id,'\\.'))) %>%
              select(matches(paste0('\\.',instances,'\\.', collapse="|")))

        # test if there are multiple rounds
        tmp.test.rounds <- tmp.block.dat %>%
              select(matches(paste0('f.',field.id,'.0.')))
        rounds.count = ncol(tmp.test.rounds)
        if (rounds.count>1){
            tmp.record=data.frame(c(-999,-999,-999))
            colnames(tmp.record)=paste0('f.',field.id)
            rownames(tmp.record)=paste0('instance.',instances)
            count.n.record=tmp.record
            tmp.dat=data.frame(pheno=rep(NA,nrow(master.dat)))
            colnames(tmp.dat)=paste0('f.',field.id)
        }else{
        # No multiple rounds
          tmp.dat=data.frame(tmp.block.dat[,1])
          colnames(tmp.dat)=paste0('f.',field.id)
          count.n.record<-data.frame(c(NA,NA,NA))
          rownames(count.n.record)=paste0('instance.',instances)
          colnames(count.n.record)=paste0('f.',field.id)
          for (i in 1:ncol(tmp.block.dat)){
            tmp.interpolate.dat=tmp.block.dat[,i]
            loc.to.interpolate=is.na(tmp.dat[,names(tmp.dat)])
            tmp.dat[,names(tmp.dat)][loc.to.interpolate]=tmp.interpolate.dat[loc.to.interpolate]
            # keep a record of N derived from interpolated data
            if (i==1){
                count.n.instance=sum(!is.na(tmp.interpolate.dat))
            }else{count.n.instance=sum(!is.na(tmp.interpolate.dat[loc.to.interpolate]))}
            count.n.record[i,1]=count.n.instance
          } 
        }
        tmp.dat$f.eid=master.dat$f.eid
        output.pheno=list('phenotype'=tmp.dat,'count'=count.n.record)
        return(output.pheno)  
}


# Load available data fields from files ***
# generate a list of files
g.path = '/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/ukb_phenotype/'
fs = list.files(path=g.path,recursive=T)%>% (function(x) x[grep('-ukb',x)]) %>% 
      (function(x) x[grep('\\.rds',x)])
rls <- fs %>%
        strsplit(.,'/') %>%
        lapply(function(x) x[1]) %>%
        unlist %>%
        strsplit(.,'-') %>%
        lapply(function(x) x[length(x)]) %>%
        unlist
        
fs=data.frame(file=fs,rls=rls,stringsAsFactors=F)
fs.list=split(fs, seq(nrow(fs)))
# extracting available fields in the files
extract_available_fields <- function(fname_rls,dat.dic,path='/exports/igmm/eddie/GenScotDepression/shen/bakup.dat/ukb_phenotype/'){
    fname=fname_rls$file
    rls.no=fname_rls$rls
    tmp.dat=readRDS(paste0(path,fname))
    ls.field=colnames(tmp.dat)[2:ncol(tmp.dat)] %>%
            strsplit(.,'\\.') %>%
            lapply(function(x) x[2]) %>%
            unlist %>%
            unique
    tmp.block.dic=filter(dat.dic,rls==rls.no)
    keep.field=ls.field[ls.field %in% tmp.block.dic$field_name]
    if (length(keep.field>0)){
        ls.field.output=data.frame(field=keep.field,file=fname,rls=rls.no,stringsAsFactors=F)    
    }else{
        ls.field.output=data.frame(field=NA,file=fname,rls=rls.no,stringsAsFactors=F)
    }

    return(ls.field.output)
}      

fields.available=pblapply(fs.list,extract_available_fields,dat.dic=fields.loose) %>% 
                    bind_rows %>% filter(.,!is.na(field)) 
fields.available$obj.name<- fields.available$file %>%
            strsplit(.,'/') %>%
            lapply(function(x) x[2]) %>%
            unlist
fields.available<-mutate(fields.available,obj.name=paste0(rls,'.',obj.name)) 

# Load data needed ***
files.to.obj=filter(fields.available,!duplicated(file))%>%
        split(., seq(nrow(.)))

for(x in files.to.obj){
  tmp.dat=readRDS(paste0(g.path,x$file))
  tmp.fields=c('f.eid',paste0('f.',fields.available$field,'\\.'))
  tmp.fields=paste0(tmp.fields,collapse='|')
  tmp.dat<-tmp.dat %>%
          select(matches(tmp.fields))
  
  eval(parse(text=paste0(x$obj.name,'=tmp.dat')))
  cat(paste0(x$obj.name,'\n'))
}

# Extract data: function defined above - extract_pheno ***
loose.field.input=fields.available[,c('field','obj.name')]
loose.field.input$field=as.numeric(loose.field.input$field)
colnames(loose.field.input)=c('field.id','dat')
loose.field.input=split(loose.field.input, seq(nrow(loose.field.input)))
loose.fields.toprocess=pblapply(loose.field.input,FUN=extract_pheno,dat.dic=fields.loose)
rm(list=ls(pattern='.rds'))
rm(list=ls(pattern='tmp.'))

for (d in loose.fields.toprocess){
    if(colnames(d$phenotype)[1]==colnames(loose.fields.toprocess[[1]]$phenotype)[1]){
        #dat.loose.field=d$phenotype[,c(2,1)]
        count.n=d$count
        count.n$type=rownames(count.n)
    }else{
        #tmp.pheno=d$phenotype
        #dat.loose.field=merge(dat.loose.field,tmp.pheno,by='f.eid',all.x=T)
        tmp.count.n=d$count
        tmp.count.n$type=rownames(tmp.count.n)
        count.n=merge(count.n,tmp.count.n,by='type',all.x=T)
    }
}

loose.fields.toprocessbind_rows %>% filter(.,!is.na(field)) 

saveRDS(dat.loose.field,file='data/dat.loose_field_noblock.rds')
saveRDS(count.n,file='data/countn.loose.field.rds')

# Combine left and right side ***
fields.loose$Field=as.character(fields.loose$Field)
left.fields=fields.loose$field_name[grep('(left)',fields.loose$Field)]
right.fields<- gsub('(left)','right',x=fields.loose$Field[grep('(left)',fields.loose$Field)]) %>%
                match(.,fields.loose$Field) %>%
                (function(x) fields.loose$field_name[x])
targetdata=dat.loose.field
left.dat=targetdata[,paste0('f.',left.fields)]
right.dat=targetdata[,paste0('f.',right.fields)]


for (i in 1:ncol(left.dat)){
  aveg.source=cbind(left.dat[,i],right.dat[,i])
  tmp.aveg.dat=data.frame(f.eid=targetdata$f.eid,rowMeans(aveg.source,na.rm=T))
  tmp.aveg.dat[rowSums(!is.na(aveg.source))==0,2]=NA
  colnames(tmp.aveg.dat)[2]=paste0(colnames(left.dat)[i],'.lnr')
  if (i==1){
      aveg.dat=tmp.aveg.dat
  }else{
      aveg.dat=merge(aveg.dat,tmp.aveg.dat,by='f.eid',all.x=T)
  }
}

rm.fields=c(right.fields,left.fields)
rm.fields=paste0('f.',rm.fields)

dat.loose.field.lnr=dat.loose.field[,!colnames(dat.loose.field) %in% rm.fields]
dat.loose.field.lnr=merge(dat.loose.field.lnr,aveg.dat,by='f.eid',all.x=T)


# Add back block data (removed from the extract phenotype function) ***
tmp.fields.block.process=colnames(count.n)[count.n[1,]==-999]
tmp.fields.block.process=tmp.fields.block.process[!is.na(tmp.fields.block.process)]

dat.loose.field.lnr=dat.loose.field.lnr[,!colnames(dat.loose.field.lnr) %in% tmp.fields.block.process]

extract_blockpheno <- function(field_dat,dat.dic,instances=0:2){
    field.id=field_dat$field.id
    master.dat=get(field_dat$dat)

     for (i in instances){
           tmp.block.dat<- master.dat %>%
              select(matches(paste0('f.',field.id,'\\.'))) %>%
              select(matches(paste0('\\.',instances,'\\.', collapse="|")))
           tmp.interpolate.dat=rowMeans(tmp.block.dat,na.rm=T)
           tmp.interpolate.dat[rowSums(!is.na(tmp.block.dat))==0]=NA
           if (i==instances[1]){
             tmp.output.dat=data.frame(tmp.interpolate.dat)
             colnames(tmp.output.dat)=paste0('f.',f.ield.id)
           }else{
             loc.to.interpolate=is.na(tmp.output.dat[,names(tmp.output.dat)])
             tmp.output.dat[,names(tmp.output.dat)][loc.to.interpolate]=tmp.interpolate.dat[loc.to.interpolate]             
           }

           # keep a record of N derived from interpolated data
           count.n.instance=sum(!is.na(tmp.interpolate.dat[loc.to.interpolate]))
           if(i==1){count.n.record=sum(!is.na(tmp.interpolate.dat))}else{count.n.record=c(count.n.record,count.n.instance)}
          } 
}

saveRDS(dat.loose.field.lnr,file='data/dat.loose_field_noblock_lnr.rds')

# new data dictionary
rm.fields.dic <- tmp.fields.block.process %>% gsub('f.','',.) %>% as.numeric %>% c(right.fields,.)
fields.loose_noblock_lnr = fields.loose[! fields.loose$field_name %in% rm.fields.dic,]
fields.loose_noblock_lnr$field_tag = as.character(fields.loose_noblock_lnr$field_name)
fields.loose_noblock_lnr$field_tag[fields.loose_noblock_lnr$field_name %in% left.fields]=
      paste0(fields.loose_noblock_lnr$field_tag[fields.loose_noblock_lnr$field_name %in% left.fields],'.lnr')
fields.loose_noblock_lnr$fields_used=fields.loose_noblock_lnr$field_tag
fields.loose_noblock_lnr$fields_used[fields.loose_noblock_lnr$field_name %in% left.fields]=
      paste0(fields.loose_noblock_lnr$fields_used[fields.loose_noblock_lnr$field_name %in% left.fields],
              fields.loose_noblock_lnr$fields_used[fields.loose_noblock_lnr$field_name %in% right.fields])


# Reorganise data based on categories ***
# Add category to the data dictionary
# Remove covariates
#fields.loose_noblock_lnr=fields.loose_noblock_lnr[!fields.loose_noblock_lnr$field_name %in% c(31,34,54,55,21003),]
#write.table(fields.loose_noblock_lnr,file='data/data_dictionary/fields_for_category.csv',quote=F,col.names=T,row.names=F,sep='\t')

fields.category.process=read.delim('data/data_dictionary/fields_for_category.csv',header=T,sep='\t',stringsAsFactors=F)


add_category<-function(kw,col,rplc,cate.dat){
    cate.dat$category[grep(kw,cate.dat[,col])]=rplc
    return(cate.dat)
}
fields.category.process$category=NA
fields.category.process<-add_category('Sociodemographics|Population characteristics|Employment|Recruitment',
                                      'Path','Sociodemographics',cate.dat=fields.category.process) %>%
        add_category('Early life factors|Family history','Path','Early-life factors',.) %>%
        add_category('Lifestyle and environment|Physical activity measurement','Path','Lifestyle measure',.) %>%
        add_category('Physical measures|Health and medical history|Medications|Operations|Medical conditions','Path','Physical measure',.) %>%  
        add_category('Abdominal MRI|Heart MRI|DXA assessment','Path','Body MRI',.) %>%  
        add_category('Assay results > Blood assays','Path','Blood biomarkers',.) %>%  
        add_category('Psychosocial factors','Path','Mental health',.)
        


# Brain imaging   --------------------------------------------------------------------------------------------
# Add data
IDP=readRDS('data/dat.imaging_chunk.rds')
fields.IDP = read.delim('data/data_dictionary/fields.final.brain_imaging_QC_cov_phenotype.txt',header=T,stringsAsFactors=F,)
fields.IDP = fields.IDP[,c(1:18,20,21,19)]

# Cognitive functions   --------------------------------------------------------------------------------------
cognition=readRDS('data/cognition.rds')
fields.cognition=read.delim('data/data_dictionary/fields.final.cognition.txt',header=T,stringsAsFactors=F,)
fields.cognition = fields.cognition[,c(1:18,20,21,19)]
colnames(fields.cognition)[20]='fields_used'

# Mental health questionaires (add online to exisiting loose phenotypes)   -----------------------------------
MHQ=readRDS('/gpfs/igmmfs01/eddie/GenScotDepression/shen/bakup.dat/ukb_phenotype/2020-10-imaging-ukb40531/MHQ.1907.ukb24262.Derived.rds')
MHQ.forprocess <- MHQ %>%
        (function(x) x[,!grepl('^SR|Age.At.MHQ|Gender|No.Info$|^MH\\.|SRConditions|NoComor|AnyComor|^Depressed.Ever$',colnames(x))])
MDD.CIDI=MHQ[,c('f.eid','Depressed.Ever')]
colnames(MDD.CIDI)[2]='MDD.CIDI'
MDD=readRDS('/gpfs/igmmfs01/eddie/GenScotDepression/shen/bakup.dat/ukb_phenotype/mdd_pipeline_MAdams/ukb8238-mdd.mdd_phenotypes.rds')
colnames(MDD)[2:4]=c('MDD.Nerves','MDD.Smith','MDD.ICD')

MDD=merge(MDD,MDD.CIDI,by='f.eid',all.x=T)

MentalHealth=merge(MDD,MHQ,by='f.eid',all=T)

# Data dictionary
fields.MHQ.MDD=data.frame(field_tag=c(colnames(MHQ)[2:ncol(MHQ)],colnames(MDD)[2:ncol(MDD)]),
                          fields_used=NA,
                          category='Mental health')
                          
# Final data             --------------------------------------------------------------------------------------

covs=readRDS('data/covariates.rds')

dat.phewas=merge(dat.loose.field.lnr,IDP,by='f.eid',all.x=T)
dat.phewas=merge(dat.phewas,cognition,by='f.eid',all.x=T)
dat.phewas=merge(dat.phewas,MentalHealth,by='f.eid',all.x=T)

corr.na <- function(corr.na.dat){
      tmp.new.output=corr.na.dat
      tmp.new.output[grep('Prefer not to answer|Do not know',tmp.new.output)]=NA
      return(tmp.new.output)
      }

groot = dat.phewas %>% pbsapply(.,corr.na)

for (i in 2:ncol(dat.phewas)){
      if(!is.numeric(dat.phewas[,i])){dat.phewas[,i]=as.numeric(dat.phewas[,i])}else{}
      if(i%%50==0){cat(paste0(i,'\t'))}else{}
}




dat.phewas=merge(dat.phewas,covs,by='f.eid',all.x=T)

saveRDS(dat.phewas,file='data/dat.phewas.rds')

dat.phewas.imaging = dat.phewas[dat.phewas$f.eid %in% IDP$f.eid,]
saveRDS(dat.phewas.imaging,file='data/dat.phewas_imaging.rds')


# Phenotype col names    ---------------------------------------------------------------------------------------

pheno.totest = rbind(fields.category.process[,c('field_tag','fields_used','category')],
                     fields.IDP[fields.IDP$category!='Brain imaging QC and covariates',c('field_tag','fields_used','category')],
                     fields.cognition[,c('field_tag','fields_used','category')],
                     fields.MHQ.MDD[,c('field_tag','fields_used','category')])
saveRDS(pheno.totest,file='pheno.totest.rds')