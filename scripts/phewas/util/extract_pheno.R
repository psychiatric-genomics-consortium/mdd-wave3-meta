library(dplyr)

extract_pheno <- function(field_dat,dat.dic,instances=0:2){
  field.id=field_dat$field.id
  master.dat=get(field_dat$dat)
  tmp.dic=filter(dat.dic,FieldID==field.id)
  # get data **
  tmp.block.dat <- master.dat %>%
    select(matches(paste0('f.',field.id,'\\.'))) %>%
    select(matches(paste0('\\.',instances,'\\.', collapse="|")))
  
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
  
  tmp.dat$f.eid=master.dat$f.eid
  output.pheno=list('phenotype'=tmp.dat,'count'=count.n.record)
  return(output.pheno)  
}