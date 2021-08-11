run_model <- function(ls.mod,x.dat_short,x.dat_long){
   # define vars
   dep = as.character(ls.mod[1])
   factor = as.character(ls.mod[2])
   covs = as.character(ls.mod[3])
   mod.type = as.character(ls.mod[4])
   
   # run model
   if (mod.type=='lme'){
      # model.expression
      fh_r=1:(nrow(x.dat_long)/2)
      sh_r=(nrow(x.dat_long)/2+1):nrow(x.dat_long)
      
      x.dat_long[fh_r,dep]=scale(x.dat_long[fh_r,dep])
      x.dat_long[sh_r,dep]=scale(x.dat_long[sh_r,dep])
      if(is.numeric(x.dat_long[,factor])){
         x.dat_long[fh_r,factor]=scale(x.dat_long[fh_r,factor])
         x.dat_long[sh_r,factor]=scale(x.dat_long[sh_r,factor])
      }
      
      mod=paste0(dep,'~',covs,'+',factor)
      
      fit=lme(as.formula(as.character(mod)),data=x.dat_long,na.action=na.exclude,random=~1|f.eid,control=lmeControl(opt = "optim"))
      table = summary(fit)$tTable
      tarv = nrow(table)
      stats = table[tarv,c(1,2,4,5)]
      mod_result = data.frame(mod_name=paste0(dep,'~',factor),t(stats))
      colnames(mod_result)[2:5]=c('beta','std','t.value','p.value')
      
   }else{
      dep.dat=x.dat_short[,dep]            
      if (length(table(dep.dat))==2){
         mod=paste0('as.factor(',dep,')~',covs,'+scale(',factor,')')
         fit=glm(as.formula(as.character(mod)),data=x.dat_short,na.action=na.exclude,family = 'binomial')
      }else{
         x.dat_short[,dep]=as.numeric(x.dat_short[,dep])
         mod=paste0('scale(',dep,')~',covs,'+scale(',factor,')')
         fit=glm(as.formula(as.character(mod)),data=x.dat_short,na.action=na.exclude)
      }            
      
      table = summary(fit)$coefficients
      tarv = nrow(table)
      stats = table[tarv,c(1:4)]
      mod_result = data.frame(mod_name=paste0(dep,'~',factor),t(stats))
      colnames(mod_result)[2:5]=c('beta','std','t.value','p.value')
   }
   
   return(mod_result)
}

reg_phewasStyle <- function (ls.models,dat_short,dat_long=NA,correctByFactor=F){

      tmp.result = pbapply(X = ls.models,MARGIN = 1,FUN = run_model,x.dat_short=dat_short,x.dat_long=dat_long)
      result.table =  matrix(unlist(tmp.result),ncol=5,byrow = T)
      result.table = data.frame(result.table)
      colnames(result.table) = c('mod_name','beta','std','t.value','p.value')
      result.table$mod_name = lapply(tmp.result, function(l) as.character(l[[1]]))
      result.table = data.frame(ls.models[,1:2],result.table,stringsAsFactors = F)
      
      ls.factor=unique(ls.models$p_batch)
      result.table$p.corrected=99999
      if (correctByFactor==T){
            for (f in ls.factor){
                  loc=grep(f,ls.models$p_batch)
                  result.table$p.corrected[loc]=p.adjust(result.table$p.value[loc],method='fdr')
            }
      }else{
            result.table$p.corrected=p.adjust(result.table$p.value,method='fdr')
            
      }
      
      return(result.table)
}