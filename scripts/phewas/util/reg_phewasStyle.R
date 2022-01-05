run_model <- function(ls.mod,mod.dat_short,mod.dat_long){
   # define vars
   mod.dep = as.character(ls.mod[1])
   mod.factor = as.character(ls.mod[2])
   mod.covs = as.character(ls.mod[3])
   mod.type = as.character(ls.mod[4])
   
   # run model
   if (mod.type=='lme'){
      # model.expression
      fh_r=1:(nrow(mod.dat_long)/2)
      sh_r=(nrow(mod.dat_long)/2+1):nrow(mod.dat_long)
      
      mod.dat_long[fh_r,mod.dep]=scale(mod.dat_long[fh_r,mod.dep])
      mod.dat_long[sh_r,mod.dep]=scale(mod.dat_long[sh_r,mod.dep])
      if(is.numeric(mod.dat_long[,mod.factor])){
         mod.dat_long[fh_r,mod.factor]=scale(mod.dat_long[fh_r,mod.factor])
         mod.dat_long[sh_r,mod.factor]=scale(mod.dat_long[sh_r,mod.factor])
      }
      
      mod=paste0(mod.dep,'~',mod.covs,'+',mod.factor)
      
      fit=lme(as.formula(as.character(mod)),data=mod.dat_long,
              na.action=na.exclude,random=~1|f.eid,control=lmeControl(opt = "optim"))
      tmp.ci = intervals(fit,which='fixed')$fixed %>% as.data.frame %>% 
         select(Lower_95CI=lower,Upper_95CI=upper) %>% 
         tail(1)
      tmp.res = summary(fit)$tTable %>% 
         as.data.frame %>% 
         select(beta=Value,std=Std.Error,t.value=`t-value`,p.value=`p-value`) %>% 
         tail(1) %>% 
         mutate(mod_name = paste0(mod.dep,'~',mod.factor),tmp.ci) %>% 
         select(mod_name, everything())
      
      
   }else{
      dep.dat=mod.dat_short[,mod.dep]            
      if (length(table(dep.dat))==2){
         mod=paste0(mod.dep,'~',mod.covs,'+scale(',mod.factor,')')
         fit=glm(as.formula(as.character(mod)),data=mod.dat_short,na.action=na.exclude,family = 'binomial')
      }else{
         mod=paste0('scale(',mod.dep,')~',mod.covs,'+scale(',mod.factor,')')
         fit=glm(as.formula(as.character(mod)),data=mod.dat_short,na.action=na.exclude)
      }            
      tmp.ci = confint(fit) %>% as.data.frame %>% 
         select(Lower_95CI=`2.5 %`,Upper_95CI=`97.5 %`) %>% 
         tail(1)
      tmp.res = summary(fit)$coefficients %>% 
         as.data.frame %>% 
         select(beta=Estimate,std=`Std. Error`,t.value=`t value`,p.value=`Pr(>|t|)`) %>% 
         tail(1)%>% 
         mutate(mod_name = paste0(mod.dep,'~',mod.factor),tmp.ci) %>% 
         select(mod_name, everything())
   }
   
   return(tmp.res)
}


reg_phewasStyle <- function (ls.models,dat_short,dat_long,correctByFactor=F){
      
      result.table = ls.models %>% split(.,seq(nrow(.))) %>% 
         pblapply(.,FUN = run_model,
                           mod.dat_short=dat_short,mod.dat_long=dat_long) %>% 
         bind_rows %>% 
         as.data.frame %>% 
         mutate(ls.models[,1:2])
      
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
      rownames(result.table)=NULL
      return(result.table)
}