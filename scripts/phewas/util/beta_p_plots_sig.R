beta_fig <- function(tmp.dat.fig,target.cat=NA,fig.output='results/phewas/figs/',is.allres=F){
  if(is.allres==T){
    cat.dat = tmp.dat.fig
  }else{
    cat.dat = tmp.dat.fig %>% filter(lv1==target.cat)
  }
  cat.dat = cat.dat %>% 
    mutate(sig.shape=ifelse(p.fdr<0.05,19,1)) %>% 
    .[order(.$b,decreasing=F),] %>% 
    mutate(ord=1:nrow(.))
  
  fig.tmp=
  ggplot(cat.dat, aes(x=reorder(Field,ord), y=b,color=method)) + 
  geom_point(position=position_dodge(width = 0.03), stat="identity") +
  geom_errorbar(aes(ymin=b-se, ymax=b+se), width=0.05,
                position=position_dodge(width = 0.02),colour="grey")+
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.x = element_line(size=0.5),
    axis.text=element_text(size=10), axis.title=element_text(size=11),
    plot.title = element_text(lineheight=1, face="bold", vjust=1, hjust=0.5,size=9),
    strip.text = element_text(size=8),
    plot.margin=unit(c(1,1,1,3),'mm'),
    legend.position = "bottom") +
  ylab("Standardised effect size") + xlab("\n\n") +
  #scale_y_continuous(limits=c(-3,3))+
  #scale_y_reverse()+
  scale_x_discrete(position='top')+
  geom_hline(yintercept=0,color = "black", size=0.3)+
  coord_flip()
  
  if(!is.na(target.cat)){
    fig.tmp=fig.tmp+
      facet_grid(rows = vars(category),scales = "free_y",space = "free_y")+
      ggtitle(target.cat)
  }
  
  if(!is.na(fig.output)&(sum(nchar(cat.dat$Field)>50)>5)){
    fig.output=fig.output %>% paste0(.,'/fig.compare.',target.cat %>% gsub(' ','',.) %>% gsub('(','',.,fixed = T) %>% gsub(')','',.),'.png')
      ggsave(plot = fig.tmp,filename = here::here(fig.output),
             device='png',units='cm',width=35,height=2+0.4*nrow(cat.dat))
  }else if(!is.na(fig.output)){
    fig.output=fig.output %>% paste0(.,'/fig.compare.',target.cat %>% gsub(' ','',.) %>% gsub('(','',.,fixed = T) %>% gsub(')','',.),'.png')
    ggsave(plot = fig.tmp,filename = here::here(fig.output),
             device='png',units='cm',width=35,height=2+0.26*nrow(cat.dat))
  }
  
  return(fig.tmp)
}


p_fig <- function(tmp.dat.fig,target.cat,fig.output='results/phewas/figs/'){
  cat.dat = tmp.dat.fig %>% filter(lv1==target.cat) %>% 
    mutate(shape=ifelse(p.fdr<0.05,19,1))
  
  fig.tmp=
  ggplot(cat.dat, aes(x=reorder(Field,ord), y=-log10(pval),color=method)) + 
  geom_point(position=position_dodge(), stat="identity") +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line.x = element_line(size=0.5),
    axis.text=element_text(size=10), axis.title=element_text(size=11),
    plot.title = element_text(lineheight=1, face="bold", vjust=1, hjust=0.5,size=9),
    strip.text = element_text(size=8),
    plot.margin=unit(c(1,1,1,3),'mm')) +
  ylab("Standardised effect size") + xlab("\n\n") +
  #scale_y_continuous(limits=c(0,3))+
  #scale_y_reverse()+
  scale_x_discrete(position='top')+
  geom_hline(yintercept=0,color = "black", size=0.3)+
  ggtitle(target.cat)+
  facet_grid(rows = vars(category),scales = "free_y",space = "free_y")+
  coord_flip()
  
  if(!is.na(fig.output)&(sum(nchar(cat.dat$Field)>50)>5)){
    fig.output=fig.output %>% paste0(.,'/fig.compare.',target.cat %>% gsub(' ','',.) %>% gsub('(','',.,fixed = T) %>% gsub(')','',.),'.png')
      ggsave(plot = fig.tmp,filename = here::here(fig.output),
             device='png',units='cm',width=28,height=2+0.4*nrow(cat.dat))
  }else if(!is.na(fig.output)){
    fig.output=fig.output %>% paste0(.,'/fig.compare.',target.cat %>% gsub(' ','',.) %>% gsub('(','',.,fixed = T) %>% gsub(')','',.),'.png')
    ggsave(plot = fig.tmp,filename = here::here(fig.output),
             device='png',units='cm',width=28,height=2+0.26*nrow(cat.dat))
  }
  
  return(fig.tmp)
}