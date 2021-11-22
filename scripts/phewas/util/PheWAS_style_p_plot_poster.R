# By shen
# 14/05/2020
# Not very well annotated. But don't ask me to annotate... I'll be depressed
# *TargetResult: result table (column names:beta, std, p.value, p.corrected)
# *color_theme_usr:Shen/a list of colours
# *shape_sig: T=a different shape for significant assoc
# *category.input: a table of categories, first column=keywords for dependent var
#  second column=category label
# *labels_annot: a table of labels, 1st col=phenotype name in the result table
#  2nd col=phenotype label
# *outputpath: path for a tiff image to be generated.
# *plot_title: plot title
# *add_category_name: T=add category labels
# *fig_size: default=c(24,20)

library(dplyr)
library(ggplot2)
library(magrittr)
library(ggpubr)
library(cowplot)
library(ggrepel)
############################## XShen style p plot ##############################

p_plot<-function(TargetResult,color_theme_usr='Shen',shape_sig=T,sig_nominal=F,
                 labels_annot,outputpath=NA,plot_title=NA,
                 add_category_name=F,fig_size=c(24,20),y_lim=40){
  
  # Add label for some regions to have a different shape
  TargetResult$p.value[TargetResult$p.value==0]=1e-308
  TargetResult$shape=1
  if (shape_sig==T){
    TargetResult$shape[TargetResult$p.corrected<0.05]=19
  }else{}
  
  TargetResult$sig = 99999
  TargetResult$sig[TargetResult$p.corrected<0.05]='*'
  TargetResult$sig[TargetResult$sig==99999]=''
  
  # Set color theme
  if(length(color_theme_usr)>1){
    cl.theme=color_theme_usr
  }else if (color_theme_usr=='Shen'){
    cl.theme=c('orangered1','orange1','royalblue3','lightseagreen','maroon',
               'salmon2','palevioletred2','olivedrab','darkslategray3',
               'slategray3','mediumpurple1','orangered1','orange1','royalblue3','lightseagreen','maroon',
               'salmon2','palevioletred2','olivedrab','darkslategray3',
               'slategray3','mediumpurple1')
  }else{}
  
  # Identify significance line
  if (sum(TargetResult$p.corrected<0.05)>0){
    max.p.corrected = max(TargetResult$p.value[TargetResult$p.corrected<0.05])
    pt=TargetResult$p.value[TargetResult$p.value==max.p.corrected]
    pt=-log10(max(pt))
    pt.color='pink'
  }else{
    pt=-log10(0.05)
    pt.color='grey'
  }
  
  if (sig_nominal==T){
    pt=-log10(0.05)
    pt.color='grey'
  }

    
  # Add labels
  if (sum(!is.na(labels_annot))>0){
    TargetResult$labels=rep(99999,nrow(TargetResult))
    for (i in 1:nrow(labels_annot)){
      TargetResult$labels[TargetResult$dependent==labels_annot[i,1]]=labels_annot[i,2]
    }
    TargetResult$labels[TargetResult$labels==99999]=''
  }else{
    TargetResult$labels=''
  }

  
  # Add a column to fix the order of phenotypes
  TargetResult$ord=1:nrow(TargetResult)
  
  # Set x axis to print categories
  axis.set <- TargetResult %>% 
    group_by(category) %>% 
    summarize(center = (max(ord) + min(ord)) / 2-8) %>% 
    as.data.frame
  
  # Make the plot
  fig=ggplot(TargetResult, aes(x=ord, y=-log10(p.value),label=labels)) +
    geom_point(size=1, shape=TargetResult$shape,stroke=1.5,
               aes(colour = category),
               alpha=0.5)+
    scale_colour_manual(values = cl.theme)+
    geom_text_repel(box.padding = unit(0.5,'lines'),segment.size = 0.2,
                    max.iter = 10000, max.overlaps = 500,
                    nudge_x = 0.8, nudge_y = 0.8,size = 5.4)+
    # geom_label( 
    #   data=TargetResult %>% filter(labels!=''), # Filter data first
    #   aes(label=labels),
    #   nudge_x = 1, nudge_y = 1, 
    #   check_overlap = T
    # )+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_text(angle = 90, size = 13,hjust=1, vjust = 1, face='bold'),
          axis.ticks.x=element_blank(),
          legend.position = 'none',
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "grey"),
          axis.line.x = element_blank(),
          plot.title=element_text(lineheight = 1,face='bold'),
          plot.margin = unit(c(1,3,1,1), "lines"))+
    geom_hline(yintercept=0 , color = "grey", size=0.5)+
    geom_hline(yintercept=pt, linetype="dashed", color = pt.color, size=1)+
    geom_hline(yintercept=308, linetype="dashed", color = 'light grey', size=0.5)+
    annotate(
      "text", label = paste0("p < 10^308"),
      x = max(TargetResult$ord)-10, y = 308, size = 7, colour = "dark grey",alpha = 0
    )+
    ylab('-log10(p)')+
    ylim(c(0,y_lim))
  if (!is.na(plot_title)){
    fig=fig+ggtitle(plot_title)
  }
  if (add_category_name==T){
    fig=fig+
      scale_x_continuous(label = axis.set$category, breaks = axis.set$center)
  }else{
    fig=fig+theme(axis.text.x=element_blank())
  }

  # Print figure to a file
  if (!is.na(outputpath)){
    ggsave(plot = fig,filename = outputpath,
           width = fig_size[1], height = fig_size[2], 
           units = 'cm', dpi = 300)
  }else{}
  
  return(fig)
}

