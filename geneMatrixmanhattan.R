library(ggplot2)
library(ggrepel)

data <- read.table("mdd_fastbat_geneMatrix.gene.fastbat",header=TRUE,sep="")

data$mid <- (data$Start + data$End) / 2
data$logp <- -log10(data$Pvalue)

chrNum <- 22
for (i in 1:22){ ndx <- which(data$Chr==i)
                 lstMrk <- max(data[ndx, "mid"])
                 if (i < chrNum) ndx2 <- which(data[,"Chr"]==i+1)
                 if (i < chrNum) data[ndx2, "mid"] <- data[ndx2, "mid"] + lstMrk
}
bpMidVec <- vector(length=chrNum)
for (i in 1:chrNum){ndx <- which(data[,"Chr"]==i)
                    posSub <- data[ndx,"mid"]
                    bpMidVec[i] <- ((max(posSub) - min(posSub))/2) + min(posSub)
}

fastbat <- cbind(as.numeric(data$Chr),data$mid,data$logp,data$Gene)
fastbat <- data.frame(fastbat)
colnames(fastbat) <- c("chr","pos","logp","Gene")
fastbat$chr<-as.numeric(fastbat$chr)
fastbat$pos<-as.numeric(fastbat$pos)
fastbat$logp<-as.numeric(fastbat$logp)

output<-fastbat
output$Gene[which(output$logp < 16)]<-""
bonf<-paste0("P = ",as.character(signif(0.05/nrow(output),digits=3)))

p <- ggplot(output,aes(x=pos, y=logp, label = Gene, colour=as.factor(chr)), alpha=1/2)  + theme(panel.background = element_blank()) + theme(axis.line = element_line(colour = "black"))
p2 <- p + geom_point(size=2,shape=20) + scale_color_manual(values=rep(c('steelblue4', 'darkorange3'), 11)) + scale_x_continuous(labels=as.character(1:chrNum), breaks=bpMidVec) + theme(legend.position='none')
p3 <- p2 + theme(axis.title = element_text(size=14)) + theme(axis.text = element_text(size=10,colour="black")) # + ggtitle(bquote(atop(.(plot.title), atop(italic(.(plot.subtitle)), "" )))) + theme(plot.title = element_text(size=25))
p4 <- p3 + labs(x="Chromosome",y="-log10(P)") + coord_cartesian(ylim = c(-0.2,50),xlim = c(min(data$mid),max(data$mid)+300000000)) + scale_y_continuous(breaks=seq(0,50,5)) + theme(panel.grid.major.x = element_blank()) + theme(panel.grid.minor.x = element_line(colour="gray45", size=0.2)) + theme(panel.grid.major.y = element_line(colour="gray45", size=0.2))
p5 <- p4 + geom_text_repel(max.overlaps = Inf, size = 2.5, nudge_x = 1, nudge_y = 1)
p6 <- p5 + geom_segment(aes(x = min(data$mid), xend = max(data$mid)+50000000), y=-log10(0.05/nrow(output)), yend=-log10(0.05/nrow(output)), linetype=2, col='black', lwd=1) + annotate("text",x=max(data$mid)+240000000,y=-log10(0.05/nrow(output)),label=bonf,size=4)
ggsave("geneMatrixmanhattan.tiff",p6,width=10,height=6.8,units="in")
