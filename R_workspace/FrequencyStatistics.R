library(gdata)
library(ggplot2)

fusion_file <- '~/Roche/Data/filteredFusions_v2.xls'
suppressWarnings(fusion_data <- read.xls(fusion_file))
freq_mat <- with(fusion_data,tapply(Percent,list(factor(Tumor),factor(gene3p)),sum))
freq_table <- NULL

for(i in 1:nrow(freq_mat)){
	for(j in 1:ncol(freq_mat)){
		if(!is.na(freq_mat[i,j])){
			freq_table <- rbind(freq_table,c(rownames(freq_mat)[i],colnames(freq_mat)[j],freq_mat[i,j]))
		}
	}
}

colnames(freq_table) <- c('Tumor','Gene','Percent')
freq_table <- as.data.frame(freq_table)
freq_table$Percent <- as.numeric(as.character(freq_table$Percent))
freq_table <- freq_table[order(freq_table$Percent,decreasing = T),]

freq_table$Gene2 <- factor(freq_table$Gene,level=unique(freq_table$Gene))


# Plot --------------------------------------------------------------------
p <- ggplot(data=freq_table,aes(Tumor,Percent,fill=Gene)) + 
  geom_bar(stat='identity',position='dodge')
png('~/Roche/20140731/Figure 1.png',width=1600,height=600,res=100)
p + theme(axis.text=element_text(size=14),axis.title=element_text(size=22))
dev.off()

freq_table$xVal <- seq(from=1,to=nrow(freq_table))
p <- ggplot(data=freq_table,aes(xVal,Percent)) + 
  geom_bar(stat='identity',width=0.8)
tick_label_x <- apply(cbind(as.character(freq_table$Tumor),as.character(freq_table$Gene)),
                      1,function(x) paste(x[1],x[2],sep=':'))
png('~/Roche/20140731/Figure 2.png',width=2400,height=600,res=100)
p + theme(axis.text.x=element_text(size=14,angle=90,vjust=0.5,hjust=1),
          axis.text.y=element_text(size=14),axis.title=element_text(size=22)) +
  scale_x_continuous(name='',breaks=seq(1,nrow(freq_table)),labels=tick_label_x)
dev.off()

png('~/Roche/20140804/Figure 1.png',width=1600,height=600,res=100)
p <- ggplot(data=freq_table,aes(Tumor,Percent,fill=Gene)) + 
  geom_bar(stat='identity')
p + theme(axis.text=element_text(size=14),axis.title=element_text(size=22)) +
  scale_x_discrete(limits=c('GBM','LGG','HNSC','LUAD','LUSC','ESCA','STAD','COAD',
                            'LIHC','READ','KIRC','KIRP','BLCA','PRAD','BRCA','CESC',
                            'UCEC','OV','PCPG','PAAD','THCA','SKCM','SARC','LAML'))
dev.off()