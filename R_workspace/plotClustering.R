library(ggplot2)
source('/Users/HSH/Roche/R_workspace/UserFunctions.R')

folder <- '/home/staff/hes14/temp'
data_files <- list.files(path=folder,pattern = '*_RData',full.names = T)
for(df in data_files){
  fname_split <- strsplit(basename(df),'_')[[1]]
  cancer <- fname_split[1]
  gene5p <- fname_split[2]
  gene3p <- fname_split[3]
  cTable <- clusterTable(df, cancer, gene5p, gene3p)
  if(class(cTable)!='data.frame') next
  fname <- paste(folder,'/',
                 paste(cancer,gene5p,gene3p,'clustering',sep='_'),
                 sep='')
  write.table(cTable,file=fname,quote=F,sep='\t',row.names=F)
}

plotDF <- NULL
files <- list.files(path=folder,pattern = '*_clustering',full.names = T)
for(f in files){
  dat <- read.delim(f,stringsAsFactors=F)
  nSample <- nrow(dat)
  fname_split <- strsplit(basename(f),'_')[[1]]
  plotDF_new <- data.frame(fusion=rep(paste(fname_split[2],fname_split[3],fname_split[1],sep='-'),nSample), num=rep(1,nSample), 
                       harbor=dat$HarborFusion, stringsAsFactors=F)
  plotDF <- rbind(plotDF,plotDF_new)
}


png('/Users/HSH/Roche/20140912/Figure 1.png',width=12*300,height=7*300,res=300)
ggplot(data=plotDF, aes(x=fusion,y=num,fill=harbor)) + geom_bar(stat='identity') +
  scale_fill_manual(values=c('#D9D9D9','#3795ED'),name='',labels=c('Without fusion','With fusion')) +
  scale_x_discrete(limits=c('FGFR3-TACC3-BLCA','FGFR3-TACC3-GBM','EGFR-SEPT14-GBM','MLL-ELL-LAML','PML-RARA-LAML','RUNX1-RUNX1T1-LAML','DNAJB1-PRKACA-LIHC','EML4-ALK-LUAD','TMPRSS2-ERG-PRAD','TMPRSS2-ETV4-PRAD','CCDC6-RET-THCA','PAX8-PPARG-THCA'), 
                   labels=c('FGFR3-TACC3\nBLCA','FGFR3-TACC3\nGBM','EGFR-SEPT14\nGBM','MLL-ELL\nLAML','PML-RARA\nLAML','RUNX1-RUNX1T1\nLAML','DNAJB1-PRKACA\nLIHC','EML4-ALK\nLUAD','TMPRSS2-ERG\nPRAD','TMPRSS2-ETV4\nPRAD','CCDC6-RET\nTHCA','PAX8-PPARG\nTHCA'),
                   name='Fusion') +
  scale_y_continuous(name='Number of patients') +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=20), axis.text.x = element_text(size=16),
        axis.title.y = element_text(vjust=1.5,size=20), axis.text.y = element_text(size=16),
        legend.text = element_text(size=16))
dev.off()