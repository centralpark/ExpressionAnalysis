args <- commandArgs(trailingOnly = T)
folder <- args[1]
folder2 <- args[2]
out_fname_age <- args[3]
out_fname_os <- args[4]


library(snowfall)

# number of cores
nCPU = 12

compare <- function(fname,response=NULL){
  dat <- read.table(fname,sep='\t',header=T)
  population_1 <- dat[dat[,3]=='Y',2]
  population_2 <- dat[dat[,3]=='N',2]
  switch(response,
         age={result <- wilcox.test(age~harbor_fusion,data=dat,conf.int=T)},
         OS = {result <- wilcox.test(OS~harbor_fusion,data=dat,conf.int=T)},
         stop('Either age or OS is accepted response.')
         )
  name <- basename(fname)
  name_split <- strsplit(name,'_')
  result$Y <- length(population_1)
  result$N <- length(population_2)
  result$cancer <- name_split[[1]][1]
  result$gene5p <- name_split[[1]][2]
  result$gene3p <- name_split[[1]][3]
  result
}

fileList <- dir(folder,full.names=T)
fileList <- fileList[!file.info(fileList)$isdir]

fileList2 <- dir(folder2,full.names=T)
fileList2 <- fileList2[!file.info(fileList2)$isdir]

sfInit(parallel=TRUE,cpus=nCPU)
result <- sfLapply(fileList,compare,response='age')
resultOS <- sfLapply(fileList2,compare,response='OS')
sfStop()

ofile <- file(out_fname_age,'w')
line <- paste('Cancer','Gene5p','Gene3p','SampleSizeWithFusion',
              'SampleSizeWithoutFusion','W','pValue','Shift',sep='\t')
writeLines(line,con=ofile)
for(i in 1:length(result)){
  result_this <- result[[i]]
  line <- with(result_this,paste(cancer,gene5p,gene3p,Y,N,statistic,p.value,estimate,
                                 sep='\t'))
  writeLines(line,con=ofile)
}
close(ofile)

ofile <- file(out_fname_os,'w')
line <- paste('Cancer','Gene5p','Gene3p','SampleSizeWithFusion',
              'SampleSizeWithoutFusion','W','pValue','Shift',sep='\t')
writeLines(line,con=ofile)
for(i in 1:length(resultOS)){
  result_this <- resultOS[[i]]
  line <- with(result_this,paste(cancer,gene5p,gene3p,Y,N,statistic,p.value,estimate,
                                 sep='\t'))
  writeLines(line,con=ofile)
}
close(ofile)
