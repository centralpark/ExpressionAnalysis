## Add two new columns to table about clustering

# args[1]: input all15.txt file
# args[2]: output file name
# args[3]: temporary folder that stores computation results returned by Roche cluster
# Output: write a new file with additional columns

args <- commandArgs(trailingOnly = T)
inputFile <- args[1]
outputFileName <- args[2]
folder <- args[3]

source('/home/staff/hes14/Roche/R_workspace/UserFunctions.R')

datIn <- read.delim(inputFile, stringsAsFactors=F)
nRow <- nrow(datIn)
datIn$pValue <- NULL
datIn$Estimate <- NULL

if(substr(folder,nchar(folder),nchar(folder))=='/') folder <- substr(folder,1,nchar(folder)-1)

for(i in 1:nRow){
  if(datIn[i,'Fusions']<3) next
  if(datIn[i,'IsReal']!=1) next
  cancer <- datIn[i,'Tumor']
  gene5p <- datIn[i,'gene5p']
  gene3p <- datIn[i,'gene3p']
  data_fname <- file.path(folder,paste(cancer,gene5p,gene3p,'RData',sep='_'))
  if(!file.exists(data_fname)){
    msg <- sprintf('Consider run gene differential analysis for fusion %s-%s in %s',gene5p,gene3p,
                   cancer)
    message(msg)
    next
  }
  value <- computeCluster(data_fname, cancer, gene5p, gene3p)
  if(is.null(value)) next
  datIn[i,'pValue'] <- value$p.value
  datIn[i,'Estimate'] <- -diff(value$estimate)
}

write.table(datIn,file=outputFileName,quote=F,sep='\t',row.names=F)