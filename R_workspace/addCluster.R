# args[1]: input all*.txt file
# args[2]: output file name
# args[3]: temporary folder to store computation results returned by Roche cluster
# args[4]: Roche user name, eg. hes14

args <- commandArgs(trailingOnly = T)
inputFile <- args[1]
outputFileName <- args[2]
clustFolder <- args[3]
user <- args[4]

source('/home/staff/hes14/Roche/R_workspace/UserFunctions.R')

# maximum number of nodes to use, do not occupy too much resources.
maxNodeNum <- 30
datIn <- read.delim(inputFile, stringsAsFactors=F)
nRow <- nrow(datIn)

if(substr(clustFolder,nchar(clustFolder),nchar(clustFolder))=='/') clustFolder <- substr(clustFolder,1,nchar(clustFolder)-1)

## To compute only for selected rows, run this block
# selectRow <- c(315427, 326897)
# for(i in selectRow){
#   # leave it if number of cases is too low
#   if(datIn[i,'Fusions']<3) next
#   cmd <- paste('qstat -u',user)
#   qstat <- system(cmd,intern=TRUE)
#   while(length(qstat)>(maxNodeNum+2)){
#     Sys.sleep(600)
#     qstat <- system(cmd,intern=TRUE)
#   }
#   run.cluster(datIn[i,'Tumor'],datIn[i,'gene5p'],datIn[i,'gene3p'],folder=clustFolder)
# }
## end block

for(i in 1:nRow){
  if(datIn[i,'IsReal']!=1) next
  # leave it if number of cases is too low
  if(datIn[i,'Fusions']<3) next
  cmd <- paste('qstat -u',user)
  qstat <- system(cmd,intern=TRUE)
  while(length(qstat)>(maxNodeNum+2)){
    Sys.sleep(600)
    qstat <- system(cmd,intern=TRUE)
  }

  fname <- file.path(clustFolder,paste(datIn[i,'Tumor'],datIn[i,'gene5p'],datIn[i,'gene3p'],'RData',sep='_'))
  if(file.exists(fname)) next

  run.cluster(datIn[i,'Tumor'],datIn[i,'gene5p'],datIn[i,'gene3p'],folder=clustFolder)
}

