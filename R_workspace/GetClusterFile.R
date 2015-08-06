args <- commandArgs(trailingOnly = T)
dataFile <- args[1]
cancer <- args[2]
gene5p <- args[3]
gene3p <- args[4]
output_fname <- args[5]

source('/Users/HSH/Roche/R_workspace/UserFunctions.R')
result <- clusterTable(dataFile, cancer, gene5p, gene3p)
write.table(result,file=output_fname,quote=F,sep='\t',row.names = F)
