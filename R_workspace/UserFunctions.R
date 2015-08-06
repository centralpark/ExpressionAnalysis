## Assign color to the labels of cluster dendrogram
colLab <- function(n,map=NULL){
  if(missing(map)) stop('Argument map must be provided.')
  if(is.leaf(n)){
    a <- attributes(n)
    fType <- map[map$X==a$label,'harbor_fusion']
    switch(fType,
           Y = (attr(n,'nodePar') <- c(a$nodePar, lab.col='red')),
           N = (attr(n,'nodePar') <- c(a$nodePar, lab.col='green')) )
  }
  n
}

## Get patient ID from copied GeneHub data
get.pid <- function(file,gene3p,diagnosis,gene5p = NULL,ordered=FALSE,colDiagnosis=5){
  if(missing(file) | missing(gene3p) | missing(diagnosis)) stop('Provide required arguments!')
  dataGH <- read.table(file,sep='\t',stringsAsFactors=F)
  dataGH$harbor_fusion <- NA
  if(!missing(gene5p)){
    select <- (dataGH[,1]==gene5p & dataGH[,2]==gene3p & 
                 (grepl(diagnosis,dataGH[,colDiagnosis])))
    dataGH[select,'harbor_fusion'] <- 'Y'
    if(ordered){
      select <- (dataGH[,1]!=gene5p & dataGH[,2]==gene3p & (grepl(diagnosis,dataGH[,colDiagnosis])))
    } else{
      select <- (dataGH[,1]==gene3p | dataGH[,2]==gene3p & 
                   dataGH[,1]!=gene5p & (grepl(diagnosis,dataGH[,colDiagnosis])))
    }
    dataGH[select,'harbor_fusion'] <- 'O'
  }
  select <- !is.na(dataGH$harbor_fusion)
  resultData <- dataGH[select,c(1,2,3,ncol(dataGH))]
  resultData[,3] <- substr(resultData[,3],1,12)
  resultData
}


prepare.data <- function(cancer,gene5p,gene3p,folder=NULL){
  if(missing(folder)) folder <- getwd()
  cmd <- paste('/apps/Python-2.7.2/bin/python /Users/HSH/Roche/workspace/GeneFusion/DESeq/preparedata.py',cancer,gene5p,gene3p,folder)
  system(cmd)
}

run.cluster <- function(cancer,gene5p,gene3p,folder=NULL){
  if(missing(folder)) folder <- getwd()
  if(substr(folder,nchar(folder),nchar(folder))=='/') folder <- substr(folder,1,nchar(folder)-1)
  prepare.data(cancer,gene5p,gene3p,folder)
  gct_fname <- file.path(folder,paste(cancer,'.gct',sep=''))
  # avoid possible read/write conflict later
  if(file.exists(gct_fname)){
    gct_fname_new <- gct_fname
    while(file.exists(gct_fname_new)) gct_fname_new <- file.path(folder,paste(cancer,sample(1:1000,1),'.gct',sep=''))
    file.copy(gct_fname, gct_fname_new)
    gct_fname <- gct_fname_new
  }
  col_fname <- file.path(folder,paste(cancer,gene5p,gene3p,'colData',sep='_'))
  data_fname <- file.path(folder,paste(cancer,gene5p,gene3p,'RData',sep='_'))
  result_fname <- file.path(folder,paste(cancer,gene5p,gene3p,'result',sep='_'))
  pbs_str <- sprintf('#!/bin/sh

#PBS -N %s
#PBS -m e
#PBS -M siheng.he@roche.com
#PBS -V
#PBS -q short
#PBS -j oe
#PBS -l walltime=100:00:00
#PBS -l select=1:ncpus=12:mem=12gb
                 
/apps64/R-3.1.1/bin/R --vanilla --args %s %s %s %s < /Users/HSH/Roche/R_workspace/DEanalysis.R',
                     paste(gene5p,gene3p,sep='_'),gct_fname,col_fname,data_fname,result_fname)
  pbs_file <- tempfile(fileext='.pbs')
  cat(pbs_str,file=pbs_file)
  cmd <- paste('qsub',pbs_file)
  system(cmd)
  unlink(pbs_file)
}


## Calculate p-value and shift of clustering of samples
computeCluster <- function(dataFile, cancer, gene5p, gene3p, all8File='/Users/HSH/Roche/Data/all8.txt'){
  require(gdata)
  require(DESeq2)
  load(dataFile)
  dataCol <- data.frame(dds@colData@listData)
  cols <- c(1,2,3,5)
  if(class(dataCol[,2])=='logical') dataCol[,2] <- 'T'
  dataCol[,cols] <- apply(dataCol[,cols], 2, function(x) as.character(x))
  dataCol <- data.frame(sample=rownames(dataCol),
                        patient=dataCol$patient,
                        tumor = dataCol$tumor,
                        harbor_fusion = dataCol$harbor_fusion,
                        stringsAsFactors = F)
  
  dataAll8 <- read.delim(all8File,stringsAsFactors=F)
  pidG3p <- dataAll8[(dataAll8$Tumor==cancer) & (dataAll8$gene3p==gene3p) & (dataAll8$Sample=='T'), 'PatientID']
  pidG3p <- unique(pidG3p)
  
  dataCol[(dataCol$patient %in% pidG3p) & (dataCol$tumor=='T') & (dataCol$harbor_fusion=='N'),
          'harbor_fusion'] <- 'O'
  
  idxY <- which(dataCol$tumor=='T' & dataCol$harbor_fusion=='Y')
  patientY <- unique(dataCol[idxY,'patient'])
  if(length(patientY)<3){
    message('Less than 3 patients harboring fusion, cluster will not be computed.')
    return(NULL)
  }
  idxN <- which(dataCol$tumor=='T' & dataCol$harbor_fusion=='N')
  resMidExpr <- res[res$baseMean>100 & res$baseMean<10000,]
  genes <- row.names(resMidExpr)[resMidExpr$padj<1e-3]
  genes <- genes[!is.na(genes)]
  # remove the two genes that form fusion
  genes <- setdiff(genes,c(gene5p,gene3p))
  if(length(genes)==0){
    warning('No significantly differentially expressed genes.')
    return(NULL)
  }
  if(length(genes)<3){
    warning('Not enough number of significantly differentially expressed genes.')
    return(NULL)
  }
  select <- c(idxY,idxN)
  ddsSel <- dds[genes,select]
  logCounts <- log2(counts(ddsSel,normalized=TRUE)+1)
  r <- cor(logCounts)
  label <- c(dataCol[idxY,1], dataCol[idxN,1])
  dimnames(r) <- list(as.vector(label),as.vector(label))
  
  numY <- length(idxY)
  rYY <- as.vector(upperTriangle(r[1:numY,1:numY]))
  rYN <- as.vector(r[1:numY,-c(1:numY)])
  return(t.test(rYY,rYN))
}


## Generate table of clustering
clusterTable <- function(dataFile, cancer, gene5p, gene3p, all8File='/Users/HSH/Roche/Data/all8.txt'){
  require(DESeq2)
  
  load(dataFile)
  dataCol <- data.frame(dds@colData@listData)
  cols <- c(1,2,3,5)
  if(class(dataCol[,2])=='logical') dataCol[,2] <- 'T'
  dataCol[,cols] <- apply(dataCol[,cols], 2, function(x) as.character(x))
  dataCol <- data.frame(sample=rownames(dataCol),
                        patient=dataCol$patient,
                        tumor = dataCol$tumor,
                        harbor_fusion = dataCol$harbor_fusion,
                        stringsAsFactors = F)
  
  dataAll8 <- read.delim(all8File,stringsAsFactors=F)
  pidG3p <- dataAll8[(dataAll8$Tumor==cancer) & (dataAll8$gene3p==gene3p) & (dataAll8$Sample=='T'), 'PatientID']
  pidG3p <- unique(pidG3p)
  
  dataCol[(dataCol$patient %in% pidG3p) & (dataCol$tumor=='T') & (dataCol$harbor_fusion=='N'),
          'harbor_fusion'] <- 'O'
  
  idxY <- which(dataCol$tumor=='T' & dataCol$harbor_fusion=='Y')
  idxN <- which(dataCol$tumor=='T' & dataCol$harbor_fusion=='N')
  resMidExpr <- res[res$baseMean>100 & res$baseMean<10000,]
  genes <- row.names(resMidExpr)[resMidExpr$padj<1e-3]
  genes <- genes[!is.na(genes)]
  # remove the two genes that form fusion
  genes <- setdiff(genes,c(gene5p,gene3p))
  if(length(genes)==0){
    warning('No significantly differentially expressed genes.')
    return(NA)
  }
  if(length(genes)<3){
    warning('Not enough significantly differentially expressed genes.')
    return(NA)
  }
  select <- c(idxY,idxN)
  ddsSel <- dds[genes,select]
  logCounts <- log2(counts(ddsSel,normalized=TRUE)+1)
  r <- cor(logCounts)
  label <- c(dataCol[idxY,1], dataCol[idxN,1])
  dimnames(r) <- list(as.vector(label),as.vector(label))
  
  d <- 1-r
  hc <- hclust(as.dist(d),method = 'complete')
  numY <- length(idxY)
  nSample <- nrow(d)
  sampleHarbor <- rep('N',nSample)
  sampleHarbor[1:numY] <- 'Y'
  sampleHarbor <- sampleHarbor[hc$order]
  sampleID <- rownames(d)
  sampleID <- sampleID[hc$order]
  result <- data.frame(Index=1:nSample,SampleID=sampleID,HarborFusion=sampleHarbor,
                       stringsAsFactors=F)
  result
}


## Cluster by gene expression only, i.e. no information about fusion
cluster.by.expr <- function(exp_file,col_file){
  require(DESeq2)
  # create count matrix
  dat <- read.table(exp_file,header=T,sep='\t')
  countdata <- dat[,-1]
  rownames(countdata) <- dat[,1]
  colnames(countdata) <- gsub('.','-',colnames(countdata),fixed=T)
  countdata <- round(countdata)
  
  coldata <- data.frame(Dummy=base::sample(c('Yes','No'),ncol(countdata),replace=T))
  rownames(coldata) <- colnames(countdata)
  
  # construct the DESeqDataSet
  ddsFullCountTable <- DESeqDataSetFromMatrix(countData = countdata, colData=coldata, design=~Dummy)
  dds <- ddsFullCountTable
  message('Estimating size factors')
  dds <- estimateSizeFactors(dds)
  logCounts <- log2(counts(dds,normalized=TRUE)+1)
  gene_mean <- rowMeans(logCounts)
  gene_mid_expr <- rownames(logCounts)[(gene_mean>log2(1e2)) & (gene_mean<log2(1e4))]
  
  df <- data.frame(Gene=gene_mid_expr,Mean=NA,StandardDeviation=NA, stringsAsFactors=F)
  df$Mean <- gene_mean[gene_mid_expr]
  df$StandardDeviation <- apply(logCounts[gene_mid_expr,],1,sd)
  df$NormalizedSd <- df$StandardDeviation * 10 / df$Mean
  
  gene_sel <- df[df$NormalizedSd>2,'Gene']
  logCounts_sel <- logCounts[gene_sel,]
  r <- cor(logCounts_sel)
  label <- colnames(logCounts_sel)
  dimnames(r) <- list(as.vector(label),as.vector(label))
  d <- 1-r
  hc <- hclust(as.dist(d),method = 'complete')
  dat_clust <- data.frame(SampleID=colnames(r), Fusion=NA, stringsAsFactors=F)
  dat_clust <- dat_clust[hc$order,]
  
  # map fusion information back to cluster result
  dat <- read.table(col_file,header=T,sep='\t',stringsAsFactors=F)
  dat_col <- dat[,-1]
  rownames(dat_col) <- dat[,1]
  for(i in 1:nrow(dat_clust)) dat_clust[i,'Fusion'] <- dat_col[dat_clust[i,'SampleID'],'Fusion']
  dat_clust
}