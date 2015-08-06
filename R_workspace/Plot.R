require(DESeq2)
require(gdata)

get.plot.df <- function(gene5p,gene3p,cancer,dataFile, ghFile=NULL, mapGHfile=NULL,
                        all8File='~/Roche/Data/all8.txt'){
  source('~/Roche/R workspace/UserFunctions.R')
  
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
  if(length(genes)==0) return(NA)
  select <- c(idxY,idxN)
  ddsSel <- dds[genes,select]
  logCounts <- log2(counts(ddsSel,normalized=TRUE)+1)
  r <- cor(logCounts)
  label <- c(dataCol[idxY,1], dataCol[idxN,1])
  dimnames(r) <- list(as.vector(label),as.vector(label))
  numY <- length(idxY)
  rYY <- as.vector(upperTriangle(r[1:numY,1:numY]))
  rYN <- as.vector(r[1:numY,-c(1:numY)])
  
  # Modify the harbor_fusion column based on GeneHub data
  if(!missing(ghFile)){
    dataGHmap <- read.delim(mapGHfile,stringsAsFactors=F)
    grepStrGH <- dataGHmap[dataGHmap[,1]==cancer,2]
    dataGH <- get.pid(ghFile,gene3p,grepStrGH,gene5p=gene5p)
    pidGH_O <- dataGH[dataGH[,4]=='O',3]
    select <- (dataCol$patient %in% pidGH_O) & (dataCol$harbor_fusion == 'N')
    dataCol[select,'harbor_fusion'] <- 'O'
    pidGH_Y <- dataGH[dataGH[,4]=='Y',3]
    pidY <- substr(rownames(r)[1:numY],1,12)
    addPatientY <- setdiff(pidGH_Y,pidY)
    if(length(addPatientY)>0){
      for(p in addPatientY){
        select <- (dataCol$patient==p) & (dataCol$tumor=='T') & (dataCol$harbor_fusion != 'Y')
        if(any(select)){
          sID <- dataCol[select,'sample'][1]
          if(mean(r[sID,1:numY]) > mean(rYN)) dataCol[select,'harbor_fusion'] <- 'Y'
        }
      }
    }
  }
  
  # reselect samples
  idxY <- which(dataCol$tumor=='T' & dataCol$harbor_fusion=='Y')
  idxN <- which(dataCol$tumor=='T' & dataCol$harbor_fusion=='N')
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
  
  plotDF <- data.frame(fusion=rep(paste(gene5p,gene3p,sep='-'),nSample), num=rep(1,nSample), 
                       harbor=sampleHarbor)
  plotDF
}
