library(DESeq2)
library('reactome.db')

# Check if the expression gct file and colData file have the same set of samples
# If not, rewrite the expression file to be consistent with the colData file
sanity.check <- function(exp_file,col_file){
  dataCol <- read.table(col_file,header=T,sep='\t',stringsAsFactors=F)
  dataCol[,1] <- gsub('.','-',dataCol[,1],fixed=T)
  sampleCol <- dataCol[,1]
  dataExpr <- read.table(exp_file,header=T,sep='\t',stringsAsFactors=F)
  sampleExpr <- colnames(dataExpr)[-1]
  sampleExpr <- gsub('.','-',sampleExpr,fixed=T)
  colnames(dataExpr) <- c('NAME',sampleExpr)
  colSel <- c(TRUE,sampleExpr %in% sampleCol)
  dataExprNew <- dataExpr[,colSel]
  write.table(dataExprNew,file=exp_file,quote=F,sep='\t',row.names=F)
  rowSel <- (sampleCol %in% sampleExpr)
  dataColNew <- dataCol[rowSel,]
  if(class(dataColNew[,3])=='logical') dataColNew[,3] <- 'T'
  write.table(dataColNew,file=col_file,quote=F,sep='\t',row.names=F)
  return
}

# Differential gene expression analysis
differentialExpressionAnalysis <- function(exp_file,col_file,save_data_file=NULL,
                                           result_file=NULL){
  # create count matrix
  dat <- read.table(exp_file,header=T,sep='\t')
  countdata <- dat[,-1]
  rownames(countdata) <- dat[,1]
  colnames(countdata) <- gsub('.','-',colnames(countdata),fixed=T)
  countdata <- round(countdata)
  # create column information
  dat <- read.table(col_file,header=T,sep='\t')
  coldata <- dat[,-1]
  rownames(coldata) <- dat[,1]
  
  # construct the DESeqDataSet
  ddsFullCountTable <- DESeqDataSetFromMatrix(
    countData = countdata,
    colData = coldata,
    design = ~harbor_fusion)
  dds <- ddsFullCountTable
  
  dds$harbor_fusion <- relevel(dds$harbor_fusion,'N')
  
  # Running the pipeline
  dds <- DESeq(dds)
  res <- results(dds)
  
  # Save the analysis results
  if(is.character(save_data_file)){
    save(dds,res,file=save_data_file)
  }
  if(is.character(result_file)){
    write.csv(res[order(res$padj),],file=result_file)
  }
  
  res
}

# Map between gene symbol and Entrez gene ID
symbol2id <- function(symbol,
                      file_path='/Users/HSH/Roche/Data/expression/gene_annotation_v1.txt'){
  dat <- read.table(file_path,header=T,sep='\t')
  symbol_all <- dat[,2]
  entrez_id <- as.character(dat[,1])
  names(entrez_id) <- symbol_all
  mapping <- as.list(entrez_id)
  IDs <- base::sapply(symbol, function (x) mapping[[x]])
  return(IDs)
}

# Perform t-test for pathway
testCategory <- function(reactomeID){
  isMember <- incm[reactomeID, ]
  data.frame(
    reactomeID = reactomeID,
    numGenes = sum(isMember),
    avgLFC = mean(res2$log2FoldChange[isMember]),
    strength = sum(res2$log2FoldChange[isMember])/sqrt(sum(isMember)),
    pvalue = t.test(res2$log2FoldChange[isMember])$p.value,
    reactomeName = reactomePATHID2NAME[[reactomeID]])}

# Pathway analysis
pathwayAnalysis <- function(object,result_file=NULL){
  object$entrez <- symbol2id(row.names(object))
  res2 <- object[object$entrez %in% keys(reactome.db,'ENTREZID') & !is.na(object$pvalue), ]
  reactomeTable <- AnnotationDbi::select(reactome.db,keys=res2$entrez,
                                         keytype='ENTREZID',
                                         columns=c('ENTREZID','REACTOMEID'))
  incm <- do.call(rbind,with(reactomeTable,
                             tapply(ENTREZID,factor(REACTOMEID),function(x) res2$entrez %in% x)))
  incm <- incm[rowSums(incm)>=5, ]
  reactomeResult <- do.call(rbind,lapply(rownames(incm),testCategory))
  reactomeResult$padjust <- p.adjust(reactomeResult$pvalue,'BH')
  reactomeResultSignif <- reactomeResult[reactomeResult$padjust < 0.05, ]
  result <- reactomeResultSignif[order(abs(reactomeResultSignif$strength),decreasing = TRUE), ]
  if(is.character(result_file)){
    write.csv(result,file=result_file)
  }
}

