#' do PCA using irlba library
#' 
#' @param vat 
#' @param pc.num default 50
#' @param use.genes run PCA using the genes. if NULL, using vat@use.genes. default NULL
doPCA <- function(vat, pc.num = 50, use.genes = NULL, ...){
  data <- getUseData(vat, use.genes=use.genes)
  pc.num <- min(pc.num, nrow(data) - 1)
  pca <- prcomp_irlba(t(data), n = pc.num, retx=TRUE, ...)
  rownames(pca$rotation) <- use.genes
  pca$gene.loadings <- pca$rotation
  pca$cell.values <- pca$x
  pca$rotation <- NULL
  pca$x <- NULL
  eval(parse(text=paste0("vat@analysis$","PCA","<-pca")))
  return(vat)
}

#' do tSNE using Rtsne library
#' 
#' @param vat vat entity
#' @param dims tSNE dimension
doTSNE <- function(vat, dims = 2, analysis.key = "PCA", use.col = 50, ...){
  data <- getAnalysisData(vat,key = analysis.key, cols =c(1:use.col))
  tsne <- Rtsne(as.matrix(data), dims = dims, ...)
  tsne <- tsne$Y
  paste("tSNE_", c(1:dims),sep="")
  colnames(tsne) <- getAnalysisColName("tSNE",c(1:dims))
  vat@analysis$tSNE$cell.values <- tsne
  return(vat)
}

