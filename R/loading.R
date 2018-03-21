#' load data from csv file
#' @param filename csv file names
#' @param gene.names set rownames=gene.names, if null, using the first column as gene names. default null
#' @param cell.names set colnames=cell.names, if null, using the header as cell names. default null
#' @param gene.in.row gene store in row, default TRUE
#' @param to.sparseMatrix if TRUE, saving data as sparseMatrix, otherwise matrix. default TRUE
#' @param ... see also read.csv
#'
#' @return data matrix or sparseMatrix data
#' @import Matrix
#'
#' @examples
#' data <- loadCSVData("./filename.csv")
loadCSVData <- function(filename, gene.names=NULL, cell.names=NULL, gene.in.row=TRUE, to.sparseMatrix=TRUE,...){
  data <- read.csv(filename, ...)
  if(!gene.in.row) data <- t(data)
  if(is.null(gene.names)){
    gene.names <- data[,1]
    data <- data[,-1]
  }
  data <- as.matrix(data)
  if(to.sparseMatrix){
    data <- Matrix::Matrix(data, sparse=TRUE)
    #data <- as(data, "sparseMatrix")
  }

  rownames(data) <- gene.names
  if(!is.null(cell.names)){
    colnames(data) <- cell.names
  }
  return(data)
}

#' load data from H5 file
#' @param filename h5 file name
#' @param genome h5group genome name, default NULL, returns first group
#' @return data sparseMatrix format
#' @import Matrix
#' @importFrom rhdf5 h5read
#'
#' @examples
#' data <- loadH5Data("./filtered_gene_bc_matrices_h5.h5",genome="mm10",filter=TRUE, log=FLASE)
loadH5Data <- function(filename, genome=NULL)
{
  if(is.null(genome)){
    dset <- h5read(filename,"/")[[1]]
  }else{
    dset <- h5read(filename, genome)
  }

  data <- Matrix::sparseMatrix(i=dset$indices + 1,p=dset$indptr,
                       x=dset$data, dims=dset$shape,
                       dimnames = list(dset$gene_names,dset$barcodes))
  return(data)
}

#' initialize VAT Entity
#' @param data.raw data matrix
#' @param title dataset name
#' @param min.genes Include cells where more than this many genes are detected, default 0
#' @param min.cells Include genes where more than this many cells are detected, default 0
#' @param cell.scale whether or not normalizing by library size, default TRUE
#' @param scale.factor multiple factor for cell.scale, default 0, using the mean for library size
#' @param log.trans whether or not log-transform (log(p+pseudocount)), default TRUE
#' @param pseducount pseducount for log-transform, default 1
#' @param verbose whether or not printing verbose, default TRUE
#'
#' @import Matrix
initVATEntity <- function(data.raw, title = "VAT", min.genes = 0, min.cells = 0,
                          cell.scale = FALSE, scale.factor = 0,
                          log.trans = TRUE, pseudocount = 1,
                          verbose = TRUE){
  if(verbose) logging("create VAT Entity...")
  entity <- new(Class="VATEntity", data.raw=as(data.raw,"sparseMatrix"), title=title)
  use.genes <- as.vector(rownames(entity@data.raw))
  use.cells <- as.vector(colnames(entity@data.raw))

  gene.nums <- apply(entity@data.raw, 2, function(x){return (sum(x>0))})  # expressed gene's number for eache cell
  umi.nums <- Matrix::colSums(entity@data.raw) #apply(entity@data.raw, 2, function(x){return (sum(x))}) # expressed gene's UMI for each cell
  cell.nums <- apply(entity@data.raw, 1, function(x){return (sum(x>0))}) # cell's number for each gene

  #filtering cells and genes
  if(verbose) logging("Filtering cells and genes...")
  if(min.genes >= 0){
    use.cells <- use.cells[which(gene.nums > min.genes)]
  }
  if(min.cells >= 0){
    use.genes <- use.genes[which(cell.nums > min.cells)]
  }
  entity@data.raw <- entity@data.raw[use.genes, use.cells]

  #set gene and cell properties
  entity@use.genes <- use.genes
  entity@gene.props <- data.frame(name=use.genes, cell.nums=cell.nums[use.genes])
  entity@cell.props <- data.frame(name=use.cells, gene.nums=gene.nums[use.cells],umi=umi.nums[use.cells])
  entity@cell.props$manual.cluster <- 0
  entity@prop.meta[["Group"]] <- c("manual.cluster")

  #Normalize the data using log and scale UMI
  if(verbose) logging("Scaling and Normalizing data...")
  entity@data <- entity@data.raw
  if(cell.scale){
    index <- c(1:ncol(entity@cell.props))
    if(scale.factor==0){
      scale.factor <- mean(entity@cell.props$umi.nums)
    }

    if(verbose) ulapply <- pbapply::pblapply
    else ulapply <- lapply
    entity@data <- ulapply(index, function(x){return(entity@data[,x] * scale.factor /entity@cell.props$umi.nums[x])})

    entity@data <- do.call("cbind", entity@data)
    entity@data <- as(entity@data,"sparseMatrix")
  }
  if(log.trans){
    entity@data <- log(entity@data + pseudocount)
    entity@data <- as(entity@data,"sparseMatrix")
  }
  if(verbose) logging("It's done!")
  return(entity)
}

#' load analysis result manually
#' @param vat VATEntity object
#' @param filename csv filename
#' @param ndims the dim of analysis result
#' @param key key value in the anlaysis list
#'
#' @return new vat entity including new analysis result
loadAnalsisFromCSV <- function(vat, filename, ndims = 2, key = "tSNE", ...){
  data <- read.csv(filename,...)
  if((ndims + 1) == ncol(data)){
    rownames <- data[,1]
    data <- data[,-1]
    data <- as.matrix(data)
    rownames(data) <- rownames
  }
  if(nrow(data)!=nrow(vat@cell.props)){
    stop("The nrow of analsis data must equal with the nrow of vat@cell.props")
  }
  data <- data[,c(1:ndims)]
  colnames(data) <- getAnalysisColName(key,c(1:ndims))
  data <- list(cell.values=data)
  eval(parse(text=paste0("vat@analysis$",key," <- data")))
  return(vat)
}

