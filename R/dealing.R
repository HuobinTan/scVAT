#' Add cell properties
#'
#' @param vat vat entity
#' @param cell.prop cell property's values
#' @param prop.name cell property's name
#' @param meta.key cell property's meta feature key
#' @param meta.value property's meta feature value
#' @export
addCellProp <- function(vat, cell.prop, prop.name, meta.key=NULL, meta.value=NULL){
  if(nrow(vat@cell.props)!=length(cell.prop)){
    stop("The cell.prop length must equal vat@cell.props's length")
  }
  vat@cell.props[,prop.name] <- cell.prop
  vat@prop.meta <- AddPropMeta(vat@prop.meta, meta.key, meta.value)
  return(vat)
}

#' Add gene properties
#' @param vat vat entity
#' @param gene.prop gene property's values
#' @param prop.name gene property's name
#' @param meta.key property's meta feature key
#' @param meta.value property's meta feature value
#' @export
addGeneProp <- function(vat, gene.prop, prop.name, meta.key=NULL, meta.value=NULL){
  if(nrow(vat@gene.props)!=length(gene.prop)){
    stop("The gene.prop's length must equal vat@gene.props's length")
  }
  vat@gene.props[,prop.name] <- gene.prop
  vat@prop.meta <- AddPropMeta(vat@prop.meta, meta.key, meta.value)
  return(vat)
}

#' Add meta feature
#' @param prop.meta property name
#' @param meta.key meta feature key
#' @param meta.value meta feature value
#' @export
addPropMeta <- function(prop.meta, meta.key, meta.value){
  if(isEmpth(meta.key)) return(prop.meta)
  old.meta.value <- prop.meta[[meta.key]]
  if(is.null(old.meta.value)){
    prop.meta[[meta.key]] <- meta.value
  }else{
    is.in <- meta.value %in% old.meta.value
    no.in.idx <- meta.value[which(!is.in)]
    prop.meta[[meta.key]] <- c(old.meta.value, meta.value[no.in.idx])
  }
  return(prop.meta)
}

#' Get the property meta features
#' @param vat vat entity
#' @param meta.key meta feature key
#' @export
getPropMeta <- function(vat, meta.key){
  return(vat@prop.meta[[meta.key]])
}

#' Get dataset
#'
#' @param vat vat entity
#' @param use.genes use gene names
#' @param use.raw return raw data if TRUE, default FALSE
#' @param drop whether or not drop
#' @export
getUseData <- function(vat, use.genes = NULL, use.raw = FALSE, drop = FALSE){
  if (is.null(use.genes)){
    use.genes <- vat@use.genes
  }else{
    if(use.raw){
      avail.genes <- rownames(vat@data.raw)
    }else{
      avail.genes <- rownames(vat@data)
    }
    not.genes.pos <- !(use.genes %in% avail.genes)
    if(sum(not.genes.pos)>0){
      stop(sprintf("%s is/are not existence", paste0(use.genes[not.genes.pos],collapse=",")))
    }
  }
  if (use.raw){
    return(vat@data.raw[use.genes,,drop=drop])
  }else{
    return(vat@data[use.genes,,drop=drop])
  }
}

#' Get analysis result
#'
#' @param vat vat entity
#' @param key analysis key
#' @param cols analysis data column vector
#' @param as.data.frame convert to data.frame if TRUE, default FALSE
#' @export
getAnalysisData <- function(vat, key = "PCA", cols=c(1:50),as.data.frame=FALSE){
  if(isEmpty(key)){
    data <- getUseData(vat)
  } else{
    data <- eval(parse(text=paste0("vat@analysis$",key)))$cell.values[,cols]
  }
  if(as.data.frame) data <-  data.frame(data)
  return(data)
}

#' get available analysis key name
#' @param vat vat entity
#' @export
getAnalysisKey <- function(vat){
  return (names(bhsc@analysis))
}

#' get availbale analysis data colnames
#' @param vat vat entity
#' @param key analysis key, default tSNE
#' @export
getAnalysisColnames <- function(vat, key="tSNE"){
  return (colnames(eval(parse(text=paste0("vat@analysis$",key,"$cell.values")))))
}

#' Get the analysis key's colnames
#'
#' @param key analysis key
#' @param dims colnames subscript
#' @export
getAnalysisColName <- function(key="tSNE",dims=c(1:2)){
  if(key=="PCA"){
    return(paste("PC", dims, sep=""))
  } else{
    return(paste(key, dims, sep="_"))
  }
}

#' Get cell properties
#'
#' @param vat vat entity
#' @param prop.name cell property's name
#' @export
getCellPropData <- function(vat, prop.name){
  if(!(prop.name %in% colnames(vat@cell.props))){
    stop(paste0("There is not the ",prop.name, " in the cell.props"))
  }
  return(vat@cell.props[,prop.name])
}

#' Get gene properties
#'
#' @param vat vat entity
#' @param prop.name gene property's name
#' @export
getGenePropData <- function(vat, prop.name){
  if(!(prop.name %in% colnames(vat@gene.props))){
    stop(paste0("There is not the ",prop.name, " in the cell.props"))
  }
  return(vat@gene.props[,prop.name])
}

#' Calculating the expression's percent for some pattern's genes
#' @param vat vat entity
#' @param gene.pattern pattern expression for matching gene names
#' @param prop.name the saved column name of cell.props
#' @export
#' @examples
#' vat <- setGenePercent(vat, gene.pattern="^mt-", prop.name="mt.percent")
setGenePercent <- function(vat, gene.pattern, prop.name){
  gene.names <- vat@gene.props$name
  filter.genes <- grep(pattern=gene.pattern, gene.names, ignore.case = TRUE)
  gene.percent <- Matrix::colSums(vat@data[filter.genes,])/vat@cell.props$umi
  vat <- addCellProp(vat, gene.percent, prop.name, meta.key = "Filter", meta.value = prop.name)
  return(vat)
}

