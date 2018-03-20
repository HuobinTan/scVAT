#' The VAT Entity Class
#' 
#' The VAT entity is the core data structure for VAT pipeline.
#' It stores all contents including data, annotations, analysis, etc.
#' 
#' @slot 
#' 
#' @name VATEntity
#' @rdname VATEntity
#' @exportClass VATEntity
#' @useDynLib scVAT

VATEntity <- methods::setClass(
  "VATEntity",
  slots = c(
    title = "character",
    data = "ANY",
    data.raw = "ANY",
    gene.props = "data.frame",
    cell.props = "data.frame",
    prop.meta = "list",
    use.genes = "vector",
    analysis = "list"
  )
)

#' show method for VATEntity
#'
#' @param object A VATEntity object
#' @name show
#' @docType methods
#' @rdname show-methods
#'
setMethod(
  f = "show",
  signature = "VATEntity",
  definition = function(object) {
    cat(
      "An object of ",class(object),
      "in dataset'",object@title, "'\n",
      nrow(object@data.raw),"genes, ",ncol(object@data.raw),"cells.\n"
    )
    invisible(NULL)
  }
)