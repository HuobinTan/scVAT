#' plot tsne maps
#' 
#' @param vat VAT Entity
#' @param title plot's title
#' @param group.id cell property for group
#' @param gradient whether or not with gradient
#' @param colors color palette
#' @param source plot's data source for shinyui
#' @param ... other parameters for plot_ly
#' 
#' @return plot.ly object
#' 
#' @example plotTSNE(vat, title="Example", group.id="cluster", gradient=FALSE)
plotTSNE <- function(vat, title=NULL, group.id = NULL, gradient = FALSE, colors = NULL, source="tsne",...){
  group.data = NULL
  if(!is.null(group.id)){
    group.data <- getCellPropData(vat, group.id)
  }
  plotAnalysis(vat=vat,dims=c(1:2), title=title,key="tSNE", color.data=group.data, gradient=gradient,colors=colors,source=source,...)
}

#' plot PC maps
plotPC <- function(vat, x.pc, y.pc, z.pc=NULL, title=NULL, color="black", source="pca",...){
  plotAnalysis(vat=vat,dims=c(x.pc, y.pc, z.pc), title=title, key="PCA",colors=color,...)
}

#' plot gene expression based on the analysis result
#' @param vat
#' @param genes
#' @param dims length(dim)=2, 3 (2D or 3D)
#' @param key analysis data key 
#' @param gradient default FALSE
#' @param colors color palette
#' @param show.title whether or not show title as gene's name
#' @param show.entity.title whether or not show entity name in the title
#' 
plotGene <- function(vat, genes, dims=c(1,2), key = "tSNE", 
                     gradient = TRUE, colors=c("lightgrey","blue"),
                     show.title=TRUE, show.entity.title=TRUE, ...){
  gene.data <- getUseData(vat,use.genes = genes,drop=FALSE)
  gene.data <- colSums(gene.data)
  if(!gradient){
    gene.data[which(gene.data > 0)] <- 1
    gene.data <- factor(gene.data,levels=c(0,1))
  } 
  title <- NULL
  if(show.title){
    title <- paste0(genes,sep=" ",collapse = "")
  }
  if(show.entity.title){
    if(is.null(title)){
      title <- vat@title
    }else{
      title <- paste0(vat@title, "'s ", title)
    }
  }
  plotAnalysis (vat, dims = dims, title=title, key = key, color.data = gene.data, 
                           gradient = gradient, colors = colors, ...)
}

#' plot some gene expressions based on the analysis result
#' @param vat
#' @param genes
#' @param nrows number of rows for plot
#' @param dims length(dim)=2, 3 (2D or 3D)
#' @param key analysis data key 
#' @param gradient default FALSE
#' @param colors color palette
#' 
plotMultipleGenes <- function(vat, genes, nrows=2, dims=c(1,2), key = "tSNE", gradient = TRUE, colors=c("lightgrey","blue"),...){
  plots <- lapply(genes, FUN=function(gene){
    p <- plotGene(vat,gene,dims=dims,key=key,gradient=gradient,colors=colors,...)
  })
  title <- paste0(vat@title,"'s ", paste(genes,sep=", ",collapse=", "))
  p <- subplot(plots,nrows=nrows) %>% layout(showlegend=FALSE) %>% layout(title = title)
  p$elementId <- NULL
  p
}
#' plot two gene expression based on the analysis result at the sametime
#' 
#' @param vat
#' @param gene1
#' @param gene2
#' @param dims length(dim)=2, 3 (2D or 3D)
#' @param key analysis data key 
#' @param gradient default FALSE
#' @param colors color palette
#' 
plotTwoGenes <- function(vat, gene1, gene2, dims=c(1,2), key = "tSNE",
                         colors=c("lightgrey","blue","orange","red"),...){
  
  if(length(colors)!=4){
    stop("The length of colors must equal 4.")
  }
  gene1.pos <- which(getUseData(vat,gene1,drop=TRUE) > 0)
  gene2.pos <- which(getUseData(vat,gene2,drop=TRUE) > 0)
  both.pos <- intersect(gene1.pos, gene2.pos)
  
  gene.data <- rep("None",times=nrow(vat@cell.props))
  gene.data[gene1.pos] <- gene1
  gene.data[gene2.pos] <- gene2
  gene12 <- paste(gene1, gene2, sep=", ") 
  gene.data[both.pos] <- gene12
  gene.data <- factor(gene.data,levels=c("None",gene1,gene2,gene12))
  title <- gene12
  
  plotAnalysis (vat, dims = dims, title=title, key = key, color.data = gene.data, 
                gradient = FALSE, colors = colors,...)
}

#plot using analysis data
#' @param vat
#' @param dims data analysis, dimension index, default c(1:2)
#' @param title plot title, default NULL, name as vat@title, if "", no title
#' @param key analysis data key value
#' @param color.data group data for color
#' @param gradient default FALSE
#' @param colors color palette, if no color.data, uses colors[1]
#' @param source plot id
#' 
#' @return p
plotAnalysis <- function(vat, dims = c(1:2), title=NULL, key="tSNE", color.data = NULL, 
                         gradient = FALSE, colors = NULL, source="src",...){
  if(!(key %in% names(vat@analysis))){
    stop(sprintf("No %s result, please run analysis program first",key))
  }
  if(length(dims) < 2 || length(dims) > 3){
    stop("Only plot 2D or 3D maps")
  }
  twoD <- (length(dims)==2)
  if(!is.null(title) && title=="") title <- vat@title
  plot.data <- getAnalysisData(vat, key, cols =dims,as.data.frame = TRUE)
  if(twoD){
    colnames(plot.data) <- c("X","Y")
  }else{
    colnames(plot.data) <- c("X","Y","Z")
  }
  if(!is.null(color.data)){
    if(!gradient & !is.factor(color.data)) color.data <- factor(color.data)
    plot.data$color <- color.data
  }
  if(is.null(color.data)){
    if(!is.null(colors)) colors = I(colors[1])
    if(twoD){
      p <- plot_ly(plot.data, x = ~X, y = ~Y, source=source, color=colors,...)
    }else{
      p <- plot_ly(plot.data, x = ~X, y = ~Y, z = ~Z, source=source, color=colors,...)
    }
  }else{
    if(twoD){
      p <- plot_ly(plot.data, x = ~X, y = ~Y, color = ~color, source=source, colors=colors,...) 
    }else{
      p <- plot_ly(plot.data, x = ~X, y = ~Y, z = ~Z, color = ~color, source=source, colors=colors,...) 
    }
  }
  axis.names <- getAnalysisColName(key, dims)
  p <- p %>% add_markers() %>% layout(dragmode="select")
  if(twoD){
    p <- layout(p,
                xaxis = list(title = axis.names[1]),
                yaxis = list(title = axis.names[2]))
  } else {
    p <- layout(p,
                scene = list(xaxis = list(title = axis.names[1]),
                             yaxis = list(title = axis.names[2]),
                             zaxis = list(title = axis.names[3]))
    )
  }
  if(isEmpty(title)){
    p <- layout(p, title=title)
  }
  p$elementId <- NULL
  p
}