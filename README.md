# scVAT (single-cell RNA Visualization Analysis Toolkit)
scVAT is a Visualization Analysis Toolkit for single-cell RNA sequence data. scVAT implements a base pipeline for scRNA, including loading, preprocessing, filtering, PCA and tSNE etc. Besides, it provides:
- A Web GUI (Graphical User Interface) for visualizing all kinds of analysis results
- Graphical clustering manually based on analysis data
- Differentiating cells

## Install
You can install scVAT from GitHub with:
```{r}
#install.package("devtools")
devtools::install_github("HuobinTan/scVAT")
```

## Loading gene-cell expression matrix
There are three different functions for loading gene-cell expression matrix:
1) Loading .csv file
2) Loading 10X H5 file
3) Loading 10X data path
```{r}
library(scVAT)
# loading matrix form csv file
# gene.index: gene name's index
# cell.index: cell name's index
# gene.in.row: gene stores in a row
csv.data <- loadCSVData (filename = "./matrix.csv", gene.index = 1, cell.index = 1, gene.in.row=TRUE)

# loading matrix from 10x H5 file
# genome: genome in the H5, such as "mm10"", "GRch38", default NULL, the first genome 
h5.data <-load10XH5(filename = "./filtered_gene_bc_matrices_h5.h5", genome = NULL)

#loading matrix from 10x outs path, including barcodes.tsv, genes.tsv and matrix.mtx
path10x.data <- load10XPath(path = "./filtered_gene_bc_matrices/mm10")
```

## Initialize VAT Entity
Building VAT object based on the last expression matrix, Besides:
- Seting paremeters to filter genes and cells (min.genes, min.cells)
- Scaling the data using library size: if cell.scale is  TRUE, p <- (p / nUMI \* scale.factor) for each cell;  if scale.factor is 0,  scale.factor =  median(nUMI)
- Log-Transform: if log.trans, then p <- log(p + pesudocount)
```{r}
#initializing VAT entity
vat <- initVATEntity(data.raw = csv.data, title = "VAT", min.genes = 0, min.cells = 0,
                          cell.scale = TRUE, scale.factor = 0,
                          log.trans = TRUE, pseudocount = 1,
                          verbose = TRUE)
```

## Do analysis
```{r}
#running PCA
#PCA will be saved in the vat@analysis$PCA
vat <- doPCA(vat, pc.num = 50, use.genes = NULL)
#plotting PC's standard deviation, and selecting PC components for downstream
plotPCASDev(vat, key="PCA", ndims = 50)
#runing tSNE based on analysis data
vat <- doTSNE(vat, dims = 2, analysis.key = "PCA", use.col = 50)

#loading another analysis data (not run)
vat <- loadAnalysisFromCSV(vat, filename = "youranalysisfile", ndims = 3, key="YourKey")
```


##Visualize analysis data
```{r}
#plotting PCA result, set x.pc, y.pc. if 3D, set z.pc
plotPC(vat, x.pc = 1, y.pc = 2, z.pc = 3)
#plotting tSNE, and group.id for grouping
plotTSNE(vat, group.id = "manual.cluster")
#plotting one gene,two genes and three genes based on analysis data (parameter "key" sets visualization analysis data, the value is same as the key of loadAnalysisFromCSV function. PC result uses "PCA", tSNE result uses "tSNE")
plotGene(vat, genes="Cdc20",dims=c(1,2),key="tSNE")
plotTwoGenes(vat, gene1="Cdc20", gene2 = "Ube2c", dims = c(1,2), key="tSNE")
plotThreeGenes(vat, gene1 = "Cdc20", gene2 = "Ube2c", gene2="Gata1", dim=c(1,2), key="tSNE")

#starting Web GUI for visualizing, manually clustering, and differintialing
#be cautious, the parameter is a string of variable name, not variable
startVATGUI("vat")
```
