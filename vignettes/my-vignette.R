## ----knitr-options, echo=FALSE, message=FALSE, warning=FALSE-------------
library(knitr)
opts_chunk$set(fig.align = 'center', fig.width = 6, fig.height = 5, dev = 'png')

## ---- message=FALSE, warning=FALSE---------------------------------------
library(scater)
library(SC3)
treutlein[1:3, 1:3]

## ------------------------------------------------------------------------
treutlein_cell_info <- data.frame(author_clusters = colnames(treutlein))
cell_inds <- paste("Cell", 1:ncol(treutlein), sep = "_")
rownames(treutlein_cell_info) <- cell_inds
treutlein_cell_exprs <- treutlein
colnames(treutlein_cell_exprs) <- cell_inds
pd <- new("AnnotatedDataFrame", data = treutlein_cell_info)
treutlein_sceset <- newSCESet(countData = treutlein_cell_exprs, phenoData = pd)

## ------------------------------------------------------------------------
treutlein_sceset <- calculateQCMetrics(treutlein_sceset)

## ------------------------------------------------------------------------
plotPCA(treutlein_sceset, colour_by = "author_clusters")

