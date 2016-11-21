## ----knitr-options, echo=FALSE, message=FALSE, warning=FALSE-------------
library(knitr)
opts_chunk$set(fig.align = 'center', fig.width = 6, fig.height = 5, dev = 'png')

## ---- message=FALSE, warning=FALSE---------------------------------------
library(scater)
library(SC3)
treutlein[1:3, 1:3]

## ------------------------------------------------------------------------
treutlein_cell_info <- data.frame(cell_id = colnames(treutlein))
cell_inds <- paste("Cell", 1:ncol(treutlein), sep = "_")
rownames(treutlein_cell_info) <- cell_inds
treutlein_cell_exprs <- treutlein
colnames(treutlein_cell_exprs) <- cell_inds
pd <- new("AnnotatedDataFrame", data = treutlein_cell_info)
treutlein_sceset <- newSCESet(fpkmData = treutlein_cell_exprs, phenoData = pd)

## ------------------------------------------------------------------------
is_exprs(treutlein_sceset) <- exprs(treutlein_sceset) > 0
treutlein_sceset <- calculateQCMetrics(treutlein_sceset)

