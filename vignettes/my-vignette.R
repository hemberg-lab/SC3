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
colnames(treutlein) <- cell_inds
pd <- new("AnnotatedDataFrame", data = treutlein_cell_info)
treutlein_sceset <- newSCESet(countData = treutlein, phenoData = pd)

## ------------------------------------------------------------------------
treutlein_sceset <- calculateQCMetrics(treutlein_sceset)

## ------------------------------------------------------------------------
plotPCA(treutlein_sceset, colour_by = "author_clusters")

## ---- message=FALSE, warning=FALSE---------------------------------------
# Note that n.cores = 1 is required for compilation of this vignette.
# Please remove this parameter when running on your computer:
# treutlein_sceset <- sc3(treutlein_sceset, ks = 2:4)
treutlein_sceset <- sc3(treutlein_sceset, ks = 2:4, n.cores = 1)

## ---- eval=FALSE---------------------------------------------------------
#  sc3_interactive(treutlein_sceset)

## ------------------------------------------------------------------------
treutlein_sceset <- sc3_summarise_results(treutlein_sceset, k = 3)

## ----eval=FALSE----------------------------------------------------------
#  sc3_export_results_xls(treutlein_sceset)

## ------------------------------------------------------------------------
plotPCA(treutlein_sceset, colour_by = "sc3_clusters")

## ---- fig.height=6-------------------------------------------------------
sc3_plot_consensus(treutlein_sceset, k = 3)

## ------------------------------------------------------------------------
sc3_plot_silhouette(treutlein_sceset, k = 3)

## ---- fig.height=6-------------------------------------------------------
sc3_plot_expression(treutlein_sceset, k = 3)

## ---- fig.height=3-------------------------------------------------------
sc3_plot_cluster_stability(treutlein_sceset, k = 3)

## ------------------------------------------------------------------------
sc3_plot_tsne(treutlein_sceset, k = 3)

## ---- fig.height=8-------------------------------------------------------
sc3_plot_de_genes(treutlein_sceset, k = 3)

## ---- fig.height=6-------------------------------------------------------
sc3_plot_markers(treutlein_sceset, k = 3)

## ---- fig.height=3-------------------------------------------------------
sc3_plot_cell_outliers(treutlein_sceset, k = 3)

## ------------------------------------------------------------------------
treutlein_sceset@sc3 <- list()

## ------------------------------------------------------------------------
# Note that n.cores = 1 is required for compilation of this vignette.
# Please remove this parameter when running on your computer:
# treutlein_sceset <- sc3_prepare(treutlein_sceset)
treutlein_sceset <- sc3_prepare(treutlein_sceset, n.cores = 1)
treutlein_sceset@sc3$processed_dataset[1:3, 1:3]
treutlein_sceset@sc3$kmeans_iter_max
treutlein_sceset@sc3$rand_seed
treutlein_sceset@sc3$kmeans_nstart
treutlein_sceset@sc3$n_dim
treutlein_sceset@sc3$n_cores
treutlein_sceset@sc3$rselenium

## ------------------------------------------------------------------------
treutlein_sceset <- sc3_estimate_k(treutlein_sceset)
treutlein_sceset@sc3$k_prediction

## ------------------------------------------------------------------------
treutlein_sceset <- sc3_set_ks(treutlein_sceset, ks = 2:4)
treutlein_sceset@sc3$ks

## ------------------------------------------------------------------------
treutlein_sceset <- sc3_calc_dists(treutlein_sceset)
names(treutlein_sceset@sc3$distances)

## ------------------------------------------------------------------------
treutlein_sceset <- sc3_calc_transfs(treutlein_sceset)
names(treutlein_sceset@sc3$transformations)

## ------------------------------------------------------------------------
treutlein_sceset <- sc3_kmeans(treutlein_sceset)
names(treutlein_sceset@sc3$kmeans)

## ------------------------------------------------------------------------
treutlein_sceset <- sc3_calc_consens(treutlein_sceset)
names(treutlein_sceset@sc3$consensus)

## ------------------------------------------------------------------------
treutlein_sceset <- sc3_calc_biology(treutlein_sceset)
head(treutlein_sceset@sc3$biology$`3`$markers)
head(treutlein_sceset@sc3$biology$`3`$de.genes)
head(treutlein_sceset@sc3$biology$`3`$cell.outl)

## ------------------------------------------------------------------------
treutlein_sceset <- sc3_summarise_results(treutlein_sceset, k = 3)
no.svm.labels <- treutlein_sceset@sc3$results$clusters$sc3_clusters

## ------------------------------------------------------------------------
treutlein_sceset@sc3 <- list()
# Note that n.cores = 1 is required for compilation of this vignette.
# Please remove this parameter when running on your computer:
# treutlein_sceset <- sc3(treutlein_sceset, ks = 2:4, svm.num.cells = 50)
treutlein_sceset <- sc3(treutlein_sceset, ks = 2:4, svm.num.cells = 50, n.cores = 1)

## ------------------------------------------------------------------------
treutlein_sceset <- sc3_summarise_results(treutlein_sceset, k = 3)
table(!is.na(treutlein_sceset@sc3$results$clusters$sc3_clusters))

## ------------------------------------------------------------------------
treutlein_sceset <- sc3_run_svm(treutlein_sceset, k = 3)
treutlein_sceset <- sc3_summarise_results(treutlein_sceset, k = 3)
svm.labels <- treutlein_sceset@sc3$results$clusters$sc3_clusters

## ------------------------------------------------------------------------
data.frame(no_SVM_сlusters = no.svm.labels, SVM_clusters = svm.labels)

## ---- eval=FALSE---------------------------------------------------------
#  RSelenium::checkForServer()

