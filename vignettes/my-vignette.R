## ----knitr-options, echo=FALSE, message=FALSE, warning=FALSE-------------
library(knitr)
opts_chunk$set(fig.align = 'center', fig.width = 6, fig.height = 5, dev = 'png')

## ---- message=FALSE, warning=FALSE---------------------------------------
library(scater)
library(SC3)
treutlein[1:3, 1:3]

## ------------------------------------------------------------------------
# cell annotation
treutlein_cell_info <- data.frame(cell_id = colnames(treutlein))
cell_inds <- paste("Cell", 1:ncol(treutlein), sep = "_")
rownames(treutlein_cell_info) <- cell_inds
pd <- new("AnnotatedDataFrame", data = treutlein_cell_info)
# cell expression
treutlein_cell_exprs <- treutlein
colnames(treutlein_cell_exprs) <- cell_inds
# SCESEt object
treutlein_sceset <- newSCESet(fpkmData = treutlein_cell_exprs, phenoData = pd)

## ------------------------------------------------------------------------
is_exprs(treutlein_sceset) <- exprs(treutlein_sceset) > 0
treutlein_sceset <- calculateQCMetrics(treutlein_sceset)

## ------------------------------------------------------------------------
plotPCA(treutlein_sceset, colour_by = "cell_id")

## ------------------------------------------------------------------------
# Note that n_cores = 1 is required for compilation of this vignette.
# Please remove this parameter when running on your computer:
# treutlein_sceset <- sc3(treutlein_sceset, ks = 2:4, biology = TRUE)
treutlein_sceset <- sc3(treutlein_sceset, ks = 2:4, biology = TRUE, n_cores = 1)

## ---- eval=FALSE---------------------------------------------------------
#  sc3_interactive(treutlein_sceset)

## ----eval=FALSE----------------------------------------------------------
#  sc3_export_results_xls(treutlein_sceset)

## ------------------------------------------------------------------------
p_data <- pData(treutlein_sceset)
head(p_data[ , grep("sc3_", colnames(p_data))])

## ------------------------------------------------------------------------
plotPCA(
    treutlein_sceset, 
    colour_by = "sc3_3_clusters", 
    size_by = "sc3_3_log2_outlier_score"
)

## ------------------------------------------------------------------------
f_data <- fData(treutlein_sceset)
head(f_data[ , grep("sc3_", colnames(f_data))])

## ------------------------------------------------------------------------
plotFeatureData(
    treutlein_sceset, 
    aes(
        x = sc3_3_markers_clusts, 
        y = sc3_3_markers_auroc, 
        colour = sc3_3_markers_padj
    )
)

## ---- fig.height=6-------------------------------------------------------
sc3_plot_consensus(treutlein_sceset, k = 3)

## ---- fig.height=6, fig.width=8------------------------------------------
sc3_plot_consensus(
    treutlein_sceset, k = 3, 
    show_pdata = c(
        "cell_id", 
        "log10_total_features",
        "sc3_3_clusters", 
        "sc3_3_log2_outlier_score"
    )
)

## ------------------------------------------------------------------------
sc3_plot_silhouette(treutlein_sceset, k = 3)

## ---- fig.height=6-------------------------------------------------------
sc3_plot_expression(treutlein_sceset, k = 3)

## ---- fig.height=6, fig.width=8------------------------------------------
sc3_plot_expression(
    treutlein_sceset, k = 3, 
    show_pdata = c(
        "cell_id", 
        "log10_total_features",
        "sc3_3_clusters", 
        "sc3_3_log2_outlier_score"
    )
)

## ---- fig.height=3-------------------------------------------------------
sc3_plot_cluster_stability(treutlein_sceset, k = 3)

## ---- fig.height=9-------------------------------------------------------
sc3_plot_de_genes(treutlein_sceset, k = 3)

## ---- fig.height=9, fig.width=8------------------------------------------
sc3_plot_de_genes(
    treutlein_sceset, k = 3, 
    show_pdata = c(
        "cell_id", 
        "log10_total_features",
        "sc3_3_clusters", 
        "sc3_3_log2_outlier_score"
    )
)

## ---- fig.height=6-------------------------------------------------------
sc3_plot_markers(treutlein_sceset, k = 3)

## ---- fig.height=6, fig.width=8------------------------------------------
sc3_plot_markers(
    treutlein_sceset, k = 3, 
    show_pdata = c(
        "cell_id", 
        "log10_total_features",
        "sc3_3_clusters", 
        "sc3_3_log2_outlier_score"
    )
)

## ------------------------------------------------------------------------
# Note that n_cores = 1 is required for compilation of this vignette.
# Please remove this parameter when running on your computer:
# treutlein_sceset <- sc3_prepare(treutlein_sceset, ks = 2:4)
treutlein_sceset <- sc3_prepare(treutlein_sceset, ks = 2:4, n_cores = 1)
str(treutlein_sceset@sc3)

## ------------------------------------------------------------------------
treutlein_sceset <- sc3_estimate_k(treutlein_sceset)
str(treutlein_sceset@sc3)

## ------------------------------------------------------------------------
treutlein_sceset <- sc3_calc_dists(treutlein_sceset)
names(treutlein_sceset@sc3$distances)

## ------------------------------------------------------------------------
treutlein_sceset <- sc3_calc_transfs(treutlein_sceset)
names(treutlein_sceset@sc3$transformations)

## ------------------------------------------------------------------------
treutlein_sceset@sc3$distances

## ------------------------------------------------------------------------
treutlein_sceset <- sc3_kmeans(treutlein_sceset)
names(treutlein_sceset@sc3$kmeans)

## ------------------------------------------------------------------------
treutlein_sceset@sc3$transformations

## ------------------------------------------------------------------------
p_data <- pData(treutlein_sceset)
head(p_data[ , grep("sc3_", colnames(p_data))])

## ------------------------------------------------------------------------
treutlein_sceset <- sc3_calc_consens(treutlein_sceset)
names(treutlein_sceset@sc3$consensus)
names(treutlein_sceset@sc3$consensus$`3`)

## ------------------------------------------------------------------------
treutlein_sceset@sc3$kmeans

## ------------------------------------------------------------------------
p_data <- pData(treutlein_sceset)
head(p_data[ , grep("sc3_", colnames(p_data))])

## ------------------------------------------------------------------------
treutlein_sceset <- sc3_calc_biology(treutlein_sceset)

## ------------------------------------------------------------------------
p_data <- pData(treutlein_sceset)
head(p_data[ , grep("sc3_", colnames(p_data))])

## ------------------------------------------------------------------------
f_data <- fData(treutlein_sceset)
head(f_data[ , grep("sc3_", colnames(f_data))])

## ------------------------------------------------------------------------
no_svm_labels <- pData(treutlein_sceset)$sc3_3_clusters

## ------------------------------------------------------------------------
# Note that n_cores = 1 is required for compilation of this vignette.
# Please remove this parameter when running on your computer:
# treutlein_sceset <- sc3(treutlein_sceset, ks = 2:4, svm.num.cells = 50)
treutlein_sceset <- sc3(treutlein_sceset, ks = 2:4, biology = TRUE, svm_num_cells = 50, n_cores = 1)

## ------------------------------------------------------------------------
p_data <- pData(treutlein_sceset)
head(p_data[ , grep("sc3_", colnames(p_data))])

## ---- message=FALSE, warning=FALSE---------------------------------------
treutlein_sceset <- sc3_run_svm(treutlein_sceset)
p_data <- pData(treutlein_sceset)
head(p_data[ , grep("sc3_", colnames(p_data))])

## ------------------------------------------------------------------------
treutlein_sceset@sc3$svm_train_inds <- NULL
treutlein_sceset <- sc3_calc_biology(treutlein_sceset)
p_data <- pData(treutlein_sceset)
head(p_data[ , grep("sc3_", colnames(p_data))])

## ------------------------------------------------------------------------
svm_labels <- pData(treutlein_sceset)$sc3_3_clusters

## ------------------------------------------------------------------------
library(mclust)
adjustedRandIndex(no_svm_labels, svm_labels)

