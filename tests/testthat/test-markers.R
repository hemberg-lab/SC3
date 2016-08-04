d <- SC3:::cell_filter(treutlein, 2000)
d <- SC3:::gene_filter(d, 0.06, 2, 0)

res <- get_marker_genes(d, colnames(d))
expect_is(res, "data.frame")
expect_equal(dim(res)[1], 95)
expect_equal(dim(res)[2], 3)
