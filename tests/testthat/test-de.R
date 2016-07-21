d <- SC3:::cell_filter(treutlein, 2000)
d <- SC3:::gene_filter(d, 0.06, 2, 0)

res <- get_de_genes(d, colnames(d))
expect_is(res, "numeric")
expect_equal(length(res), 87)
