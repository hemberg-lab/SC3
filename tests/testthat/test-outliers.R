
d <- SC3:::cell_filter(treutlein, 2000)
d <- SC3:::gene_filter(d, 0.06)

res <- SC3:::outl_cells_main(d[ , colnames(d) == "4"], 0.9999)
expect_is(res, "numeric")
expect_equal(length(res), 13)
expect_equal(unique(names(res)), "4")

res <- SC3:::outl_cells_main(d[ , c(1, 2, 4)], 0.9999)
expect_is(res, "numeric")
expect_equal(length(res), 3)
expect_equal(unique(names(res)), c("4", "2"))
