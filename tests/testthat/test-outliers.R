d <- SC3:::cell_filter(treutlein, 2000)
d <- SC3:::gene_filter(d, 0.06)

d1 <- d[ , colnames(d) == "4"]
res <- SC3:::get_outl_cells(d1, colnames(d1))
expect_is(res, "numeric")
expect_equal(length(res), 13)
expect_equal(unique(names(res)), "4")

d2 <- d[ , c(1, 2, 4)]
res <- SC3:::get_outl_cells(d2, colnames(d2))
expect_is(res, "numeric")
expect_equal(length(res), 3)
expect_equal(unique(names(res)), c("4", "2"))
