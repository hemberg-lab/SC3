
expect_equal(ncol(SC3:::cell_filter(treutlein, 2000)), 50)
expect_equal(nrow(SC3:::cell_filter(treutlein, 2000)), 23271)

expect_equal(ncol(SC3:::gene_filter(treutlein, 0.06)), 80)
expect_equal(nrow(SC3:::gene_filter(treutlein, 0.06)), 7417)
