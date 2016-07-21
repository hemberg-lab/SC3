expect_equal(ncol(SC3:::cell_filter(treutlein, 2000)), 50)
expect_equal(nrow(SC3:::cell_filter(treutlein, 2000)), 23271)

expect_equal(ncol(SC3:::gene_filter(treutlein, 0.06, 2, 0)), 80)
expect_equal(nrow(SC3:::gene_filter(treutlein, 0.06, 2, 0)), 7417)

expect_equal(ncol(SC3:::gene_filter(treutlein, 0.06, 0, 0)), 80)
expect_equal(nrow(SC3:::gene_filter(treutlein, 0.06, 0, 0)), 7929)

expect_equal(ncol(SC3:::gene_filter(treutlein, 0.06, 3, 0)), 80)
expect_equal(nrow(SC3:::gene_filter(treutlein, 0.06, 3, 0)), 7257)
