set.seed(1)
A <- matrix(runif(100), ncol = 10)
res <- readRDS("res.rds")

B <- norm_laplacian(A)
expect_equal(res, B)