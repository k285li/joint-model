source("share-func-inc.R")
cores.num <- 10


# parameter setting -------------------------------------------------------
# Number of exposures
K <- 10
# Parameters for within-exposure correlation block
a <- 0.1
b <- 0.05
# 4x4 block for each exposure
Sigma_block <- matrix(c(
  a, b, b, b,
  b, a, b, b,
  b, b, a, b,
  b, b, b, a
), nrow = 4, byrow = TRUE)

# 4x4 off-diagonal block between exposures
Sigma_off <- diag(c(0.01, 0.01, 0.01, 0.01))

# Initialize list of blocks
block_list <- vector("list", K)

for (i in 1:K) {
  row_list <- vector("list", K)
  
  for (j in 1:K) {
    if (i == j) {
      row_list[[j]] <- Sigma_block
    } else {
      row_list[[j]] <- Sigma_off
    }
  }
  
  block_list[[i]] <- do.call(cbind, row_list)
}

# Final 4Kx4K matrix
Sigma <- do.call(rbind, block_list)
# eigen(Sigma)$values # if is not positive definite, reduce the absolute value of the off-diagonal elements that you suspect are too large, or increase the value of the diagonal elements.




# GVA simulation --------------------------------------------------------------
N = 1000

for (target_event in c(0.1, 0.3)) {
  sim.inc.f(sim.num = 1000, N = N, K = K, visit_times = c(0, 1, 2, 3), gamma1 = rep(0.5, K),
            target_pi = c(0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.1, 0.1), ######
            delta = matrix(rep(c(0.3, 0.2), K), nrow = K, byrow = TRUE),
            sigma_eps = rep(0.5, K), Sigma = Sigma, z_means = rep(0, K), z_sds = rep(1, K),
            target_event = target_event, ######
            weib_shape = 1.1, omega = 0.2,
            beta_X = c(0, 0, 0, 0, 0, 0, 0, 0, log(1.5), log(1.5)), ######
            cores.num = cores.num)
}











