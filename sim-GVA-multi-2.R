source("share-func.R")
cores.num <- 50


# parameter setting -------------------------------------------------------
# Number of exposures
K <- 2
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



# GVA simulation --------------------------------------------------------------
N = 1000
# N = 2000


for (target_event in c(0.5, 0.3, 0.1)) {
  for (target_pi_val in c(0.9, 0.5, 0.3, 0.1)) {
    for (beta_val in c(0, log(1.5), log(2))) {
      sim.inc.f(
        sim.num = 1000,
        N = N,
        K = K,
        visit_times = c(0, 1, 2, 3),
        gamma1 = rep(0.5, K),
        target_pi = rep(target_pi_val, K),
        delta = matrix(rep(c(0.3, 0.2), K), nrow = K, byrow = TRUE),
        sigma_eps = rep(0.5, K),
        Sigma = Sigma,
        z_means = rep(0, K),
        z_sds = rep(1, K),
        target_event = target_event,
        weib_shape = 1.1,
        omega = 0.2,
        beta_X = rep(beta_val, K),
        cores.num = cores.num
      )
    }
  }
}


for (target_event in c(0.5, 0.3, 0.1)) {
  sim.inc.f(sim.num = 1000, N = N, K = K, visit_times = c(0, 1, 2, 3), gamma1 = rep(0.5, K),
            target_pi = c(0.9, 0.1),
            delta = matrix(rep(c(0.3, 0.2), K), nrow = K, byrow = TRUE), sigma_eps = rep(0.5, K), Sigma = Sigma, z_means = rep(0, K), z_sds = rep(1, K),
            target_event = target_event,
            weib_shape = 1.1, omega = 0.2,
            beta_X = c(0, 0), 
            cores.num = cores.num)
  
  sim.inc.f(sim.num = 1000, N = N, K = K, visit_times = c(0, 1, 2, 3), gamma1 = rep(0.5, K),
            target_pi = c(0.9, 0.1),
            delta = matrix(rep(c(0.3, 0.2), K), nrow = K, byrow = TRUE), sigma_eps = rep(0.5, K), Sigma = Sigma, z_means = rep(0, K), z_sds = rep(1, K),
            target_event = target_event,
            weib_shape = 1.1, omega = 0.2,
            beta_X = c(log(2), log(2)), 
            cores.num = cores.num)
  
  sim.inc.f(sim.num = 1000, N = N, K = K, visit_times = c(0, 1, 2, 3), gamma1 = rep(0.5, K),
            target_pi = c(0.9, 0.1),
            delta = matrix(rep(c(0.3, 0.2), K), nrow = K, byrow = TRUE), sigma_eps = rep(0.5, K), Sigma = Sigma, z_means = rep(0, K), z_sds = rep(1, K),
            target_event = target_event,
            weib_shape = 1.1, omega = 0.2,
            beta_X = c(0, log(1.5)), 
            cores.num = cores.num)
  
  sim.inc.f(sim.num = 1000, N = N, K = K, visit_times = c(0, 1, 2, 3), gamma1 = rep(0.5, K),
            target_pi = c(0.9, 0.1),
            delta = matrix(rep(c(0.3, 0.2), K), nrow = K, byrow = TRUE), sigma_eps = rep(0.5, K), Sigma = Sigma, z_means = rep(0, K), z_sds = rep(1, K),
            target_event = target_event,
            weib_shape = 1.1, omega = 0.2,
            beta_X = c(log(1.5), 0),
            cores.num = cores.num)
}


# Summary --------------------------------------------------------------


file.out <- "summary.semi.weib.bias.out"
if ( file.exists(as.character(file.out)) ) { unlink(as.character(file.out)) }

for (N in c(1000, 2000)) { 
  for (target_event in c(0.5, 0.3, 0.1)) {
    for (target_pi_val in c(0.9, 0.5, 0.3, 0.1)) { 
      for (beta_val in c(0, log(1.5), log(2))) {
        summary.semi.gva.bias.f(N=N, K=K, target_pi=rep(target_pi_val, K), target_event=target_event, sigma_eps = rep(0.5, K), 
                                omega=0.2, weib_shape=1.1, beta_X = rep(beta_val, K), 
                                out_dir="results", file.out=file.out)
      }
    }
    
    summary.semi.gva.bias.f(N=N, K=K, target_pi=c(0.9, 0.1), target_event=target_event, sigma_eps = rep(0.5, K), 
                            omega=0.2, weib_shape=1.1, beta_X = c(0, 0), 
                            out_dir="results", file.out=file.out)
    summary.semi.gva.bias.f(N=N, K=K, target_pi=c(0.9, 0.1), target_event=target_event, sigma_eps = rep(0.5, K), 
                            omega=0.2, weib_shape=1.1, beta_X = c(log(2), log(2)), 
                            out_dir="results", file.out=file.out)
    summary.semi.gva.bias.f(N=N, K=K, target_pi=c(0.9, 0.1), target_event=target_event, sigma_eps = rep(0.5, K), 
                            omega=0.2, weib_shape=1.1, beta_X = c(0, log(1.5)), 
                            out_dir="results", file.out=file.out)
    summary.semi.gva.bias.f(N=N, K=K, target_pi=c(0.9, 0.1), target_event=target_event, sigma_eps = rep(0.5, K), 
                            omega=0.2, weib_shape=1.1, beta_X = c(log(1.5), 0), 
                            out_dir="results", file.out=file.out)
    
  }
}




file.out <- "summary.semi.weib.se.out"
if ( file.exists(as.character(file.out)) ) { unlink(as.character(file.out)) }

for (N in c(1000, 2000)) { 
  for (target_event in c(0.5, 0.3, 0.1)) {
    for (target_pi_val in c(0.9, 0.5, 0.3, 0.1)) { 
      for (beta_val in c(0, log(1.5), log(2))) {
        summary.semi.gva.se.f(N=N, K=K, target_pi=rep(target_pi_val, K), target_event=target_event, sigma_eps = rep(0.5, K), 
                                omega=0.2, weib_shape=1.1, beta_X = rep(beta_val, K), 
                                out_dir="results", file.out=file.out)
      }
    }
    summary.semi.gva.se.f(N=N, K=K, target_pi=c(0.9, 0.1), target_event=target_event, sigma_eps = rep(0.5, K), 
                            omega=0.2, weib_shape=1.1, beta_X = c(0, 0), 
                            out_dir="results", file.out=file.out)
    summary.semi.gva.se.f(N=N, K=K, target_pi=c(0.9, 0.1), target_event=target_event, sigma_eps = rep(0.5, K), 
                            omega=0.2, weib_shape=1.1, beta_X = c(log(2), log(2)), 
                            out_dir="results", file.out=file.out)
    summary.semi.gva.se.f(N=N, K=K, target_pi=c(0.9, 0.1), target_event=target_event, sigma_eps = rep(0.5, K), 
                            omega=0.2, weib_shape=1.1, beta_X = c(0, log(1.5)), 
                            out_dir="results", file.out=file.out)
    summary.semi.gva.se.f(N=N, K=K, target_pi=c(0.9, 0.1), target_event=target_event, sigma_eps = rep(0.5, K), 
                            omega=0.2, weib_shape=1.1, beta_X = c(log(1.5), 0), 
                            out_dir="results", file.out=file.out)
  }
}



file.out <- "summary.semi.weib.beta.out"
if ( file.exists(as.character(file.out)) ) { unlink(as.character(file.out)) }

for (N in c(1000, 2000)) { 
  for (target_event in c(0.5, 0.3, 0.1)) {
    for (target_pi_val in c(0.9, 0.5, 0.3, 0.1)) { 
      for (beta_val in c(0, log(1.5), log(2))) {
        summary.semi.gva.betaonly.f(N=N, K=K, target_pi=rep(target_pi_val, K), target_event=target_event, sigma_eps = rep(0.5, K), 
                                    omega=0.2, weib_shape=1.1, beta_X = rep(beta_val, K), 
                                    out_dir="results", file.out=file.out)
      }
    }
    summary.semi.gva.betaonly.f(N=N, K=K, target_pi=c(0.9, 0.1), target_event=target_event, sigma_eps = rep(0.5, K), 
                            omega=0.2, weib_shape=1.1, beta_X = c(0, 0), 
                            out_dir="results", file.out=file.out)
    summary.semi.gva.betaonly.f(N=N, K=K, target_pi=c(0.9, 0.1), target_event=target_event, sigma_eps = rep(0.5, K), 
                            omega=0.2, weib_shape=1.1, beta_X = c(log(2), log(2)), 
                            out_dir="results", file.out=file.out)
    summary.semi.gva.betaonly.f(N=N, K=K, target_pi=c(0.9, 0.1), target_event=target_event, sigma_eps = rep(0.5, K), 
                            omega=0.2, weib_shape=1.1, beta_X = c(0, log(1.5)), 
                            out_dir="results", file.out=file.out)
    summary.semi.gva.betaonly.f(N=N, K=K, target_pi=c(0.9, 0.1), target_event=target_event, sigma_eps = rep(0.5, K), 
                            omega=0.2, weib_shape=1.1, beta_X = c(log(1.5), 0), 
                            out_dir="results", file.out=file.out)
  }
}

