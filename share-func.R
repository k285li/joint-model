library(MASS)     
library(survival) 
library(cubature)
library(dplyr)
library(Matrix)
library(numDeriv)
library(mvtnorm)
library(parallel)
library(lme4)
library(tidyr)

# ---- Set up ----
expit <- function(x) {
  1 / (1 + exp(-x))
}

baseline_cum_increment <- function(t0, t1, weib_lambda = NULL, weib_shape = NULL) {
  t0 <- pmax(t0, 0)
  t1 <- pmax(t1, t0)
  weib_lambda * (t1^weib_shape - t0^weib_shape)
}


baseline_loghaz <- function(t, weib_lambda = NULL, weib_shape = NULL) {
  log(pmax(weib_lambda, 1e-12)) + log(pmax(weib_shape, 1e-12)) + (weib_shape - 1) * log(t)
}


# Solve gamma0 so that E[ logistic(gamma0 + gamma1*Z + (b0 + b1*Z)) ] = target
# Z ~ N(muZ, sdZ^2), (b0, b1) ~ N(0, Sigma_b2x2). Uses 2-D cubature (over Z and a single Normal).
solve_gamma0_occurrence <- function(
    gamma1,             # scalar: fixed-effect slope for Z
    muZ, sdZ,           # mean/sd of Z
    Sigma_b2x2,         # 2x2 covariance of (b0, b1)
    target = 0.9,
    root_interval = c(-12, 12),
    tol = 1e-6
){
  stopifnot(is.numeric(gamma1), length(gamma1) == 1,
            is.numeric(muZ), is.numeric(sdZ), sdZ > 0,
            all(dim(Sigma_b2x2) == c(2, 2)))
  
  s11 <- Sigma_b2x2[1,1]; s12 <- Sigma_b2x2[1,2]; s22 <- Sigma_b2x2[2,2]
  inv_logit <- function(x) 1/(1 + exp(-x))
  
  # E[ pi ] as a function of gamma0
  expected_p <- function(gamma0){
    f <- function(u){
      eval_pt <- function(pt){
        z <- muZ + sdZ * qnorm(pt[1])      # Z ~ N(muZ, sdZ^2)
        w <- qnorm(pt[2])                   # w ~ N(0,1)
        s2 <- s11 + 2*z*s12 + (z^2)*s22     # Var(b0 + z*b1) = a' Σ a, a=(1, z)
        s  <- sqrt(pmax(s2, 0))
        eta <- gamma0 + gamma1*z + s*w
        inv_logit(eta)
      }
      if (is.matrix(u)) apply(u, 1L, eval_pt) else eval_pt(u)
    }
    res <- cubature::adaptIntegrate(f, c(0,0), c(1,1))
    as.numeric(res$integral)
  }
  
  g <- function(g0) expected_p(g0) - target
  
  a <- root_interval[1]; b <- root_interval[2]
  fa <- g(a); fb <- g(b); tries <- 0
  while (fa*fb > 0 && tries < 6) {
    a <- a - 5; b <- b + 5
    fa <- g(a); fb <- g(b); tries <- tries + 1
  }
  if (fa*fb > 0) stop("Could not bracket root; widen 'root_interval' or check inputs.")
  
  uniroot(g, c(a, b), tol = tol)$root
}



calibrate_weibull_lambda <- function(sim.dat, beta_X, omega,
                                     target, shape,
                                     gfun = function(x) ifelse(x > 0, log(x), 0),
                                     bracket = c(1e-10, 1),
                                     max_expand = 40) {
  
  stopifnot(is.data.frame(sim.dat),
            is.numeric(beta_X), length(beta_X) >= 1,
            is.numeric(omega), length(omega) == 1,
            is.numeric(target), target > 0, target < 1,
            is.numeric(shape), length(shape) == 1, shape > 0)
  
  psi_col <- if ("psi_true" %in% names(sim.dat)) "psi_true" else
    if ("psi" %in% names(sim.dat)) "psi" else
      stop("sim.dat must contain 'psi_true' or 'psi'.")
  
  sim_df <- tibble::as_tibble(sim.dat)
  Kdat <- dplyr::n_distinct(sim_df$exposure)
  if (length(beta_X) == 1L) beta_X <- rep(beta_X, Kdat)
  stopifnot(length(beta_X) == Kdat)
  
  beta_tbl <- data.frame(exposure = seq_len(Kdat), beta_x = as.numeric(beta_X))
  
  by_visit <- sim_df %>%
    dplyr::left_join(beta_tbl, by = "exposure") %>%
    dplyr::mutate(gpsi = gfun(.data[[psi_col]]),
                  term = beta_x * gpsi) %>%
    dplyr::group_by(id, visit) %>%
    dplyr::summarise(theta = omega * dplyr::first(W) + sum(term),
                     time_start = dplyr::first(time),
                     .groups = "drop") %>%
    dplyr::arrange(id, time_start) %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(time_stop = dplyr::lead(time_start)) %>%
    dplyr::filter(!is.na(time_stop) & time_stop > time_start) %>%
    dplyr::ungroup()
  
  Si <- by_visit %>%
    dplyr::mutate(risk = exp(theta) * (time_stop^shape - time_start^shape)) %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(S = sum(risk), .groups = "drop") %>%
    dplyr::pull(S)
  
  if (!all(is.finite(Si)) || all(Si <= .Machine$double.eps)) {
    stop("All subject-level cumulative risks are ~0; cannot calibrate Weibull lambda.")
  }
  
  f <- function(lambda) mean(1 - exp(-lambda * Si)) - target
  
  lo <- max(bracket[1], .Machine$double.eps)
  hi <- bracket[2]
  fhi <- f(hi)
  tries <- 0L
  while (fhi <= 0 && tries < max_expand) {
    hi <- hi * 2
    fhi <- f(hi)
    tries <- tries + 1L
  }
  if (fhi <= 0) stop("Could not bracket a solution for Weibull lambda.")
  
  uniroot(f, c(lo, hi), tol = 1e-10)$root
}



logp <- function(x) {
  ifelse(x>0, log(x), 0)
}

g_incr <- function(x, a = a_shift){
  log(pmax(x, 0) + a)
} 


make_pd <- function(S, eps = 1e-8) {
  ev <- eigen(S, symmetric = TRUE)
  if (min(ev$values) > eps) return(S)
  Q <- ev$vectors; d <- pmax(ev$values, eps)
  Q %*% diag(d, nrow = length(d)) %*% t(Q)
}

# ---- Data generation ----
# expected incement as the exposure in th survival model
simulate_joint_data <- function(N, K, visit_times,
                                # Fixed effects (k = 1..K); each is length-2: (intercept, slope)
                                gamma_list,       # occurrence on (1, Z)
                                delta_list,       # size on (1, Z)
                                sigma_eps,        # length K (SD on g=log scale)
                                # Random-effects covariance for (b_int, b_slope, c_int, c_slope) across k
                                Sigma,
                                # Per-exposure Normal covariate distribution at visits
                                z_means, z_sds,
                                # Baseline covariates (affect survival)
                                omega, weib_lambda, weib_shape, 
                                beta_X,       # beta_X length K
                                seed) {
  set.seed(seed)
  stopifnot(length(gamma_list) == K,
            length(delta_list) == K,
            length(sigma_eps) == K,
            length(z_means) == K,
            length(z_sds) == K,
            length(beta_X) == K,
            all(dim(Sigma) == c(4L * K, 4L * K)))
  
  tvec <- sort(unique(visit_times))
  stopifnot(length(tvec) >= 2L, all(diff(tvec) > 0))
  nJ <- length(tvec)
  
  # --- ensure Sigma is PD ---
  Sigma <- make_pd(Sigma)
  
  # --- subject-level random effects: (b_int, b_slope, c_int, c_slope) per k ---
  RE <- MASS::mvrnorm(n = N, mu = rep(0, 4 * K), Sigma = Sigma)
  colnames(RE) <- as.vector(t(outer(1:K, c("b_int","b_slope","c_int","c_slope"),
                                    function(k, nm) paste0("k", k, "_", nm))))
  
  # --- baseline covariate W ---
  W <- rnorm(N, 0, 1)
  
  # storage
  out_list <- vector("list", N)
  
  for (i in seq_len(N)) {
    # per-visit covariates: Z[j,k] at visit j 
    Z <- sapply(seq_len(K), function(k) rnorm(nJ, mean = z_means[k], sd = z_sds[k]))
    colnames(Z) <- paste0("Z", seq_len(K))
    Z <- as.matrix(Z)
    
    # split REs
    b_int   <- sapply(1:K, function(k) RE[i, paste0("k",k,"_b_int")])
    b_slope <- sapply(1:K, function(k) RE[i, paste0("k",k,"_b_slope")])
    c_int   <- sapply(1:K, function(k) RE[i, paste0("k",k,"_c_int")])
    c_slope <- sapply(1:K, function(k) RE[i, paste0("k",k,"_c_slope")])
    
    # storage for realized increments and sums
    X <- dX <- A <- R <- matrix(0, nrow = nJ, ncol = K)
    
    # --- simulate increments (occurrence then size) ---
    for (j in seq_len(nJ)) {
      if (j > 1L) {
        for (k in seq_len(K)) {
          eta_k <- (gamma_list[[k]][1] + b_int[k]) + (gamma_list[[k]][2] + b_slope[k]) * Z[j, k]
          A[j, k] <- rbinom(1, 1, plogis(eta_k))
          if (A[j, k] == 1L) {
            mu_k <- (delta_list[[k]][1] + c_int[k]) + (delta_list[[k]][2] + c_slope[k]) * Z[j, k]
            dX[j, k] <- exp(rnorm(1, mean = mu_k, sd = sigma_eps[k]))
          }
        }
      }
      if (j == 1L) {
        X[j, ] <- 0
      } else {
        X[j, ] <- X[j - 1L, ] + dX[j, ]
      }
      
      for (k in seq_len(K)) {
        R[j, k] <- as.integer(X[j, k] > 0)
      }
    }
    
    Xtilde_true <- matrix(0, nrow = nJ, ncol = K)
    Psi_true    <- matrix(0, nrow = nJ, ncol = K) 
    dmu_mat <- pi_mat <- matrix(0, nrow = nJ, ncol = K)
    
    for (j in 2:nJ) {
      for (k in 1:K) {
        eta_b <- (gamma_list[[k]][1] + b_int[k]) + (gamma_list[[k]][2] + b_slope[k]) * Z[j, k]
        pi_jk <- plogis(eta_b)
        mu_c  <- (delta_list[[k]][1] + c_int[k]) + (delta_list[[k]][2] + c_slope[k]) * Z[j, k]
        dmu   <- exp(mu_c + 0.5 * sigma_eps[k]^2)
        
        psi_jk <- pi_jk * dmu
        Psi_true[j, k]    <- psi_jk
        Xtilde_true[j, k] <- Xtilde_true[j-1, k] + psi_jk
        dmu_mat[j,k] <- dmu
        pi_mat[j,k] <- pi_jk
      }
    }
    
    # --- Survival on intervals: Weibull baseline h0(t) = lambda * nu * t^(nu - 1) ---
    T_event <- max(tvec)
    event <- 0L
    
    for (j in 1:(nJ - 1L)) {
      t_start <- tvec[j]
      t_end   <- tvec[j + 1L]
      
      eta_j <- sum(beta_X * logp(Psi_true[j, ])) + omega * W[i]
      Hinc  <- exp(eta_j) * weib_lambda * (t_end^weib_shape - t_start^weib_shape)
      p     <- 1 - exp(-Hinc)
      
      if (runif(1) < p) {
        u <- runif(1)
        Hdraw <- -log(1 - u * (1 - exp(-Hinc))) 
        T_raw <- (t_start^weib_shape + Hdraw / (weib_lambda * exp(eta_j)))^(1 / weib_shape)
        
        eps <- min(1e-6, 0.5 * (t_end - t_start))
        T_event <- max(t_start + eps, min(T_raw, t_end - eps))
        event <- 1L
        break
      }
    }
    
    # --- long format rows: stack visits × exposures ---
    df_i <- do.call(rbind, lapply(seq_len(K), function(k) {
      data.frame(
        id = i, exposure = k, visit = seq_len(nJ), time = tvec,
        Z    = Z[, k],           # end-of-interval value
        R = R[, k], A = A[, k], dX = dX[, k], X = X[, k],
        dmu = dmu_mat[, k],
        pi = pi_mat[, k],
        Xtilde = Xtilde_true[, k],
        psi    = Psi_true[, k],
        re_b_int   = b_int[k],   re_b_slope = b_slope[k],
        re_c_int   = c_int[k],   re_c_slope = c_slope[k],
        gamma0 = gamma_list[[k]][1], gamma1 = gamma_list[[k]][2],
        delta0 = delta_list[[k]][1], delta1 = delta_list[[k]][2],
        sigma_eps = sigma_eps[k],
        W = W[i], T_event = T_event, event = event
      )
    }))
    out_list[[i]] <- df_i
  }
  
  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  class(out) <- c("data.frame", "joint_sim")
  attr(out, "visit_times") <- tvec
  attr(out, "Sigma")       <- Sigma
  attr(out, "RE")          <- RE
  attr(out, "beta_X")      <- beta_X
  attr(out, "omega")      <- omega
  attr(out, "K")           <- K
  attr(out, "baseline")    <- list(type = "weibull", lambda = weib_lambda, shape = weib_shape)
  out
}


# ---------------- Utilities: indexing & selectors ----------------
re_index <- function(k, qb, qc) {
  pREk <- qb + qc
  a <- pREk*(k-1L)+1L
  list(
    b = a + seq_len(qb) - 1L,       # occurrence RE indices
    c = a + qb + seq_len(qc) - 1L   # size RE indices
  )
}

occ_row  <- function(Zjk, k, K, qb, qc) {
  pREk <- qb + qc
  pRE <- K * pREk
  idx  <- re_index(k, qb, qc)$b
  V    <- numeric(pRE)
  V[idx[1:2]] <- c(1, Zjk)  
  V
}


size_row <- function(Zbar_jk, k, K, qb, qc) {
  pREk <- qb + qc
  pRE <- K * pREk
  idx  <- re_index(k, qb, qc)$c
  V    <- numeric(pRE)
  V[idx[1:2]] <- c(1, Zbar_jk)  # assumes qc=2: (1, Zbar)
  V
}

logdet <- function(M) as.numeric(determinant(M, TRUE)$modulus)

symPD  <- function(S) {
  S <- 0.5*(S + t(S))
  ev <- eigen(S, symmetric = TRUE)
  ev$vectors %*% diag(pmax(ev$values, 1e-10), nrow(S)) %*% t(ev$vectors)
}



robust_newton_step <- function(H, s,
                               base_ridge = 1e-6, max_ridge = 1e8, max_tries = 8L,
                               eig_rel_tol = 1e-10, step_clip = Inf) {
  p <- nrow(H)
  ridge <- base_ridge
  for (t in 0:max_tries) {
    Ht <- (H + t(H)) * 0.5 + ridge * diag(p)
    cholHt <- tryCatch(chol(Ht), error = function(e) NULL)
    if (!is.null(cholHt)) {
      step <- backsolve(cholHt, forwardsolve(t(cholHt), s))
      n2 <- sqrt(sum(step^2))
      if (is.finite(step_clip) && is.finite(n2) && n2 > step_clip) step <- step * (step_clip / n2)
      return(list(step = as.numeric(step), ridge = ridge, method = "chol"))
    }
    ridge <- min(ridge * 10, max_ridge)
  }
  He <- (H + t(H)) * 0.5
  ee <- eigen(He, symmetric = TRUE)
  lam <- ee$values; U <- ee$vectors
  lam_max <- max(lam, 0)
  keep <- lam > max(eig_rel_tol * lam_max, .Machine$double.eps)
  if (!any(keep)) return(list(step = rep(0, p), ridge = NA_real_, method = "null"))
  step <- U[, keep, drop = FALSE] %*% ((t(U[, keep, drop = FALSE]) %*% s) / lam[keep])
  step <- as.numeric(step)
  n2 <- sqrt(sum(step^2))
  if (is.finite(step_clip) && is.finite(n2) && n2 > step_clip) step <- step * (step_clip / n2)
  list(step = step, ridge = NA_real_, method = "eigen")
}

# ---------------- Expected paths & gradients ----------------
build_Xtilde_paths <- function(df_i, m_i, S_i, K, qb, qc, gamma, delta, sigma_eps, tvec) {
  nJ <- length(tvec)
  Xtilde <- matrix(0, nJ, K)
  
  for (k in seq_len(K)) {
    for (j in seq_len(nJ)) {
      row_j <- which(df_i$visit == j & df_i$exposure == k)
      if (!length(row_j)) next
      Zj   <- df_i$Z[row_j]
      Zbj  <- df_i$Z[row_j]
      Vb   <- occ_row(Zj,  k, K, qb, qc)
      Vc   <- size_row(Zbj, k, K, qb, qc)
      
      a_fix_b <- sum(gamma[k,] * c(1, Zj))
      mu_b    <- a_fix_b + sum(Vb * m_i)
      # v_b     <- as.numeric(t(Vb) %*% S_i %*% Vb)
      # pi_hat  <- plogis(mu_b / sqrt(1 + (pi/8)*pmax(v_b,0)))
      pi_m    <- plogis(mu_b)
      
      a_fix_c <- sum(delta[k,] * c(1, Zbj))
      mu_c    <- a_fix_c + sum(Vc * m_i)
      v_c   <- as.numeric(t(Vc) %*% S_i %*% Vc)
      
      # Expos <- pi_m * exp(mu_c + 0.5*sigma_eps[k]^2 + 0.5 * v_c)
      Expos <- pi_m * exp(mu_c + 0.5*sigma_eps[k]^2)
      if (j > 1L) Xtilde[j, k] <- Xtilde[j-1, k] + Expos
    }
  }
  Xtilde
}




build_Psi_paths <- function(df_i, m_i, S_i, K, qb, qc, gamma, delta, sigma_eps, tvec) {
  nJ  <- length(tvec)
  Psi <- matrix(0, nJ, K)
  for (k in seq_len(K)) {
    for (j in 2:nJ) {  # first interval uses zero by convention
      row_j <- which(df_i$visit == j & df_i$exposure == k)
      if (!length(row_j)) next
      Zj <- df_i$Z[row_j]
      Zbj  <- df_i$Z[row_j]
      Vb <- occ_row(Zj,  k, K, qb, qc)
      Vc <- size_row(Zbj, k, K, qb, qc)
      
      mu_b <- sum(gamma[k, ] * c(1, Zj)) + sum(Vb * m_i)
      pi_m <- plogis(mu_b)
      
      mu_c <- sum(delta[k, ] * c(1, Zbj)) + sum(Vc * m_i)
      dmu  <- exp(mu_c + 0.5 * sigma_eps[k]^2)
      
      Psi[j, k] <- pi_m * dmu
    }
  }
  Psi
}



build_r_s_at_visits <- function(df_i, m_i, S_i, K, qb, qc,
                                gamma, delta, sigma_eps, tvec) {
  nJ  <- length(tvec)
  pRE <- K * (qb + qc)
  
  r_list <- vector("list", K)
  s_list <- vector("list", K)
  
  for (k in seq_len(K)) {
    r_k <- numeric(nJ)
    s_k <- matrix(0, nrow = nJ, ncol = pRE)
    
    for (j in seq_len(nJ)) {
      if (j == 1L) { r_k[j] <- 0; next }  # first interval: convention
      row_j <- which(df_i$visit == j & df_i$exposure == k)
      if (!length(row_j)) next
      
      Zj <- df_i$Z[row_j]
      Vb <- occ_row(Zj,  k, K, qb, qc)
      Vc <- size_row(Zj, k, K, qb, qc)
      
      mu_b <- sum(gamma[k, ] * c(1, Zj)) + sum(Vb * m_i)
      pi_m <- plogis(mu_b)
      
      mu_c <- sum(delta[k, ] * c(1, Zj)) + sum(Vc * m_i)
      dmu  <- exp(mu_c + 0.5 * sigma_eps[k]^2)
      
      psi_jk <- pi_m * dmu
      r_k[j] <- logp(psi_jk)                 # g(psi) = log psi with zero guard
      s_k[j, ] <- (1 - pi_m) * Vb + Vc       # ∂ log psi / ∂ d_i
    }
    
    r_list[[k]] <- r_k
    s_list[[k]] <- s_k
  }
  
  list(r = r_list, s = s_list)
}


# ---------------- Survival assembly (A_i, G_i, H_i) ----------------
surv_AGH <- function(df_i, m_i, S_i, K, qb, qc, beta, omega,
                     weib_lambda, weib_shape, tvec, r_s, inflate = FALSE) {
  nJ  <- length(tvec)
  pRE <- K * (qb + qc)
  
  a_mat <- matrix(0, nJ, pRE)
  alpha <- numeric(nJ)
  Wi    <- df_i$W[1]; Wi <- if (is.matrix(Wi) || is.data.frame(Wi)) as.numeric(Wi) else Wi
  
  for (j in seq_len(nJ)) {
    s_sum <- numeric(pRE); adj <- 0
    for (k in seq_len(K)) {
      s_kj <- r_s$s[[k]][j, ]; r_kj <- r_s$r[[k]][j]
      s_sum <- s_sum + beta[k] * s_kj
      adj   <- adj   + beta[k] * (r_kj - sum(s_kj * m_i))
    }
    a_mat[j, ] <- s_sum
    alpha[j]   <- sum(Wi * omega) + adj
  }
  
  Ti <- df_i$T_event[1]; Delta_i <- df_i$event[1]
  Ai <- 0; Gi <- numeric(pRE); Hi <- matrix(0, pRE, pRE)
  
  for (j in 1:(nJ-1)) {
    t0 <- tvec[j]; t1 <- tvec[j+1]; if (t0 >= Ti) break
    aj <- a_mat[j, ]; alphj <- alpha[j]
    
    t1i <- min(t1, Ti)
    base_inc <- baseline_cum_increment(t0, t1i, weib_lambda = weib_lambda, weib_shape = weib_shape)
    if (!is.finite(base_inc) || base_inc <= 0) next
    
    w_core <- base_inc * exp(alphj + sum(aj * m_i))
    w      <- if (inflate) w_core * exp(0.5 * as.numeric(t(aj) %*% S_i %*% aj)) else w_core
    
    Ai <- Ai + w
    Gi <- Gi + w * aj
    Hi <- Hi + w * tcrossprod(aj)
  }
  
  jstar <- max(which(tvec < Ti))
  a_Ti  <- if (is.finite(jstar) && jstar >= 1) a_mat[jstar, ] else a_mat[1, ]
  list(Ai = Ai, Gi = Gi, Hi = Hi, a_Ti = a_Ti, a_mat = a_mat, alpha = alpha)
}



# ---------------- One-subject Q,h (occ + size + survival) ----------------
subject_Qh <- function(df_i, m_i, S_i, K, qb, qc, gamma, delta, sigma_eps,
                       beta, omega, weib_lambda, weib_shape, tvec) {
  
  r_s <- build_r_s_at_visits(df_i, m_i, S_i, K, qb, qc, gamma, delta, sigma_eps, tvec)
  
  AGH <- surv_AGH(df_i, m_i, S_i, K, qb, qc, beta, omega,
                  weib_lambda, weib_shape, tvec, r_s, inflate = FALSE)
  
  pRE <- K * (qb + qc)
  Q_occ <- matrix(0, pRE, pRE)
  h_occ <- numeric(pRE)
  
  # Occurrence: JJ-bound quadratics
  for (k in seq_len(K)) {
    rows_k <- which(df_i$exposure == k & df_i$visit > 1L)
    for (row in rows_k) {
      Zj <- df_i$Z[row]
      Aij <- df_i$A[row]
      kappa <- Aij - 0.5
      Vb <- occ_row(Zj, k, K, qb, qc)
      a_fix <- sum(gamma[k,] * c(1, Zj))
      mu_d  <- sum(Vb * m_i)
      v_d   <- as.numeric(t(Vb) %*% S_i %*% Vb)
      xi    <- sqrt((a_fix + mu_d)^2 + v_d)
      lam   <- if (xi < 1e-10) 1/8 else tanh(xi/2)/(4*xi)
      Q_occ <- Q_occ + 2*lam * tcrossprod(Vb)
      h_occ <- h_occ + (kappa - 2*lam*a_fix) * Vb
    }
  }
  
  # Size: Gaussian on g-scale
  Q_siz <- matrix(0, pRE, pRE)
  h_siz <- numeric(pRE)
  for (k in seq_len(K)) {
    rows_k <- which(df_i$exposure == k & df_i$A == 1 & is.finite(df_i$dX) & df_i$dX > 0)
    if (!length(rows_k)) next
    for (row in rows_k) {
      Zbj <- df_i$Z[row]
      y <- logp(df_i$dX[row])
      Vc  <- size_row(Zbj, k, K, qb, qc)
      a_fix <- sum(delta[k,] * c(1, Zbj))
      w  <- 1/(sigma_eps[k]^2)
      Q_siz <- Q_siz + w * tcrossprod(Vc)
      h_siz <- h_siz + w * Vc * (y - a_fix)
    }
  }
  
  # Survival
  Q_surv <- AGH$Hi
  h_surv <- df_i$event[1] * AGH$a_Ti - AGH$Gi + Q_surv %*% m_i
  
  list(Q = Q_occ + Q_siz + Q_surv,
       h = h_occ + h_siz + h_surv,
       r_s = r_s,
       AGH = AGH)
}


# occ + size
subject_Qh_2stage <- function(df_i, m_i, S_i, K, qb, qc, gamma, delta, sigma_eps) {
  
  pRE <- K * (qb + qc)
  Q_occ <- matrix(0, pRE, pRE)
  h_occ <- numeric(pRE)
  
  # Occurrence: JJ-bound quadratics
  for (k in seq_len(K)) {
    rows_k <- which(df_i$exposure == k & df_i$visit > 1L)
    for (row in rows_k) {
      Zj <- df_i$Z[row]
      Aij <- df_i$A[row]
      kappa <- Aij - 0.5
      Vb <- occ_row(Zj, k, K, qb, qc)
      a_fix <- sum(gamma[k,] * c(1, Zj))
      mu_d  <- sum(Vb * m_i)
      v_d   <- as.numeric(t(Vb) %*% S_i %*% Vb)
      xi    <- sqrt((a_fix + mu_d)^2 + v_d)
      lam   <- if (xi < 1e-10) 1/8 else tanh(xi/2)/(4*xi)
      Q_occ <- Q_occ + 2*lam * tcrossprod(Vb)
      h_occ <- h_occ + (kappa - 2*lam*a_fix) * Vb
    }
  }
  
  # Size: Gaussian on g-scale
  Q_siz <- matrix(0, pRE, pRE)
  h_siz <- numeric(pRE)
  for (k in seq_len(K)) {
    rows_k <- which(df_i$exposure == k & df_i$A == 1 & is.finite(df_i$dX) & df_i$dX > 0)
    if (!length(rows_k)) next
    for (row in rows_k) {
      Zbj <- df_i$Z[row]
      y <- logp(df_i$dX[row])
      Vc  <- size_row(Zbj, k, K, qb, qc)
      a_fix <- sum(delta[k,] * c(1, Zbj))
      w  <- 1/(sigma_eps[k]^2)
      Q_siz <- Q_siz + w * tcrossprod(Vc)
      h_siz <- h_siz + w * Vc * (y - a_fix)
    }
  }
  
  list(Q = Q_occ + Q_siz,
       h = h_occ + h_siz)
}



# ---------------- Global updates for gamma and delta ----------------

update_gamma_JJ_WLS <- function(dat, ids, Ms, K, qb, qc, gamma) {
  # For each exposure k, build WLS system: w = 2λ(ξ), z = (A-1/2)/w - (V_b^T m_i)
  for (k in seq_len(K)) {
    Xfixed <- NULL; w_vec <- NULL; z_vec <- NULL
    for (ii in seq_along(ids)) {
      df_i <- dat[dat$id == ids[ii] & dat$exposure == k & dat$visit > 1L, , drop = FALSE]
      if (!nrow(df_i)) next
      mi <- Ms[[ii]]$m
      Si <- Ms[[ii]]$S
      for (row in seq_len(nrow(df_i))) {
        Zj <- df_i$Z[row]
        Aij <- df_i$A[row]
        Vb <- occ_row(Zj, k, K, qb, qc)
        a_fix <- sum(gamma[k,] * c(1, Zj))         # current fixed part
        mu_d  <- sum(Vb * mi)
        v_d   <- as.numeric(t(Vb) %*% Si %*% Vb)
        xi    <- sqrt((a_fix + mu_d)^2 + v_d)
        lam   <- if (xi < 1e-10) 1/8 else tanh(xi/2)/(4*xi)
        w     <- 2*lam
        z     <- (Aij - 0.5)/w - mu_d              # pseudo-response
        Xfixed <- rbind(Xfixed, c(1, Zj))
        w_vec  <- c(w_vec, w)
        z_vec  <- c(z_vec, z)
      }
    }
    if (!is.null(Xfixed)) {
      Wm  <- Matrix::Diagonal(x = pmax(w_vec, 1e-10))
      XtW <- t(Xfixed) %*% Wm
      lhs <- XtW %*% Xfixed
      rhs <- XtW %*% z_vec
      gamma[k, ] <- as.numeric(tryCatch(solve(lhs, rhs),
                                        error = function(e) MASS::ginv(as.matrix(lhs)) %*% rhs))
    }
  }
  gamma
}



update_delta_WLS_sigma <- function(dat, ids, Ms, K, qb, qc, delta, sigma_eps) {
  # For each exposure k, regress y* = log(dX) - Vc^T m_i on (1, Zbar) with weight 1/sigma^2
  for (k in seq_len(K)) {
    d_k <- dat[dat$exposure == k & dat$A == 1 & is.finite(dat$dX) & dat$dX > 0, , drop = FALSE]
    if (nrow(d_k) < 3) next
    
    y_star <- numeric(0)
    qterm  <- numeric(0)
    Xk     <- cbind(1, d_k$Z)
    
    for (ii in seq_along(ids)) {
      rows <- which(d_k$id == ids[ii])
      if (!length(rows)) next
      mi <- Ms[[ii]]$m; Si <- Ms[[ii]]$S
      for (t in rows) {
        Vc <- size_row(d_k$Z[t], k, K, qb, qc)
        y_star <- c(y_star, logp(d_k$dX[t]) - sum(Vc * mi))
        qterm  <- c(qterm,  as.numeric(t(Vc) %*% Si %*% Vc))
      }
    }
    
    if (length(y_star) >= 3) {
      w   <- rep(1 / (sigma_eps[k]^2), length(y_star))
      Wm  <- Matrix::Diagonal(x = w)
      XtW <- t(Xk) %*% Wm
      lhs <- XtW %*% Xk
      rhs <- XtW %*% y_star
      ck  <- tryCatch(solve(lhs, rhs), error = function(e) MASS::ginv(as.matrix(lhs)) %*% rhs)
      delta[k, ] <- as.numeric(ck)
      resid <- y_star - drop(Xk %*% ck)
      sigma_eps[k] <- sqrt(mean(pmax(resid^2 + qterm, 1e-12)))
    }
  }
  list(delta = delta, sigma_eps = sigma_eps)
}



# ---------------- Survival scores (omega, beta) ----------------
surv_scores <- function(dat, ids, Ms, K, qb, qc,
                        gamma, delta, sigma_eps,
                        beta, omega, weib_lambda, weib_shape, tvec,
                        # numeric stabilization
                        stabilize = FALSE, clip_eta_lin = 10, clip_asSa = 10) {
  
  Wi0 <- dat$W[match(ids[1], dat$id)][1]
  pW  <- if (is.matrix(Wi0) || is.data.frame(Wi0)) ncol(as.data.frame(Wi0)) else length(as.numeric(Wi0))
  
  score_beta  <- rep(0, K)
  H_beta      <- matrix(0, K, K)
  score_omega <- numeric(pW)
  H_omega     <- matrix(0, pW, pW)
  
  nJ <- length(tvec)
  
  for (ii in seq_along(ids)) {
    df_i <- dat[dat$id == ids[ii], , drop = FALSE]
    m_i  <- Ms[[ii]]$m
    S_i <- Ms[[ii]]$S
    
    r_s <- build_r_s_at_visits(df_i, m_i, S_i, K, qb, qc, gamma, delta, sigma_eps, tvec)
    
    Ti    <- df_i$T_event[1]
    Delta <- df_i$event[1]
    
    Wi <- df_i$W[1]
    Wi <- if (is.matrix(Wi) || is.data.frame(Wi)) as.numeric(Wi) else as.numeric(Wi)
    
    logw <- rep(-Inf, nJ - 1L)
    Phi  <- matrix(0, nrow = nJ - 1L, ncol = K)
    
    for (j in 1:(nJ - 1L)) {
      t0 <- tvec[j]; t1 <- tvec[j + 1L]
      if (t0 >= Ti) break
      
      t1i <- min(t1, Ti)
      
      base_inc <- baseline_cum_increment(t0, t1i, weib_lambda = weib_lambda, weib_shape = weib_shape)
      if (!is.finite(base_inc) || base_inc <= 0) next
      
      gX_j <- vapply(seq_len(K), function(k) r_s$r[[k]][j], 0.0)
      
      a_j <- numeric(length(m_i))
      for (k in seq_len(K)) a_j <- a_j + beta[k] * r_s$s[[k]][j, ]
      
      adj <- 0
      for (k in seq_len(K)) {
        s_kj <- r_s$s[[k]][j, ]
        adj  <- adj + beta[k] * (gX_j[k] - sum(s_kj * m_i))
      }
      alpha_j <- sum(Wi * omega) + adj
      
      eta_lin <- alpha_j + sum(a_j * m_i)
      aS_a    <- as.numeric(t(a_j) %*% S_i %*% a_j); if (!is.finite(aS_a)) aS_a <- 0
      aS_a    <- max(aS_a, 0)
      
      if (stabilize) {
        eta_lin <- max(min(eta_lin, clip_eta_lin), -clip_eta_lin)
        aS_a    <- min(aS_a, clip_asSa)
      }
      
      logw[j] <- log(base_inc) + eta_lin + 0.5 * aS_a
      
      Phi[j, ] <- vapply(seq_len(K), function(k) {
        s_kj <- r_s$s[[k]][j, ]
        gX_j[k] + as.numeric(t(a_j) %*% S_i %*% s_kj)
      }, 0.0)
    }
    
    finite_idx <- is.finite(logw)
    if (any(finite_idx)) {
      L <- max(logw[finite_idx])
      w_scaled <- exp(logw[finite_idx] - L)
      PhiF <- Phi[finite_idx, , drop = FALSE]
      
      H_i <- t(PhiF * sqrt(w_scaled)) %*% (PhiF * sqrt(w_scaled))
      s_i <- colSums(w_scaled * PhiF)
      
      H_beta     <- H_beta     + exp(L) * H_i
      score_beta <- score_beta - exp(L) * s_i
      
      wsum <- sum(w_scaled)
      H_omega     <- H_omega     + exp(L) * wsum * tcrossprod(Wi)
      score_omega <- score_omega - exp(L) * wsum * Wi
    }
    
    if (Delta == 1L) {
      jstar <- findInterval(Ti, tvec, left.open = TRUE)
      if (!is.finite(jstar) || jstar < 1) jstar <- 1
      gX_T <- vapply(seq_len(K), function(k) r_s$r[[k]][jstar], 0.0)
      score_beta  <- score_beta  + gX_T
      score_omega <- score_omega + Wi
    }
  }
  
  list(
    score_beta  = score_beta,
    H_beta      = H_beta,
    score_omega = score_omega,
    H_omega     = H_omega
  )
}



# ---------------- ELBO form ----------------
elbo_full <- function(dat, ids, Ms, K, qb, qc, gamma, delta, sigma_eps, Sigma,
                      beta, omega, weib_lambda, weib_shape, tvec) {
  N    <- length(ids)
  pRE  <- K * (qb + qc)
  qdim <- pRE
  
  Sigma_inv <- solve(Sigma)
  logdet_Sigma <- logdet(Sigma)
  
  total <- 0
  
  for (ii in seq_len(N)) {
    df_i <- dat[dat$id == ids[ii], , drop = FALSE]
    m_i  <- Ms[[ii]]$m; S_i <- Ms[[ii]]$S
    
    # ---------- (A) Occurrence: direct expectation (delta method) ----------
    elbo_occ <- 0
    for (k in seq_len(K)) {
      rows_k <- which(df_i$exposure == k & df_i$visit > 1L)
      if (!length(rows_k)) next
      for (row in rows_k) {
        Zj  <- df_i$Z[row]
        Aij <- df_i$A[row]
        kappa <- Aij - 0.5
        
        Vb <- occ_row(Zj, k, K, qb, qc)
        alpha <- sum(gamma[k, ] * c(1, Zj))
        mu_d  <- sum(Vb * m_i)
        v_d   <- as.numeric(t(Vb) %*% S_i %*% Vb)
        mu    <- alpha + mu_d
        
        xi  <- sqrt(mu^2 + v_d)
        lam <- if (xi < 1e-10) 1/8 else tanh(xi/2)/(4 * xi)
        
        quad_lin <- - lam * (mu_d^2 + v_d) + (kappa - 2 * lam * alpha) * mu_d
        c_occ <- kappa * alpha - lam * alpha^2 + (lam * xi^2 + 0.5 * xi - log1p(exp(xi)))
        elbo_occ <- elbo_occ + quad_lin + c_occ
      }
    }
    
    # ---------- (B) Size: Gaussian log-likelihood ----------
    elbo_size <- 0
    for (k in seq_len(K)) {
      rows_k <- which(df_i$exposure == k & df_i$A == 1 & is.finite(df_i$dX) & df_i$dX > 0)
      if (!length(rows_k)) next
      sig2 <- sigma_eps[k]^2
      for (row in rows_k) {
        y   <- log(df_i$dX[row])
        Zbj <- df_i$Z[row]
        Vc    <- size_row(Zbj, k, K, qb, qc)
        alpha <- sum(delta[k, ] * c(1, Zbj))
        mu_d  <- sum(Vc * m_i)
        v_d   <- as.numeric(t(Vc) %*% S_i %*% Vc)
        mu    <- alpha + mu_d
        elbo_size <- elbo_size + (-y) + (-0.5 * log(2 * pi * sig2)) - 0.5 * ((y - mu)^2 + v_d) / sig2
      }
    }
    
    # ---------- (C) Survival: quadratic approximation ----------
    r_s <- build_r_s_at_visits(df_i, m_i, S_i, K, qb, qc, gamma, delta, sigma_eps, tvec)
    AGH <- surv_AGH(df_i, m_i, S_i, K, qb, qc, beta, omega,
                    weib_lambda, weib_shape, tvec, r_s, inflate = TRUE)
    
    Ti     <- df_i$T_event[1]
    Delta  <- df_i$event[1]
    Ai_surv <- AGH$Ai
    elbo_surv <- -Ai_surv
    
    if (Delta == 1L) {
      jstar <- max(which(tvec < Ti)); if (!is.finite(jstar) || jstar < 1) jstar <- 1
      a_T   <- AGH$a_Ti
      alpha_T <- AGH$alpha[jstar]
      
      log_base <- baseline_loghaz(Ti, weib_lambda = weib_lambda, weib_shape = weib_shape)
      
      elbo_surv <- elbo_surv + log_base + alpha_T + sum(a_T * m_i)
    }
    
    # ---------- (D) Prior ----------
    prior_term <- -0.5 * qdim * log(2 * pi) - 0.5 * logdet_Sigma -
      0.5 * sum(Sigma_inv * S_i) - 0.5 * as.numeric(t(m_i) %*% Sigma_inv %*% m_i)
    
    # ---------- (E) Entropy ----------
    logdet_Si <- logdet(S_i)
    entropy_term <- 0.5 * logdet_Si + 0.5 * qdim * (1 + log(2 * pi))
    
    total <- total + elbo_occ + elbo_size + elbo_surv + prior_term + entropy_term
  }
  
  as.numeric(total)
}




elbo_2stage <- function(dat, ids, Ms, K, qb, qc, gamma, delta, sigma_eps, Sigma) {
  N    <- length(ids)
  pRE  <- K * (qb + qc)
  qdim <- pRE
  
  Sigma_inv <- solve(Sigma)
  logdet_Sigma <- logdet(Sigma)
  
  total <- 0
  
  for (ii in seq_len(N)) {
    df_i <- dat[dat$id == ids[ii], , drop = FALSE]
    m_i  <- Ms[[ii]]$m; S_i <- Ms[[ii]]$S
    
    # ---------- (A) Occurrence: direct expectation (delta method) ----------
    elbo_occ <- 0
    for (k in seq_len(K)) {
      rows_k <- which(df_i$exposure == k & df_i$visit > 1L)
      if (!length(rows_k)) next
      for (row in rows_k) {
        Zj  <- df_i$Z[row]
        Aij <- df_i$A[row]
        kappa <- Aij - 0.5
        
        Vb <- occ_row(Zj, k, K, qb, qc)
        alpha <- sum(gamma[k, ] * c(1, Zj))
        mu_d  <- sum(Vb * m_i)
        v_d   <- as.numeric(t(Vb) %*% S_i %*% Vb)
        mu    <- alpha + mu_d
        
        xi  <- sqrt(mu^2 + v_d)
        lam <- if (xi < 1e-10) 1/8 else tanh(xi/2)/(4 * xi)
        
        quad_lin <- - lam * (mu_d^2 + v_d) + (kappa - 2 * lam * alpha) * mu_d
        c_occ <- kappa * alpha - lam * alpha^2 + (lam * xi^2 + 0.5 * xi - log1p(exp(xi)))
        elbo_occ <- elbo_occ + quad_lin + c_occ
      }
    }
    
    # ---------- (B) Size: Gaussian log-likelihood ----------
    elbo_size <- 0
    for (k in seq_len(K)) {
      rows_k <- which(df_i$exposure == k & df_i$A == 1 & is.finite(df_i$dX) & df_i$dX > 0)
      if (!length(rows_k)) next
      sig2 <- sigma_eps[k]^2
      for (row in rows_k) {
        y   <- log(df_i$dX[row])
        Zbj <- df_i$Z[row]
        Vc    <- size_row(Zbj, k, K, qb, qc)
        alpha <- sum(delta[k, ] * c(1, Zbj))
        mu_d  <- sum(Vc * m_i)
        v_d   <- as.numeric(t(Vc) %*% S_i %*% Vc)
        mu    <- alpha + mu_d
        elbo_size <- elbo_size + (-y) + (-0.5 * log(2 * pi * sig2)) - 0.5 * ((y - mu)^2 + v_d) / sig2
      }
    }
    
    # ---------- (D) Prior ----------
    prior_term <- -0.5 * qdim * log(2 * pi) - 0.5 * logdet_Sigma -
      0.5 * sum(Sigma_inv * S_i) - 0.5 * as.numeric(t(m_i) %*% Sigma_inv %*% m_i)
    
    # ---------- (E) Entropy ----------
    logdet_Si <- logdet(S_i)
    entropy_term <- 0.5 * logdet_Si + 0.5 * qdim * (1 + log(2 * pi))
    
    total <- total + elbo_occ + elbo_size + prior_term + entropy_term
  }
  
  as.numeric(total)
}


fit_weibull_baseline_intervals <- function(tstart, tstop, status, risk_weight,
                                           init = c(lambda = 0.1, shape = 1),
                                           shape_bounds = c(0.05, 10)) {
  eps <- .Machine$double.eps
  
  stopifnot(length(tstart) == length(tstop),
            length(tstop) == length(status),
            length(status) == length(risk_weight),
            all(tstop > tstart),
            all(status %in% c(0, 1)),
            all(is.finite(risk_weight)),
            all(risk_weight >= 0))
  
  init_lambda <- if (is.null(names(init))) init[1] else init[["lambda"]]
  init_shape  <- if (is.null(names(init))) init[2] else init[["shape"]]
  
  D <- sum(status)
  if (D <= 0) {
    return(list(lambda = max(init_lambda, eps),
                shape = max(init_shape,  eps),
                vcov = matrix(NA_real_, 2, 2)))
  }
  
  event_times <- tstop[status == 1]
  sum_logT <- sum(log(pmax(event_times, eps)))
  
  prof_loglik <- function(log_shape) {
    shape <- exp(log_shape)
    inc <- pmax(tstop^shape - tstart^shape, 0)
    den <- sum(risk_weight * inc)
    if (!is.finite(den) || den <= 0) return(-Inf)
    lambda <- D / den
    D * log(lambda) + D * log(shape) + (shape - 1) * sum_logT - lambda * den
  }
  
  opt <- optimize(function(ls) -prof_loglik(ls), interval = log(shape_bounds))
  log_shape_hat <- opt$minimum
  shape_hat <- exp(log_shape_hat)
  
  inc_hat <- pmax(tstop^shape_hat - tstart^shape_hat, 0)
  den_hat <- sum(risk_weight * inc_hat)
  lambda_hat <- D / den_hat
  
  ll_logpar <- function(par) {
    log_lambda <- par[1]
    log_shape  <- par[2]
    lambda <- exp(log_lambda)
    shape  <- exp(log_shape)
    inc <- pmax(tstop^shape - tstart^shape, 0)
    D * log_lambda + D * log_shape + (shape - 1) * sum_logT -
      lambda * sum(risk_weight * inc)
  }
  
  H <- tryCatch(
    stats::optimHess(c(log(lambda_hat), log_shape_hat),
                     fn = function(par) -ll_logpar(par)),
    error = function(e) matrix(NA_real_, 2, 2)
  )
  
  vcov_log <- tryCatch(solve(H), error = function(e) matrix(NA_real_, 2, 2))
  J <- diag(c(lambda_hat, shape_hat), 2L)
  vcov <- J %*% vcov_log %*% J
  
  list(lambda = lambda_hat, shape = shape_hat, vcov = vcov)
}


update_weibull_baseline <- function(dat, ids, Ms, K, qb, qc,
                                    gamma, delta, sigma_eps,
                                    beta, omega, tvec,
                                    init = c(lambda = 0.1, shape = 1),
                                    shape_bounds = c(0.05, 10)) {
  tstart <- tstop <- status <- risk_weight <- numeric(0)
  
  for (ii in seq_along(ids)) {
    df_i <- dat[dat$id == ids[ii], , drop = FALSE]
    m_i  <- Ms[[ii]]$m
    S_i  <- Ms[[ii]]$S
    
    r_s <- build_r_s_at_visits(df_i, m_i, S_i, K, qb, qc, gamma, delta, sigma_eps, tvec)
    
    AGH <- surv_AGH(df_i, m_i, S_i, K, qb, qc, beta, omega,
                    weib_lambda = 1, weib_shape = 1,
                    tvec = tvec, r_s = r_s, inflate = TRUE)
    
    Ti <- df_i$T_event[1]
    Delta_i <- df_i$event[1]
    
    for (j in 1:(length(tvec) - 1L)) {
      tj0 <- tvec[j]
      if (tj0 >= Ti) break
      
      tj1 <- min(tvec[j + 1L], Ti)
      aj <- AGH$a_mat[j, ]
      alphj <- AGH$alpha[j]
      infl <- exp(0.5 * as.numeric(t(aj) %*% S_i %*% aj))
      rw <- exp(alphj + sum(aj * m_i)) * infl
      
      tstart <- c(tstart, tj0)
      tstop  <- c(tstop, tj1)
      risk_weight <- c(risk_weight, rw)
      status <- c(status, as.integer(Delta_i == 1L && Ti > tj0 && Ti <= tvec[j + 1L]))
    }
  }
  
  fit_weibull_baseline_intervals(
    tstart = tstart,
    tstop = tstop,
    status = status,
    risk_weight = risk_weight,
    init = init,
    shape_bounds = shape_bounds
  )
}


# ---------------- Main driver ----------------
fit_gva_joint <- function(dat, init_gamma, init_delta, init_sigma_eps, init_beta, init_omega, 
                          init_weib_lambda, init_weib_shape = 1,
                          max_iter = 200L, tol = 1e-6, verbose = TRUE, ridge = 1e-8,
                          weib_shape_bounds = c(0.05, 10)){
  ids  <- unique(dat$id)
  N    <- length(ids)
  K    <- attr(dat, "K"); if (is.null(K)) K <- length(unique(dat$exposure))
  tvec <- attr(dat, "visit_times"); stopifnot(length(tvec) >= 2)
  qb   <- attr(dat, "qb"); if (is.null(qb)) qb <- 2L
  qc   <- attr(dat, "qc"); if (is.null(qc)) qc <- 2L
  pRE  <- K * (qb + qc)
  
  # Globals (init)
  gamma <- init_gamma
  delta <- init_delta
  sigma_eps <- init_sigma_eps
  beta  <- init_beta
  omega <- init_omega
  
  Wi0   <- dat$W[match(ids[1], dat$id)][1]
  pW    <- if (is.matrix(Wi0) || is.data.frame(Wi0)) ncol(as.data.frame(Wi0)) else 1L
  
  
  weib_lambda <- init_weib_lambda
  weib_shape  <- init_weib_shape
  
  Sigma <- attr(dat, "Sigma"); if (is.null(Sigma)) Sigma <- diag(pRE)
  Ms <- lapply(seq_len(N), function(...) list(m = numeric(pRE), S = diag(pRE)))
  
  elbo_hist <- numeric(0)
  
  for (it in seq_len(max_iter)) {
    Sigma_inv <- tryCatch(solve(Sigma), error = function(e) MASS::ginv(Sigma))
    # ---- Local updates ----
    for (ii in seq_len(N)) {
      df_i <- dat[dat$id == ids[ii], , drop = FALSE]
      m_i  <- Ms[[ii]]$m; S_i <- Ms[[ii]]$S
      
      Qh <- subject_Qh(df_i, m_i, S_i, K, qb, qc, gamma, delta, sigma_eps,
                       beta, omega, weib_lambda, weib_shape, tvec)
      
      Qi <- Sigma_inv + Qh$Q
      Si_new <- tryCatch(solve(Qi), error = function(e) MASS::ginv(Qi))
      mi_new <- Si_new %*% Qh$h
      
      Ms[[ii]]$m <- as.numeric(mi_new)
      Ms[[ii]]$S <- as.matrix(symPD(Si_new))
    }
    
    # ---- Global Σ ----
    Sigma_new <- matrix(0, pRE, pRE)
    for (ii in seq_len(N)) {
      Sigma_new <- Sigma_new + (Ms[[ii]]$S + tcrossprod(Ms[[ii]]$m))
    }
    Sigma <- symPD(Sigma_new / N)
    
    # ---- Global γ (occurrence) ----
    gamma <- update_gamma_JJ_WLS(dat, ids, Ms, K, qb, qc, gamma)
    
    # ---- Global δ (size) and σ ----
    d_up <- update_delta_WLS_sigma(dat, ids, Ms, K, qb, qc, delta, sigma_eps)
    delta <- d_up$delta
    sigma_eps <- d_up$sigma_eps
    
    # ---- Baseline update ----
    wb_up <- update_weibull_baseline(
      dat = dat, ids = ids, Ms = Ms,
      K = K, qb = qb, qc = qc,
      gamma = gamma, delta = delta, sigma_eps = sigma_eps,
      beta = beta, omega = omega, tvec = tvec,
      init = c(lambda = weib_lambda, shape = weib_shape),
      shape_bounds = weib_shape_bounds
    )
    weib_lambda <- wb_up$lambda
    weib_shape  <- wb_up$shape
    
    # ---- (β, ω) Newton step ----
    sc <- surv_scores(dat, ids, Ms, K, qb, qc, gamma, delta, sigma_eps,
                      beta, omega, weib_lambda, weib_shape, tvec,
                      stabilize = FALSE, clip_eta_lin = 10, clip_asSa = 10)
    
    b_up <- robust_newton_step(sc$H_beta,  sc$score_beta,  base_ridge = 1e-6, step_clip = 2.0)
    o_up <- robust_newton_step(sc$H_omega, sc$score_omega, base_ridge = 1e-6, step_clip = 2.0)
    
    beta  <- as.numeric(beta  + b_up$step)
    omega <- as.numeric(omega + o_up$step)
    
    # ---- Monitor: full ELBO ----
    elb <- elbo_full(dat = dat, ids = ids, Ms = Ms, K = K, qb = qb, qc = qc,
                     gamma = gamma, delta = delta, sigma_eps = sigma_eps, Sigma = Sigma,
                     beta = beta, omega = omega,
                     weib_lambda = weib_lambda, weib_shape = weib_shape, tvec = tvec)
    elbo_hist <- c(elbo_hist, elb)
    
    # beta_str <- paste0("(", paste(round(beta, 4), collapse = ", "), ")")
    # if (verbose) {
    #   message(sprintf("%sIter %3d  ELBO=%.3f  lambda=%.4f  nu=%.4f  omega=%.4f  beta=%s",
    #               "", it, elb, weib_lambda, weib_shape, omega, beta_str))
    # }
    
    if (it >= 5) {
      rel_change <- abs((elbo_hist[it] - elbo_hist[it-1]) /(abs(elbo_hist[it-1]) + 1e-6))
      if (rel_change < tol) {
        # if (verbose) message("", "Converged (ELBO stabilized).")
        break
      }
    }
  }
  
  list(
    gamma = gamma,
    delta = delta,
    sigma_eps = sigma_eps,
    Sigma = Sigma,
    beta = beta,
    omega = omega,
    baseline = list(type = "weibull", lambda = weib_lambda, shape = weib_shape),
    locals = Ms,
    elbo = elbo_hist,
    settings = list(max_iter = max_iter, tol = tol)
  )
}


fit_gva_2stage <- function(dat, init_gamma, init_delta, init_sigma_eps,
                           max_iter = 200L, tol = 1e-6, verbose = TRUE, ridge = 1e-8){
  ids  <- unique(dat$id)
  N    <- length(ids)
  K    <- attr(dat, "K"); if (is.null(K)) K <- length(unique(dat$exposure))
  tvec <- attr(dat, "visit_times"); stopifnot(length(tvec) >= 2)
  qb   <- attr(dat, "qb"); if (is.null(qb)) qb <- 2L
  qc   <- attr(dat, "qc"); if (is.null(qc)) qc <- 2L
  pRE  <- K * (qb + qc)
  
  # Globals (init)
  gamma <- init_gamma
  delta <- init_delta
  sigma_eps <- init_sigma_eps
  
  Sigma <- attr(dat, "Sigma"); if (is.null(Sigma)) Sigma <- diag(pRE)
  Ms <- lapply(seq_len(N), function(...) list(m = numeric(pRE), S = diag(pRE)))
  
  elbo_hist <- numeric(0)
  
  for (it in seq_len(max_iter)) {
    # ---- Local updates ----
    for (ii in seq_len(N)) {
      df_i <- dat[dat$id == ids[ii], , drop = FALSE]
      m_i  <- Ms[[ii]]$m; S_i <- Ms[[ii]]$S
      
      Qh <- subject_Qh_2stage(df_i, m_i, S_i, K, qb, qc, gamma, delta, sigma_eps)
      
      
      Qi <- solve(Sigma) + Qh$Q
      Si_new <- tryCatch(solve(Qi), error = function(e) MASS::ginv(Qi))
      mi_new <- Si_new %*% Qh$h
      
      Ms[[ii]]$m <- as.numeric(mi_new)
      Ms[[ii]]$S <- as.matrix(symPD(Si_new))
    }
    
    # ---- Global Σ ----
    Sigma_new <- matrix(0, pRE, pRE)
    for (ii in seq_len(N)) {
      Sigma_new <- Sigma_new + (Ms[[ii]]$S + tcrossprod(Ms[[ii]]$m))
    }
    Sigma <- symPD(Sigma_new / N)
    
    # ---- Global γ (occurrence) ----
    gamma <- update_gamma_JJ_WLS(dat, ids, Ms, K, qb, qc, gamma)
    
    # ---- Global δ (size) and σ ----
    d_up <- update_delta_WLS_sigma(dat, ids, Ms, K, qb, qc, delta, sigma_eps)
    delta <- d_up$delta
    sigma_eps <- d_up$sigma_eps
    
    
    # ---- Monitor: full ELBO ----
    elb <- elbo_2stage(dat = dat, ids = ids, Ms = Ms, K = K, qb = qb, qc = qc,
                       gamma = gamma, delta = delta, sigma_eps = sigma_eps, Sigma = Sigma)
    elbo_hist <- c(elbo_hist, elb)
    
    if (verbose) {
      cat(sprintf("Iter %3d  ELBO=%.3f\n",
                  it, elb))
    }
    
    if (it >= 5) {
      rel_change <- abs((elbo_hist[it] - elbo_hist[it-1]) /(abs(elbo_hist[it-1]) + 1e-6))
      if (rel_change < tol) {
        if (verbose) cat("Converged (ELBO stabilized).\n")
        break
      }
    }
  }
  
  list(
    gamma = gamma,
    delta = delta,
    sigma_eps = sigma_eps,
    Sigma = Sigma,
    locals = Ms,
    elbo = elbo_hist,
    settings = list(max_iter = max_iter, tol = tol)
  )
}



# ---------- Helpers ----------
normalize_Kvec <- function(x, K, name) {
  if (length(x) == 1L) x <- rep(x, K)
  if (length(x) != K) {
    stop(sprintf("'%s' must have length 1 or K = %d.", name, K))
  }
  as.numeric(x)
}

normalize_Kx2 <- function(x, K, name) {
  if (is.list(x)) {
    if (length(x) != K) {
      stop(sprintf("'%s' list must have length K = %d.", name, K))
    }
    out <- do.call(rbind, lapply(seq_len(K), function(k) {
      xk <- as.numeric(x[[k]])
      if (length(xk) != 2L) {
        stop(sprintf("'%s[[%d]]' must have length 2.", name, k))
      }
      xk
    }))
    return(matrix(out, nrow = K, ncol = 2L))
  }
  
  if (is.vector(x) && length(x) == 2L) {
    return(matrix(rep(as.numeric(x), K), nrow = K, ncol = 2L, byrow = TRUE))
  }
  
  x <- as.matrix(x)
  if (all(dim(x) == c(K, 2L))) return(x)
  if (all(dim(x) == c(2L, K))) return(t(x))
  
  stop(sprintf("'%s' must be length 2, Kx2, 2xK, or a list of K length-2 vectors.", name))
}

block_idx <- function(k, block_size = 4L) {
  ((k - 1L) * block_size + 1L):(k * block_size)
}

coef_names2 <- function(prefix, K) {
  paste0(prefix, rep(seq_len(K), each = 2L), rep(0:1, times = K))
}

format_row <- function(x) {
  paste(format(x, scientific = FALSE, trim = TRUE, digits = 6), collapse = " ")
}

limit_blas_threads <- function(n = 1L) {
  n <- as.integer(n)
  
  # good to set these before forking
  Sys.setenv(
    OMP_NUM_THREADS        = as.character(n),
    OPENBLAS_NUM_THREADS   = as.character(n),
    MKL_NUM_THREADS        = as.character(n),
    VECLIB_MAXIMUM_THREADS = as.character(n),
    BLIS_NUM_THREADS       = as.character(n)
  )
  
  invisible(NULL)
}

# Simulation --------------------------------------------------------------
sim.inc.f <- function(sim.num, N, K, visit_times,
                      gamma1,          # length K or scalar
                      target_pi,       # length K or scalar
                      target_event,
                      delta,           # Kx2, 2xK, length-2, or list of K length-2 vectors
                      sigma_eps,       # length K or scalar
                      Sigma,
                      z_means, z_sds,  # length K or scalar
                      weib_shape,
                      omega,
                      beta_X,          # length K or scalar
                      cores.num, 
                      out_dir = "results") {
  
  gamma1    <- normalize_Kvec(gamma1, K, "gamma1")
  target_pi <- normalize_Kvec(target_pi, K, "target_pi")
  sigma_eps <- normalize_Kvec(sigma_eps, K, "sigma_eps")
  z_means   <- normalize_Kvec(z_means, K, "z_means")
  z_sds     <- normalize_Kvec(z_sds, K, "z_sds")
  beta_X    <- normalize_Kvec(beta_X, K, "beta_X")
  delta     <- normalize_Kx2(delta, K, "delta")
  
  stopifnot(is.matrix(Sigma), all(dim(Sigma) == c(4L * K, 4L * K)))
  
  # Solve gamma0 for each exposure
  gamma0 <- vapply(seq_len(K), function(k) {
    idx <- block_idx(k, 4L)
    Sigma_bk <- Sigma[idx, idx, drop = FALSE][1:2, 1:2, drop = FALSE]
    solve_gamma0_occurrence(
      gamma1 = gamma1[k],
      muZ = z_means[k],
      sdZ = z_sds[k],
      Sigma_b2x2 = Sigma_bk,
      target = target_pi[k]
    )
  }, numeric(1))
  
  gamma_mat  <- cbind(gamma0, gamma1)
  gamma_list <- lapply(seq_len(K), function(k) gamma_mat[k, ])
  delta_list <- lapply(seq_len(K), function(k) delta[k, ])
  
  # Calibrate Weibull lambda
  sim.dat.test <- simulate_joint_data(
    N = N, K = K, visit_times = visit_times,
    gamma_list = gamma_list,
    delta_list = delta_list,
    sigma_eps = sigma_eps,
    Sigma = Sigma,
    z_means = z_means,
    z_sds = z_sds,
    omega = omega,
    weib_lambda = 0.1,
    weib_shape = weib_shape,
    beta_X = beta_X,
    seed = 2L
  )
  
  weib_lambda <- calibrate_weibull_lambda(
    sim.dat = sim.dat.test,
    beta_X = beta_X,
    omega = omega,
    target = target_event,
    shape = weib_shape
  )
  
  
  # dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  file.out <- file.path(
    out_dir,
    paste0(
      "GVA_INC_NN", N,
      "_K", K,
      "_pi_", paste(format(target_pi, trim = TRUE), collapse = "-"),
      "_event_", format(target_event, trim = TRUE),
      "_shape_", format(weib_shape, trim = TRUE),
      "_omega_", format(omega, trim = TRUE),
      "_beta_", paste(format(round(beta_X,2), trim = TRUE), collapse = "-"),
      ".dat"
    )
  )
  
  if ( file.exists(as.character(file.out)) ) { unlink(as.character(file.out)) }
  
  header <- c(
    "nsim", "N", "K",
    paste0("target_pi", seq_len(K)),
    "target_event",
    coef_names2("true_gamma", K),
    coef_names2("true_delta", K),
    paste0("true_sigma_eps", seq_len(K)),
    "true_omega", "true_weib_shape", "true_weib_lambda",
    paste0("true_beta", seq_len(K)),
    coef_names2("est_gamma", K),
    coef_names2("est_delta", K),
    paste0("est_sigma_eps", seq_len(K)),
    "est_omega", "est_weib_shape", "est_weib_lambda",
    paste0("est_beta", seq_len(K))
  )
  
  cat(header, sep=" ", "\n", append=TRUE, file=file.out)
  
  one_sim <- function(nsim) {
    sim.dat <- simulate_joint_data(N = N, K = K, visit_times = visit_times,
                                   gamma_list = gamma_list, delta_list = delta_list,
                                   sigma_eps = sigma_eps, Sigma = Sigma,
                                   z_means = z_means, z_sds = z_sds,
                                   omega = omega, weib_lambda = weib_lambda, weib_shape = weib_shape, beta_X = beta_X,
                                   seed = nsim)
    
    # mean_pi_sim <- with(sim.dat[sim.dat$visit > 1, ], tapply(A, exposure, mean))
    # mean(sim.dat[sim.dat$visit == 1, "event"])
    # mean(sim.dat[sim.dat$visit == 1 & sim.dat$event == 1, "T_event"])
    
    fit <- try(fit_gva_joint(dat = sim.dat,
                             init_gamma = gamma_mat + matrix(runif(2 * K, 0, 0.2), nrow = K, ncol = 2L),
                             init_delta = delta + matrix(runif(2 * K, 0, 0.2), nrow = K, ncol = 2L),
                             init_sigma_eps = sigma_eps + runif(K, 0, 0.2),
                             init_beta = beta_X + runif(K, 0, 0.2),
                             init_omega = omega + runif(1, 0, 0.2),
                             init_weib_lambda = weib_lambda + runif(1, 0, 0.2),
                             init_weib_shape = weib_shape + runif(1, 0, 0.2),
                             max_iter = 200L, tol = 1e-4, verbose = TRUE, ridge = 1e-8),
               silent = TRUE)
    
    bad_fit <- inherits(fit, "try-error") ||
      is.null(fit$beta) ||
      !all(is.finite(c(
        as.vector(fit$gamma),
        as.vector(fit$delta),
        fit$sigma_eps,
        fit$omega,
        fit$baseline$shape,
        fit$baseline$lambda,
        fit$beta
      )))
    
    if (bad_fit) {
      message(sprintf("[attempt %03d] fit_gva_joint() failed; skipping.", nsim))
      return(NULL)
    }
    
    
    
    out_row <- c(
      nsim, N, K,
      target_pi,
      target_event,
      as.vector(t(gamma_mat)),
      as.vector(t(delta)),
      sigma_eps,
      omega, weib_shape, weib_lambda,
      beta_X,
      as.vector(t(fit$gamma)),
      as.vector(t(fit$delta)),
      fit$sigma_eps,
      fit$omega, fit$baseline$shape, fit$baseline$lambda,
      fit$beta
    )
    
    cat(out_row, sep=" ", "\n", append=TRUE, file=file.out)
    return()
  }
  
  if (cores.num > 1L) {
    limit_blas_threads(1L)
  }
  
  mclapply(1:sim.num, FUN = one_sim, mc.cores = cores.num)
  
  return()
  
}




sim.simple.f <- function(sim.num, N, K, visit_times, gamma11, gamma21, target_pi1, target_pi2, target_event, delta1, delta2, sigma_eps, Sigma, 
                         z_means, z_sds, omega, beta_X, a_shift, cores.num){
  Sigma1 <- Sigma[1:4, 1:4]
  Sigma2 <- Sigma[5:8, 5:8]
  Sigma_b1 <- Sigma1[1:2, 1:2]  # For exposure 1 (b1_int, b1_slope)
  Sigma_b2 <- Sigma2[1:2, 1:2] 
  # solve for gamma0 such that marginal probability = target_pi
  gamma10  <- solve_gamma0_occurrence(gamma1 = gamma11, muZ = z_means[1], sdZ = z_sds[1], Sigma_b2x2 = Sigma_b1, target = target_pi1)
  gamma20  <- solve_gamma0_occurrence(gamma1 = gamma21, muZ = z_means[2], sdZ = z_sds[2], Sigma_b2x2 = Sigma_b2, target = target_pi2)
  
  gamma1 <- c(gamma10, gamma11)
  gamma2 <- c(gamma20, gamma21)
  
  # solve h0 to reach target event rate
  sim.dat.test <- simulate_joint_data(N=N, K=K, visit_times=visit_times, gamma_list = list(gamma1, gamma2), delta_list = list(delta1,  delta2),
                                      sigma_eps = sigma_eps, Sigma = Sigma, z_means = z_means, z_sds   = z_sds,  
                                      omega = omega, h0 = 0.2, beta_X = beta_X, seed=2)
  
  h0 <- calibrate_h0(sim.dat = sim.dat.test, beta_X  = beta_X, omega  = omega, target  = target_event,     # aim for ~80% events on average
                     gfun = logp)
  
  
  # file output
  file.out <- paste("results/simple/", "SIMPLE", "_a", a_shift, "_NN", N, "_K",  K, "_gamma11", round(gamma11, 2), "_gamma21", round(gamma21, 2), 
                    "_targetpi1", target_pi1, "_targetpi2", target_pi2, "_targetevent", target_event, "_omega", omega, 
                    "_beta1_", round(beta_X[1], 2), "_beta2_", round(beta_X[2], 2), ".dat", sep="")
  
  if ( file.exists(as.character(file.out)) ) { unlink(as.character(file.out)) }
  
  cat("nsim", "N", "K", "target_pi1", "target_pi2", 
      "true_gamma10", "true_gamma11", "true_gamma20", "true_gamma21", "true_delta10", "true_delta11", "true_delta20", "true_delta21", 
      "true_omega", "true_h0", "true_beta1", "true_beta2",
      "est_omega",  "est_beta1", "est_beta2", "est_h0",
      "se_omega",  "se_beta1", "se_beta2", "se_h0",
      sep=" ", "\n", append=TRUE, file=file.out)
  
  parallel.sim.f <- function(nsim, N, K, visit_times, gamma1, gamma2, target_pi1, target_pi2, target_event, delta1, delta2, sigma_eps, Sigma, 
                             z_means, z_sds, omega, h0, beta_X){
    
    # --- simulate (use a distinct seed so each replicate differs) ---
    sim.dat <- simulate_joint_data(N=N, K=K, visit_times=visit_times,
                                   gamma_list = list(gamma1, gamma2),
                                   delta_list = list(delta1,  delta2),
                                   sigma_eps = sigma_eps, Sigma = Sigma,
                                   z_means = z_means, z_sds   = z_sds,   
                                   omega = omega, h0 = h0, beta_X = beta_X,
                                   seed = nsim)
    
    sim.dat <- data.frame(sim.dat)
    # Subject-level event/censoring info (unique per id)
    subj <- sim.dat %>% group_by(id) %>% summarise(T_event = first(T_event),
                                                   event   = first(event),
                                                   W       = first(W),
                                                   .groups = "drop")
    
    # Wide dX at each visit time: dXk1, dXk2, ..., dXkK
    dx_wide <- sim.dat %>%
      dplyr::select(id, time, exposure, dX) %>%
      distinct() %>%
      mutate(exposure = as.integer(exposure)) %>%
      tidyr::pivot_wider(
        names_from  = exposure,
        values_from = dX,
        names_prefix = "dXk",
        values_fill = 0
      ) %>%
      arrange(id, time)
    
    # Merge in W and event time
    tv0 <- dx_wide %>% left_join(subj, by = "id") %>% arrange(id, time)
    
    # Convert visits to start/stop intervals and place event at T_event
    # Covariates are carried forward from the visit at tstart.
    
    tv <- tv0 %>%
      group_by(id) %>%
      arrange(time, .by_group = TRUE) %>%
      mutate(
        tstart = time,
        tstop_nominal = lead(time),
        event_time = first(T_event),
        event_ind  = first(event)
      ) %>%
      filter(!is.na(tstop_nominal)) %>%          # drop last visit (no next interval)
      filter(tstart < event_time) %>%            # drop intervals starting after event/censor
      mutate(
        # if event happens before the next planned visit, truncate interval end at event time
        tstop  = pmin(tstop_nominal, event_time),
        status = as.integer(event_ind == 1 & event_time > tstart & event_time <= tstop_nominal)
      ) %>%
      filter(tstop > tstart) %>%
      ungroup()
    
    
    direct.fit <- coxph(Surv(tstart, tstop, status) ~ W + g_incr(dXk1, a=a_shift) + g_incr(dXk2, a=a_shift) + cluster(id), data = tv, ties = "efron")
    # print(summary(direct.fit))
    
    
    lp <- predict(direct.fit, newdata = tv, type = "lp") # Linear predictor for each interval-row
    dt <- tv$tstop - tv$tstart # Interval lengths
    D <- sum(tv$status) # Total number of events
    den <- sum(dt * exp(lp)) # Denominator: total weighted time at risk
    
    # Constant baseline hazard estimate (per unit time of 'time' variable)
    lambda_hat <- D / den
    se_lambda <- ifelse(D > 0, sqrt(D) / den, NA)
    
    cat(nsim, N, K, target_pi1, target_pi2, 
        gamma1[1], gamma1[2], gamma2[1], gamma2[2], delta1[1], delta1[2], delta2[1], delta2[2], 
        omega, h0, beta_X[1], beta_X[2], 
        direct.fit$coefficients, lambda_hat, 
        sqrt(diag(direct.fit$var)), se_lambda,
        sep=" ", "\n", append=TRUE, file=file.out)
    
    return()
    
  }
  
  mclapply(1:sim.num, parallel.sim.f, N, K, visit_times, gamma1, gamma2, target_pi1, target_pi2, target_event, delta1, delta2, sigma_eps, Sigma, 
           z_means, z_sds, omega, h0, beta_X, mc.cores=cores.num)
  
  
  return()
}



# with true expected exposure
sim.simple.true.f <- function(sim.num, N, K, visit_times, gamma11, gamma21, target_pi1, target_pi2, target_event, delta1, delta2, sigma_eps, Sigma, 
                              z_means, z_sds, omega, beta_X, cores.num){
  Sigma1 <- Sigma[1:4, 1:4]
  Sigma2 <- Sigma[5:8, 5:8]
  Sigma_b1 <- Sigma1[1:2, 1:2]  # For exposure 1 (b1_int, b1_slope)
  Sigma_b2 <- Sigma2[1:2, 1:2] 
  # solve for gamma0 such that marginal probability = target_pi
  gamma10  <- solve_gamma0_occurrence(gamma1 = gamma11, muZ = z_means[1], sdZ = z_sds[1], Sigma_b2x2 = Sigma_b1, target = target_pi1)
  gamma20  <- solve_gamma0_occurrence(gamma1 = gamma21, muZ = z_means[2], sdZ = z_sds[2], Sigma_b2x2 = Sigma_b2, target = target_pi2)
  
  gamma1 <- c(gamma10, gamma11)
  gamma2 <- c(gamma20, gamma21)
  
  # solve h0 to reach target event rate
  sim.dat.test <- simulate_joint_data(N=N, K=K, visit_times=visit_times, gamma_list = list(gamma1, gamma2), delta_list = list(delta1,  delta2),
                                      sigma_eps = sigma_eps, Sigma = Sigma, z_means = z_means, z_sds   = z_sds,  
                                      omega = omega, h0 = 0.2, beta_X = beta_X, seed=2)
  
  h0 <- calibrate_h0(sim.dat = sim.dat.test, beta_X  = beta_X, omega  = omega, target  = target_event,     # aim for ~80% events on average
                     gfun = logp)
  
  
  # file output
  file.out <- paste("results/simple1/", "SIMPLE", "_NN", N, "_K",  K, "_gamma11", round(gamma11, 2), "_gamma21", round(gamma21, 2), 
                    "_targetpi1", target_pi1, "_targetpi2", target_pi2, "_targetevent", target_event, "_omega", omega, 
                    "_beta1_", round(beta_X[1], 2), "_beta2_", round(beta_X[2], 2), ".dat", sep="")
  
  if ( file.exists(as.character(file.out)) ) { unlink(as.character(file.out)) }
  
  cat("nsim", "N", "K", "target_pi1", "target_pi2", 
      "true_gamma10", "true_gamma11", "true_gamma20", "true_gamma21", "true_delta10", "true_delta11", "true_delta20", "true_delta21", 
      "true_omega", "true_h0", "true_beta1", "true_beta2",
      "est_omega",  "est_beta1", "est_beta2", "est_h0",
      "se_omega",  "se_beta1", "se_beta2", "se_h0",
      sep=" ", "\n", append=TRUE, file=file.out)
  
  parallel.sim.f <- function(nsim, N, K, visit_times, gamma1, gamma2, target_pi1, target_pi2, target_event, delta1, delta2, sigma_eps, Sigma, 
                             z_means, z_sds, omega, h0, beta_X){
    
    # --- simulate (use a distinct seed so each replicate differs) ---
    sim.dat <- simulate_joint_data(N=N, K=K, visit_times=visit_times,
                                   gamma_list = list(gamma1, gamma2),
                                   delta_list = list(delta1,  delta2),
                                   sigma_eps = sigma_eps, Sigma = Sigma,
                                   z_means = z_means, z_sds   = z_sds,   
                                   omega = omega, h0 = h0, beta_X = beta_X,
                                   seed = nsim)
    
    sim.dat <- data.frame(sim.dat)
    # Subject-level event/censoring info (unique per id)
    subj <- sim.dat %>% group_by(id) %>% summarise(T_event = first(T_event),
                                                   event   = first(event),
                                                   W       = first(W),
                                                   .groups = "drop")
    
    # Wide dX at each visit time: dXk1, dXk2, ..., dXkK
    psi_wide <- sim.dat %>%
      dplyr::select(id, time, exposure, psi) %>%
      distinct() %>%
      mutate(exposure = as.integer(exposure)) %>%
      tidyr::pivot_wider(
        names_from  = exposure,
        values_from = psi,
        names_prefix = "psi",
        values_fill = 0
      ) %>%
      arrange(id, time)
    
    # Merge in W and event time
    tv0 <- psi_wide %>% left_join(subj, by = "id") %>% arrange(id, time)
    
    # Convert visits to start/stop intervals and place event at T_event
    # Covariates are carried forward from the visit at tstart.
    
    tv <- tv0 %>%
      group_by(id) %>%
      arrange(time, .by_group = TRUE) %>%
      mutate(
        tstart = time,
        tstop_nominal = lead(time),
        event_time = first(T_event),
        event_ind  = first(event)
      ) %>%
      filter(!is.na(tstop_nominal)) %>%          # drop last visit (no next interval)
      filter(tstart < event_time) %>%            # drop intervals starting after event/censor
      mutate(
        # if event happens before the next planned visit, truncate interval end at event time
        tstop  = pmin(tstop_nominal, event_time),
        status = as.integer(event_ind == 1 & event_time > tstart & event_time <= tstop_nominal)
      ) %>%
      filter(tstop > tstart) %>%
      ungroup()
    
    
    direct.fit <- coxph(Surv(tstart, tstop, status) ~ W + logp(psi1) + logp(psi2) + cluster(id), data = tv, ties = "efron")
    # print(summary(direct.fit))
    
    
    lp <- predict(direct.fit, newdata = tv, type = "lp") # Linear predictor for each interval-row
    dt <- tv$tstop - tv$tstart # Interval lengths
    D <- sum(tv$status) # Total number of events
    den <- sum(dt * exp(lp)) # Denominator: total weighted time at risk
    
    # Constant baseline hazard estimate (per unit time of 'time' variable)
    lambda_hat <- D / den
    se_lambda <- ifelse(D > 0, sqrt(D) / den, NA)
    
    cat(nsim, N, K, target_pi1, target_pi2, 
        gamma1[1], gamma1[2], gamma2[1], gamma2[2], delta1[1], delta1[2], delta2[1], delta2[2], 
        omega, h0, beta_X[1], beta_X[2], 
        direct.fit$coefficients, lambda_hat, 
        sqrt(diag(direct.fit$var)), se_lambda,
        sep=" ", "\n", append=TRUE, file=file.out)
    
    return()
    
  }
  
  mclapply(1:sim.num, parallel.sim.f, N, K, visit_times, gamma1, gamma2, target_pi1, target_pi2, target_event, delta1, delta2, sigma_eps, Sigma, 
           z_means, z_sds, omega, h0, beta_X, mc.cores=cores.num)
  
  
  return()
}




sim.2stage.f <- function(sim.num, N, K, visit_times, gamma11, gamma21, target_pi1, target_pi2, target_event, delta1, delta2, sigma_eps, Sigma, 
                         z_means, z_sds, omega, beta_X, cores.num){
  Sigma1 <- Sigma[1:4, 1:4]
  Sigma2 <- Sigma[5:8, 5:8]
  Sigma_b1 <- Sigma1[1:2, 1:2]  # For exposure 1 (b1_int, b1_slope)
  Sigma_b2 <- Sigma2[1:2, 1:2] 
  # solve for gamma0 such that marginal probability = target_pi
  gamma10  <- solve_gamma0_occurrence(gamma1 = gamma11, muZ = z_means[1], sdZ = z_sds[1], Sigma_b2x2 = Sigma_b1, target = target_pi1)
  gamma20  <- solve_gamma0_occurrence(gamma1 = gamma21, muZ = z_means[2], sdZ = z_sds[2], Sigma_b2x2 = Sigma_b2, target = target_pi2)
  
  gamma1 <- c(gamma10, gamma11)
  gamma2 <- c(gamma20, gamma21)
  
  # solve h0 to reach target event rate
  sim.dat.test <- simulate_joint_data(N=N, K=K, visit_times=visit_times, gamma_list = list(gamma1, gamma2), delta_list = list(delta1,  delta2),
                                      sigma_eps = sigma_eps, Sigma = Sigma, z_means = z_means, z_sds   = z_sds,  
                                      omega = omega, h0 = 0.2, beta_X = beta_X, seed=2)
  
  h0 <- calibrate_h0(sim.dat = sim.dat.test, beta_X  = beta_X, omega  = omega, target  = target_event,     # aim for ~80% events on average
                     gfun = logp)
  
  
  # file output
  file.out <- paste("results/simple/", "TWOSTAGE", "_NN", N, "_K",  K, "_gamma11", round(gamma11, 2), "_gamma21", round(gamma21, 2), 
                    "_targetpi1", target_pi1, "_targetpi2", target_pi2, "_targetevent", target_event, "_omega", omega, 
                    "_beta1_", round(beta_X[1], 2), "_beta2_", round(beta_X[2], 2), ".dat", sep="")
  
  if ( file.exists(as.character(file.out)) ) { unlink(as.character(file.out)) }
  
  cat("nsim", "N", "K", "target_pi1", "target_pi2", 
      "true_gamma10", "true_gamma11", "true_gamma20", "true_gamma21", "true_delta10", "true_delta11", "true_delta20", "true_delta21", 
      "true_omega", "true_h0", "true_beta1", "true_beta2",
      "est_gamma10"," est_gamma11", "est_gamma20", "est_gamma21", "est_delta10", "est_delta11",  "est_delta20", "est_delta21", 
      "est_sigma_eps1", "est_sigma_eps2",
      "est_omega",  "est_beta1", "est_beta2", "est_h0",
      "se_omega",  "se_beta1", "se_beta2", "se_h0", 
      sep=" ", "\n", append=TRUE, file=file.out)
  
  parallel.sim.f <- function(nsim, N, K, visit_times, gamma1, gamma2, target_pi1, target_pi2, target_event, delta1, delta2, sigma_eps, Sigma, 
                             z_means, z_sds, omega, h0, beta_X){
    
    # --- simulate (use a distinct seed so each replicate differs) ---
    sim.dat <- simulate_joint_data(N=N, K=K, visit_times=visit_times,
                                   gamma_list = list(gamma1, gamma2),
                                   delta_list = list(delta1,  delta2),
                                   sigma_eps = sigma_eps, Sigma = Sigma,
                                   z_means = z_means, z_sds   = z_sds,   
                                   omega = omega, h0 = h0, beta_X = beta_X,
                                   seed = nsim)
    
    # --- fit ---
    fit <- try(fit_gva_2stage(dat=sim.dat, init_gamma=matrix(c(gamma1, gamma2)+runif(1,0, 0.2), nrow = K, byrow = T), 
                              init_delta=matrix(c(delta1, delta2)+runif(1,0, 0.2), nrow = K, byrow = T), 
                              init_sigma_eps=sigma_eps+runif(1,0, 0.2),
                              max_iter = 200L, tol = 1e-4, verbose = TRUE, ridge = 1e-8),
               silent = TRUE)
    
    sim.dat$pi_hat <- NA_real_
    sim.dat$dmu_hat <- NA_real_
    
    for (r in seq_len(nrow(sim.dat))) {
      id <- sim.dat$id[r]
      k  <- sim.dat$exposure[r]
      
      m_i <- fit$locals[[id]]$m
      b_ik <- m_i[re_index(k, 2, 2)$b]
      c_ik <- m_i[re_index(k, 2, 2)$c]
      
      row <- sim.dat[r, , drop = FALSE]
      Z <- row$Z
      
      # occurrence prob
      eta_pi <- drop(c(1, Z) %*% fit$gamma[k,] + c(1, Z) %*% b_ik)
      sim.dat$pi_hat[r] <- expit(eta_pi)
      
      # positive increment mean (plug-in posterior mean for c_ik)
      eta_mu <- drop(c(1, Z) %*% fit$delta[k,] + c(1, Z) %*% c_ik)
      sim.dat$dmu_hat[r] <- exp(eta_mu + 0.5 * fit$sigma_eps[k])
    }
    
    
    
    
    sim.dat <- data.frame(sim.dat)
    sim.dat0 <- sim.dat %>% mutate(psi_hat = pi_hat*dmu_hat)
    
    # Subject-level event/censoring info (unique per id)
    subj <- sim.dat0 %>% group_by(id) %>% summarise(T_event = first(T_event),
                                                    event   = first(event),
                                                    W       = first(W),
                                                    .groups = "drop")
    
    # Wide psi_hat at each visit time
    dx_wide <- sim.dat0 %>%
      dplyr::select(id, time, exposure, psi_hat) %>%
      distinct() %>%
      mutate(exposure = as.integer(exposure)) %>%
      tidyr::pivot_wider(
        names_from  = exposure,
        values_from = psi_hat,
        names_prefix = "psi_hatk",
        values_fill = 0
      ) %>%
      arrange(id, time)
    
    # Merge in W and event time
    tv0 <- dx_wide %>% left_join(subj, by = "id") %>% arrange(id, time)
    
    # Convert visits to start/stop intervals and place event at T_event
    # Covariates are carried forward from the visit at tstart.
    
    tv <- tv0 %>%
      group_by(id) %>%
      arrange(time, .by_group = TRUE) %>%
      mutate(
        tstart = time,
        tstop_nominal = lead(time),
        event_time = first(T_event),
        event_ind  = first(event)
      ) %>%
      filter(!is.na(tstop_nominal)) %>%          # drop last visit (no next interval)
      filter(tstart < event_time) %>%            # drop intervals starting after event/censor
      mutate(
        # if event happens before the next planned visit, truncate interval end at event time
        tstop  = pmin(tstop_nominal, event_time),
        status = as.integer(event_ind == 1 & event_time > tstart & event_time <= tstop_nominal)
      ) %>%
      filter(tstop > tstart) %>%
      ungroup()
    
    
    direct.fit <- coxph(Surv(tstart, tstop, status) ~ W + log(psi_hatk1) + log(psi_hatk2) + cluster(id), data = tv, ties = "efron")
    # print(summary(direct.fit))
    
    
    lp <- predict(direct.fit, newdata = tv, type = "lp") # Linear predictor for each interval-row
    dt <- tv$tstop - tv$tstart # Interval lengths
    D <- sum(tv$status) # Total number of events
    den <- sum(dt * exp(lp)) # Denominator: total weighted time at risk
    
    # Constant baseline hazard estimate (per unit time of 'time' variable)
    lambda_hat <- D / den
    se_lambda <- ifelse(D > 0, sqrt(D) / den, NA)
    
    
    
    
    cat(nsim, N, K, target_pi1, target_pi2, 
        gamma1[1], gamma1[2], gamma2[1], gamma2[2], delta1[1], delta1[2], delta2[1], delta2[2], 
        omega, h0, beta_X[1], beta_X[2], 
        fit$gamma[1, ], fit$gamma[2, ], fit$delta[1, ], fit$delta[2, ],
        fit$sigma_eps,
        direct.fit$coefficients, lambda_hat, 
        sqrt(diag(direct.fit$var)), se_lambda,
        sep=" ", "\n", append=TRUE, file=file.out)
    
    return()
    
  }
  
  
  mclapply(1:sim.num, parallel.sim.f, N, K, visit_times, gamma1, gamma2, target_pi1, target_pi2, target_event, delta1, delta2, sigma_eps, Sigma, 
           z_means, z_sds, omega, h0, beta_X, mc.cores=cores.num)
  
  
  return()
}




dat.plot.f <- function(seed, N, K, visit_times, gamma11, gamma21, target_pi1, target_pi2, target_event, delta1, delta2, sigma_eps, Sigma, 
                       z_means, z_sds, omega, beta_X){
  Sigma1 <- Sigma[1:4, 1:4]
  Sigma2 <- Sigma[5:8, 5:8]
  Sigma_b1 <- Sigma1[1:2, 1:2]  # For exposure 1 (b1_int, b1_slope)
  Sigma_b2 <- Sigma2[1:2, 1:2] 
  # solve for gamma0 such that marginal probability = target_pi
  gamma10  <- solve_gamma0_occurrence(gamma1 = gamma11, muZ = z_means[1], sdZ = z_sds[1], Sigma_b2x2 = Sigma_b1, target = target_pi1)
  gamma20  <- solve_gamma0_occurrence(gamma1 = gamma21, muZ = z_means[2], sdZ = z_sds[2], Sigma_b2x2 = Sigma_b2, target = target_pi2)
  
  gamma1 <- c(gamma10, gamma11)
  gamma2 <- c(gamma20, gamma21)
  
  # solve h0 to reach target event rate
  sim.dat.test <- simulate_joint_data(N=N, K=K, visit_times=visit_times, gamma_list = list(gamma1, gamma2), delta_list = list(delta1,  delta2),
                                      sigma_eps = sigma_eps, Sigma = Sigma, z_means = z_means, z_sds   = z_sds,  
                                      omega = omega, h0 = 0.2, beta_X = beta_X, seed=2)
  
  h0 <- calibrate_h0(sim.dat = sim.dat.test, beta_X  = beta_X, omega  = omega, target  = target_event,     # aim for ~80% events on average
                     gfun = logp)
  
  
  # --- simulate (use a distinct seed so each replicate differs) ---
  sim.dat <- simulate_joint_data(N=N, K=K, visit_times=visit_times,
                                 gamma_list = list(gamma1, gamma2),
                                 delta_list = list(delta1,  delta2),
                                 sigma_eps = sigma_eps, Sigma = Sigma,
                                 z_means = z_means, z_sds   = z_sds,   
                                 omega = omega, h0 = h0, beta_X = beta_X,
                                 seed = seed)
  
  sim.dat <- data.frame(sim.dat)
  dat.pre <- sim.dat %>% filter(exposure == 1) %>% filter(visit >1)
  
  return(dat.pre)
}



# Summary -----------------------------------------------------------------
fmt3 <- function(x) {
  format(round(as.numeric(x), 3), nsmall = 3, trim = TRUE)
}

collapse_num <- function(x, digits = 2L) {
  x <- as.numeric(x)
  paste(format(round(x, digits), nsmall = digits, trim = TRUE), collapse = ",")
}


gva_param_bases <- function(K) {
  list(
    beta      = paste0("beta", seq_len(K)),
    gamma     = coef_names2("gamma", K),
    delta     = coef_names2("delta", K),
    sigma_eps = paste0("sigma_eps", seq_len(K)),
    scalar    = c("omega", "weib_shape", "weib_lambda")
  )
}

get_true_sigma_eps <- function(sim.out, K, sigma_eps = NULL) {
  true_cols <- paste0("true_sigma_eps", seq_len(K))
  if (all(true_cols %in% names(sim.out))) {
    return(vapply(true_cols, function(nm) mean(sim.out[[nm]], na.rm = TRUE), numeric(1)))
  }
  if (is.null(sigma_eps)) {
    stop("true_sigma_eps columns not found; please supply 'sigma_eps'.")
  }
  normalize_Kvec(sigma_eps, K, "sigma_eps")
}

mean_bias_cols <- function(sim.out, base_names, true_override = NULL) {
  est_names <- paste0("est_", base_names)
  
  if (!all(est_names %in% names(sim.out))) {
    stop("Missing estimated columns in sim.out: ",
         paste(est_names[!est_names %in% names(sim.out)], collapse = ", "))
  }
  
  if (is.null(true_override)) {
    true_names <- paste0("true_", base_names)
    if (!all(true_names %in% names(sim.out))) {
      stop("Missing true columns in sim.out: ",
           paste(true_names[!true_names %in% names(sim.out)], collapse = ", "))
    }
    out <- vapply(seq_along(base_names), function(j) {
      mean(sim.out[[est_names[j]]] - sim.out[[true_names[j]]], na.rm = TRUE)
    }, numeric(1))
  } else {
    if (length(true_override) == 1L) true_override <- rep(true_override, length(base_names))
    if (length(true_override) != length(base_names)) {
      stop("'true_override' must have length 1 or the same length as 'base_names'.")
    }
    out <- vapply(seq_along(base_names), function(j) {
      mean(sim.out[[est_names[j]]] - true_override[j], na.rm = TRUE)
    }, numeric(1))
  }
  
  names(out) <- base_names
  out
}

sd_est_cols <- function(sim.out, base_names) {
  est_names <- paste0("est_", base_names)
  if (!all(est_names %in% names(sim.out))) {
    stop("Missing estimated columns in sim.out: ",
         paste(est_names[!est_names %in% names(sim.out)], collapse = ", "))
  }
  out <- vapply(est_names, function(nm) sd(sim.out[[nm]], na.rm = TRUE), numeric(1))
  names(out) <- base_names
  out
}

write_summary_line <- function(fields, file.out = NULL,
                               sep = " & ", end = "\\\\ \n") {
  line <- paste0(paste(fields, collapse = sep), end)
  if (is.null(file.out)) {
    cat(line)
  } else {
    cat(line, file = file.out, append = TRUE)
  }
}

summary.semi.gva.bias.f <- function(N, K, target_pi, target_event, sigma_eps, omega, weib_shape, beta_X, out_dir = "results", file.out){
  
  file.in <- file.path(out_dir,
                       paste0(
                         "GVA_INC_NN", N,
                         "_K", K,
                         "_pi_", paste(format(target_pi, trim = TRUE), collapse = "-"),
                         "_event_", format(target_event, trim = TRUE),
                         "_shape_", format(weib_shape, trim = TRUE),
                         "_omega_", format(omega, trim = TRUE),
                         "_beta_", paste(format(round(beta_X,2), trim = TRUE), collapse = "-"),
                         ".dat"
                       )
  )
  
  sim.out <- read.table(file.in, header=TRUE)
  
  bases <- gva_param_bases(K)
  true_sigma <- get_true_sigma_eps(sim.out, K, sigma_eps)
  
  bias_vals <- c(
    mean_bias_cols(sim.out, bases$beta),
    mean_bias_cols(sim.out, "omega"),
    mean_bias_cols(sim.out, "weib_shape"),
    mean_bias_cols(sim.out, "weib_lambda"),
    mean_bias_cols(sim.out, bases$gamma),
    mean_bias_cols(sim.out, bases$delta),
    mean_bias_cols(sim.out, bases$sigma_eps, true_override = true_sigma)
  )
  
  fields <- c(
    format(target_event, trim = TRUE),
    collapse_num(target_pi, digits = 2),
    collapse_num(beta_X, digits = 2),
    fmt3(bias_vals)
  )
  
  write_summary_line(fields, file.out = file.out, sep = " & ", end = "\\\\ \n")
  invisible(bias_vals)
}



summary.semi.gva.se.f <- function(N, K, target_pi, target_event, sigma_eps, omega, weib_shape, beta_X, out_dir = "results", file.out){
  file.in <- file.path(out_dir,
                       paste0(
                         "GVA_INC_NN", N,
                         "_K", K,
                         "_pi_", paste(format(target_pi, trim = TRUE), collapse = "-"),
                         "_event_", format(target_event, trim = TRUE),
                         "_shape_", format(weib_shape, trim = TRUE),
                         "_omega_", format(omega, trim = TRUE),
                         "_beta_", paste(format(round(beta_X,2), trim = TRUE), collapse = "-"),
                         ".dat"
                       )
  )
  
  sim.out <- read.table(file.in, header=TRUE)
  
  bases <- gva_param_bases(K)
  
  se_vals <- c(
    sd_est_cols(sim.out, bases$beta),
    sd_est_cols(sim.out, "omega"),
    sd_est_cols(sim.out, "weib_shape"),
    sd_est_cols(sim.out, "weib_lambda"),
    sd_est_cols(sim.out, bases$gamma),
    sd_est_cols(sim.out, bases$delta),
    sd_est_cols(sim.out, bases$sigma_eps)
  )
  
  fields <- c(
    format(target_event, trim = TRUE),
    collapse_num(target_pi, digits = 2),
    collapse_num(beta_X, digits = 2),
    fmt3(se_vals)
  )
  
  write_summary_line(fields, file.out = file.out, sep = " & ", end = "\\\\ \n")
  invisible(se_vals)
}



summary.semi.gva.f <- function(N, K, target_pi, target_event, sigma_eps, omega, weib_shape, beta_X, out_dir = "results", file.out){
  file.in <- file.path(out_dir,
                       paste0(
                         "GVA_INC_NN", N,
                         "_K", K,
                         "_pi_", paste(format(target_pi, trim = TRUE), collapse = "-"),
                         "_event_", format(target_event, trim = TRUE),
                         "_shape_", format(weib_shape, trim = TRUE),
                         "_omega_", format(omega, trim = TRUE),
                         "_beta_", paste(format(round(beta_X,2), trim = TRUE), collapse = "-"),
                         ".dat"
                       )
  )
  
  sim.out <- read.table(file.in, header=TRUE)
  
  bases <- gva_param_bases(K)
  true_sigma <- get_true_sigma_eps(sim.out, K, sigma_eps)
  
  bias_vals <- c(
    mean_bias_cols(sim.out, bases$beta),
    mean_bias_cols(sim.out, "omega"),
    mean_bias_cols(sim.out, "weib_shape"),
    mean_bias_cols(sim.out, "weib_lambda"),
    mean_bias_cols(sim.out, bases$gamma),
    mean_bias_cols(sim.out, bases$delta),
    mean_bias_cols(sim.out, bases$sigma_eps, true_override = true_sigma)
  )
  
  se_vals <- c(
    sd_est_cols(sim.out, bases$beta),
    sd_est_cols(sim.out, "omega"),
    sd_est_cols(sim.out, "weib_shape"),
    sd_est_cols(sim.out, "weib_lambda"),
    sd_est_cols(sim.out, bases$gamma),
    sd_est_cols(sim.out, bases$delta),
    sd_est_cols(sim.out, bases$sigma_eps)
  )
  
  cells <- paste0(fmt3(abs(bias_vals)), " (", fmt3(se_vals), ")")
  
  fields <- c(
    format(target_event, trim = TRUE),
    collapse_num(target_pi, digits = 2),
    collapse_num(beta_X, digits = 2),
    cells
  )
  
  write_summary_line(fields, file.out = file.out, sep = " & ", end = "\\\\ \n")
  invisible(list(abs_bias = abs(bias_vals), se = se_vals))
  
}






summary.semi.gva.betaonly.f <- function(N, K, target_pi, target_event, sigma_eps, omega, weib_shape, beta_X, out_dir = "results", file.out){
  file.in <- file.path(out_dir,
                       paste0(
                         "GVA_INC_NN", N,
                         "_K", K,
                         "_pi_", paste(format(target_pi, trim = TRUE), collapse = "-"),
                         "_event_", format(target_event, trim = TRUE),
                         "_shape_", format(weib_shape, trim = TRUE),
                         "_omega_", format(omega, trim = TRUE),
                         "_beta_", paste(format(round(beta_X,2), trim = TRUE), collapse = "-"),
                         ".dat"
                       )
  )
  
  sim.out <- read.table(file.in, header=TRUE)
  
  beta_names <- paste0("beta", seq_len(K))
  beta_bias  <- mean_bias_cols(sim.out, beta_names)
  beta_se    <- sd_est_cols(sim.out, beta_names)
  
  beta_cells <- paste0(fmt3(beta_bias), " (", fmt3(beta_se), ")")
  
  line_fields <- c(
    format(target_event, trim = TRUE),
    collapse_num(target_pi, digits = 2),
    collapse_num(beta_X, digits = 2),
    beta_cells
  )
  
  write_summary_line(line_fields, file.out = file.out, sep = " & ", end = "\\\\ \n")
  
  invisible(list(beta_bias = beta_bias, beta_se = beta_se))
}
