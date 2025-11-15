# BAR12 design for dofferent number of doses among manufacturing processes

library(rstan)

stage1_prior = stan_model(file = "/Users/chenghanyang/Downloads/Stage1_prior.stan")
stagek_prior = stan_model(file = "/Users/chenghanyang/Downloads/Stagek_prior.stan")

# Inverse Logistic
invlogit = function(x) {
  x = pmax(-100, pmin(100, x))
  1 / (1 + exp(-x))
}

# Outcome Model
outcome = function(idx_dose, n_size, prob_T, prob_E) {
  
  epsilon = rnorm(n_size, 0, 1) # Common random effect for T and E for a patient
  logit_pi_T = log(prob_T[idx_dose] / (1 - prob_T[idx_dose]))
  logit_pi_E = log(prob_E[idx_dose] / (1 - prob_E[idx_dose]))
  
  tmp_tox = logit_pi_T + epsilon
  tmp_eff = logit_pi_E + epsilon
  tmp_tox_prob = invlogit(tmp_tox)
  tmp_eff_prob = invlogit(tmp_eff)
  Y_T = as.numeric(runif(n_size) < tmp_tox_prob)
  Y_E = as.numeric(runif(n_size) < tmp_eff_prob)
  list(Y_T = Y_T, Y_E = Y_E)
}

reshape_matrix = function(mat) as.vector(t(mat))

# BAR12 updated
BAR12_design = function(dosages = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                        p_T_sim = list(c(0.01,0.02,0.03,0.03,0.06,0.08),
                                       c(0.01,0.02,0.03,0.05,0.08,0.28)),
                        p_E_sim = list(c(0.02,0.03,0.12,0.19,0.41,0.56),
                                       c(0.06,0.11,0.15,0.21,0.58,0.58)),
                        num_stages = 2,
                        ncohort = c(6, 8), cohortsize = c(3, 3), startdose = 1,
                        target_tox = 0.3, target_eff = 0.2,
                        cutoff_tox = 0.60, cutoff_eff = 0.85,
                        w00 = 40, w01 = 100, w10 = 0, w11 = 60,
                        max_allocate_dose = 100,
                        n_mc_epsilon = 1000,
                        kappa = 1.0,
                        L = 2,
                        seed = 1,
                        print.out = FALSE) {
  
  set.seed(seed)
  
  ## ---------- Define global dose grid & stage-specific available doses ----------
  
  dosages_input <- dosages
  is_list_dose  <- is.list(dosages_input)
  
  if (is_list_dose) {
    if (length(dosages_input) != num_stages) {
      stop("When 'dosages' is a list, its length must equal 'num_stages'.")
    }
    
    # Global union dose grid (sorted)
    dosages_global <- sort(unique(unlist(dosages_input)))
    num_doses <- length(dosages_global)
    
    # Global indices available in each stage
    available_set <- vector("list", num_stages)
    for (k in seq_len(num_stages)) {
      idx <- match(dosages_input[[k]], dosages_global)
      if (any(is.na(idx))) {
        stop(sprintf("Stage %d: some doses not found in pooled dose grid.", k))
      }
      available_set[[k]] <- sort(idx)
      
      if (length(p_T_sim[[k]]) != length(dosages_input[[k]]) ||
          length(p_E_sim[[k]]) != length(dosages_input[[k]])) {
        stop(sprintf("Stage %d: p_T_sim[[%d]] and p_E_sim[[%d]] must match length(dosages[[%d]]).",
                     k, k, k, k))
      }
    }
    
    # Map true probabilities onto the global grid (unused doses are NA and never sampled)
    p_T_global <- vector("list", num_stages)
    p_E_global <- vector("list", num_stages)
    for (k in seq_len(num_stages)) {
      pt <- rep(NA_real_, num_doses)
      pe <- rep(NA_real_, num_doses)
      pt[available_set[[k]]] <- p_T_sim[[k]]
      pe[available_set[[k]]] <- p_E_sim[[k]]
      p_T_global[[k]] <- pt
      p_E_global[[k]] <- pe
    }
    
  } else {
    dosages_global <- dosages_input
    num_doses <- length(dosages_global)
    
    if (length(p_T_sim) != num_stages || length(p_E_sim) != num_stages) {
      stop("p_T_sim and p_E_sim must be lists of length 'num_stages'.")
    }
    for (k in seq_len(num_stages)) {
      if (length(p_T_sim[[k]]) != num_doses ||
          length(p_E_sim[[k]]) != num_doses) {
        stop(sprintf("Stage %d: p_T_sim[[%d]] and p_E_sim[[%d]] must have length = length(dosages).",
                     k, k, k))
      }
    }
    
    available_set <- replicate(num_stages, seq_len(num_doses), simplify = FALSE)
    p_T_global <- p_T_sim
    p_E_global <- p_E_sim
  }
  
  # Use this global dose grid for all subsequent computations
  dosages <- dosages_global
  
  ## ---------- Initialize recording objects ----------
  
  npts = ncohort * cohortsize
  
  d_tol = y_T_tol = y_E_tol = NULL
  y_T = matrix(0, nrow = num_stages, ncol = num_doses)
  y_E = matrix(0, nrow = num_stages, ncol = num_doses)
  n   = matrix(0, nrow = num_stages, ncol = num_doses)
  earlystop = 0
  elimi = matrix(0, nrow = num_stages, ncol = num_doses)
  
  pi_T_hat = matrix(0, nrow = num_stages, ncol = num_doses)
  pi_E_hat = matrix(0, nrow = num_stages, ncol = num_doses)
  
  p_T_monitor = matrix(0, nrow = num_stages, ncol = num_doses)
  p_E_monitor = matrix(0, nrow = num_stages, ncol = num_doses)
  
  p_T1_E1 = p_T1_E0 = p_T0_E1 = p_T0_E0 = matrix(0, nrow = num_stages, ncol = num_doses)
  too_toxic_dose = NULL
  futile_dose = NULL
  
  weights_samples = matrix(1, ncol = num_doses, nrow = num_stages)
  t_rounds = 1
  
  ## Starting dose index: when 'dosages' is a list, convert Stage 1 local index to global index
  if (is_list_dose) {
    avail1 <- available_set[[1]]
    if (startdose < 1L || startdose > length(avail1)) {
      stop("When 'dosages' is a list, 'startdose' should be between 1 and length(dosages[[1]]).")
    }
    d <- avail1[startdose]
  } else {
    d <- startdose
  }
  
  ## ---------- Main loop over stages ----------
  
  for (k in seq_len(num_stages)) {
    
    avail_k <- available_set[[k]]
    
    # Ensure that the carried-over dose index from previous stage is available in this stage
    d <- intersect(d, avail_k)
    if (length(d) == 0L) {
      # Fallback: choose the middle dose among those available in this stage
      d <- avail_k[ceiling(length(avail_k) / 2)]
    }
    
    ncohort_k = ncohort[k]
    too_toxic_dose = NULL
    futile_dose = NULL
    d_tol_k = y_T_tol_k = y_E_tol_k = NULL
    
    i = 0
    while (i < ncohort_k) {
      
      ## ------- Randomization and allocation -------
      
      if (length(d) > 1L) {
        
        if ((sum(n[k,]) + cohortsize[k] * L) > npts[k]) {
          nremain = npts[k] - sum(n[k,])
          if (nremain > 0) {
            w <- Utility_mean[d]^kappa
            w <- w / sum(w)
            d_sample = sample(d, size = nremain,
                              prob = w,
                              replace = TRUE)
            d_sample = sort(d_sample)
          } else {
            break
          }
        } else {
          w <- Utility_mean[d]^kappa
          w <- w / sum(w)
          d_sample = sample(d, size = cohortsize[k] * L,
                            prob = w,
                            replace = TRUE)
          d_sample = sort(d_sample)
        }
        
        for (idx in seq_along(d_sample)) {
          dose_idx <- d_sample[idx]
          # Here dose_idx is guaranteed to be in avail_k, so p_T_global / p_E_global are not NA
          Y_out = outcome(idx_dose = dose_idx, n_size = 1,
                          prob_T = p_T_global[[k]],
                          prob_E = p_E_global[[k]])
          
          d_tol_k = c(d_tol_k, dose_idx)
          y_T_tol_k = c(y_T_tol_k, Y_out$Y_T)
          y_E_tol_k = c(y_E_tol_k, Y_out$Y_E)
          
          y_T[k, dose_idx] = y_T[k, dose_idx] + Y_out$Y_T
          y_E[k, dose_idx] = y_E[k, dose_idx] + Y_out$Y_E
          n[k, dose_idx]   = n[k, dose_idx]   + 1
        }
        i = i + length(d_sample) / cohortsize[k]
        
      } else {
        
        dose_idx <- d
        Y_out = outcome(idx_dose = dose_idx, n_size = cohortsize[k],
                        prob_T = p_T_global[[k]],
                        prob_E = p_E_global[[k]])
        
        d_tol_k   = c(d_tol_k, rep(dose_idx, cohortsize[k]))
        y_T_tol_k = c(y_T_tol_k, Y_out$Y_T)
        y_E_tol_k = c(y_E_tol_k, Y_out$Y_E)
        
        if ((sum(n[k,]) + cohortsize[k]) > npts[k]) {
          nremain = npts[k] - sum(n[k,])
          if (nremain > 0) {
            y_T[k, dose_idx] = y_T[k, dose_idx] + sum(Y_out$Y_T[1:nremain])
            y_E[k, dose_idx] = y_E[k, dose_idx] + sum(Y_out$Y_E[1:nremain])
            n[k, dose_idx]   = n[k, dose_idx]   + nremain
          }
          break
        } else {
          y_T[k, dose_idx] = y_T[k, dose_idx] + sum(Y_out$Y_T)
          y_E[k, dose_idx] = y_E[k, dose_idx] + sum(Y_out$Y_E)
          n[k, dose_idx]   = n[k, dose_idx]   + cohortsize[k]
        }
        
        i = i + 1
      }
      
      ## ------- Stan posterior estimation -------
      
      if (k == 1L) {
        
        stan_data = list(
          N   = sum(n[k,]),
          Y_T = y_T_tol_k,
          Y_E = y_E_tol_k,
          d   = dosages[d_tol_k]
        )
        
        stan_fit = rstan::sampling(stage1_prior,
                                   data   = stan_data,
                                   chains = 4,
                                   iter   = 10000,
                                   warmup = 5000,
                                   thin   = 2,
                                   seed   = 123,
                                   verbose = FALSE,
                                   refresh = 0)
        
      } else {
        
        stan_data = list(
          N   = sum(n[k,]),
          Y_T = y_T_tol_k,
          Y_E = y_E_tol_k,
          d   = dosages[d_tol_k],
          m_alpha0 = m_alpha0,
          precision_alpha0 = precision_alpha0,
          m_alpha1 = m_alpha1,
          precision_alpha1 = precision_alpha1,
          m_beta0  = m_beta0,
          precision_beta0  = precision_beta0,
          m_betaL  = m_betaL,
          precision_betaL  = precision_betaL,
          m_betaR  = m_betaR,
          precision_betaR  = precision_betaR,
          m_beta2  = m_beta2,
          precision_beta2  = precision_beta2,
          m_beta_peak = m_beta_peak,
          precision_beta_peak = precision_beta_peak
        )
        
        stan_fit = rstan::sampling(stagek_prior,
                                   data   = stan_data,
                                   chains = 4,
                                   iter   = 10000,
                                   warmup = 5000,
                                   thin   = 2,
                                   seed   = 123,
                                   verbose = FALSE,
                                   refresh = 0)
      }
      
      samples = rstan::extract(stan_fit)
      alpha0_samples        = samples$alpha0
      alpha1_samples        = samples$alpha1
      beta0_samples         = samples$beta0
      betaL_samples         = samples$betaL
      betaR_samples         = samples$betaR
      beta2_samples         = samples$beta2
      beta_peak_samples     = samples$beta_peak
      sigma2_epsilon_samples = samples$sigma2_epsilon
      
      utility_samples = matrix(NA_real_, nrow = length(alpha0_samples), ncol = num_doses)
      
      for (j in seq_len(num_doses)) {
        
        pi_T_marginal_samples_j = numeric(length(alpha0_samples))
        pi_E_marginal_samples_j = numeric(length(alpha0_samples))
        
        for (s in seq_along(alpha0_samples)) {
          
          epsilon = rnorm(n_mc_epsilon, 0, sd = sqrt(sigma2_epsilon_samples[s]))
          pi_T_j = invlogit(alpha0_samples[s] +
                              exp(alpha1_samples)[s] * dosages[j] +
                              epsilon)
          pi_E_j = invlogit(beta0_samples[s] +
                              exp(betaL_samples)[s] * dosages[j] +
                              betaR_samples[s] * pmax(0, dosages[j] - beta_peak_samples[s]) +
                              beta2_samples[s] * dosages[j]^2 +
                              epsilon)
          
          pi_T_marginal_samples_j[s] = mean(pi_T_j)
          pi_E_marginal_samples_j[s] = mean(pi_E_j)
          
          utility_samples[s, j] = mean(
            w11 * pi_T_j * pi_E_j +
              w00 * (1 - pi_T_j) * (1 - pi_E_j) +
              w01 * (1 - pi_T_j) * pi_E_j +
              w10 * pi_T_j * (1 - pi_E_j)
          )
        }
        
        pi_T_hat[k, j] = mean(pi_T_marginal_samples_j)
        pi_E_hat[k, j] = mean(pi_E_marginal_samples_j)
        
        p_T_monitor[k, j] = mean(pi_T_marginal_samples_j > target_tox)
        p_E_monitor[k, j] = mean(pi_E_marginal_samples_j < target_eff)
      }
      
      n_tmp = as.numeric(colSums(weights_samples * n))
      Utility_mean = colMeans(utility_samples)
      
      ## ------- Safety / futility rules -------
      
      elmi_monitor_T = which(p_T_monitor[k, ] > cutoff_tox)
      elmi_monitor_T = elmi_monitor_T[elmi_monitor_T %in% which(n[k, ] > 0)]
      elmi_monitor_E = which(p_E_monitor[k, ] > cutoff_eff)
      elmi_monitor_E = elmi_monitor_E[elmi_monitor_E %in% which(n[k, ] > 0)]
      
      if (length(elmi_monitor_T) != 0L) {
        min_tox_dose = min(elmi_monitor_T)
        too_toxic_dose = sort(unique(c(too_toxic_dose, min_tox_dose:num_doses)))
        if (n[k, min_tox_dose] > 0) elimi[k, too_toxic_dose] = 1
      }
      if (length(elmi_monitor_E) != 0L) {
        futile_dose = sort(unique(c(futile_dose, elmi_monitor_E)))
        elimi[k, futile_dose] = 3
      }
      elimi[k, n[k, ] >= max_allocate_dose] = 2
      
      if (all(elimi[k, avail_k] >= 1)) {
        earlystop = 1
        break
      }
      
      ## ------- Choose next dose index d -------
      
      if (k == 1L) {
        
        # No skipping doses: only look upward within the available doses in this stage
        pos_in_avail <- match(d, avail_k)
        if (!is.na(pos_in_avail) && pos_in_avail < length(avail_k)) {
          next_dose <- avail_k[pos_in_avail + 1L]
        } else {
          next_dose <- NA_integer_
        }
        
        if (!is.na(next_dose) &&
            n[1, next_dose] == 0 &&
            elimi[1, next_dose] == 0) {
          best_dose_k <- next_dose
        } else {
          admissible_dose_set <- intersect(which(elimi[k, ] == 0), avail_k)
          if (length(admissible_dose_set) == 0L) {
            break
          }
          best_dose_k <- admissible_dose_set[
            which.max(Utility_mean[admissible_dose_set])
          ]
        }
        
        if (d != best_dose_k) {
          if (elimi[k, best_dose_k] == 0) {
            d <- best_dose_k
          } else if (elimi[k, best_dose_k] >= 1L && elimi[k, d] >= 1L) {
            admissible_dose_set <- intersect(which(elimi[k, ] == 0), avail_k)
            if (length(admissible_dose_set) != 0L) {
              d <- admissible_dose_set[
                which.max(Utility_mean[admissible_dose_set])
              ]
            } else {
              break
            }
          }
        } else { # d == best_dose_k
          if (elimi[k, d] >= 1L) {
            admissible_dose_set <- intersect(which(elimi[k, ] == 0), avail_k)
            if (length(admissible_dose_set) != 0L) {
              d <- admissible_dose_set[
                which.max(Utility_mean[admissible_dose_set])
              ]
            } else {
              break
            }
          }
        }
        
      } else {
        
        admissible_dose_set <- intersect(which(elimi[k, ] == 0), avail_k)
        
        # In stage k, doses that have been used and are not labeled as 1 or 3 are considered admissible
        d_n <- intersect(which(n[k, ] > 0), avail_k)
        tmp_left  <- d_n - 1L
        tmp_right <- d_n + 1L
        tmp_left  <- tmp_left[tmp_left >= 1L & tmp_left <= num_doses]
        tmp_right <- tmp_right[tmp_right >= 1L & tmp_right <= num_doses]
        d_neigh   <- union(tmp_left, tmp_right)
        d_n <- sort(unique(c(d_n, d_neigh)))
        d_n <- intersect(d_n, avail_k)
        
        if (length(d_n) != 0L && length(admissible_dose_set) != 0L) {
          admissible_dose_set <- intersect(admissible_dose_set, d_n)
        } else {
          admissible_dose_set <- integer(0L)
        }
        
        if (length(admissible_dose_set) > 0L) {
          ordered_indices <- order(Utility_mean[admissible_dose_set], decreasing = TRUE)
          n_to_select <- min(L, length(admissible_dose_set))
          d <- sort(admissible_dose_set[ordered_indices[1:n_to_select]])
        } else {
          earlystop = 1
          break
        }
      }
      
      if (print.out) {
        cat("Stage:", k, "\n")
        cat("Npts: ", paste(n[k, ], collapse = ", "), "\n")
        cat("Utility: ", paste(round(Utility_mean, 3), collapse = ", "), "\n")
        cat("pi_T: ", paste(round(pi_T_hat[k, ], 3), collapse = ", "), "\n")
        cat("pi_E: ", paste(round(pi_E_hat[k, ], 3), collapse = ", "), "\n")
        cat("Monitor_Tox: ", paste(round(p_T_monitor[k, ], 3), collapse = ", "), "\n")
        cat("Monitor_Eff: ", paste(round(p_E_monitor[k, ], 3), collapse = ", "), "\n\n")
      }
      
    } # end while (cohorts)
    
    ## Aggregate into global records
    d_tol   = c(d_tol,   d_tol_k)
    y_T_tol = c(y_T_tol, y_T_tol_k)
    y_E_tol = c(y_E_tol, y_E_tol_k)
    
    ## Update hyper-parameters for the next stage
    m_alpha0 = mean(alpha0_samples); precision_alpha0 = 1 / var(alpha0_samples)
    m_alpha1 = mean(alpha1_samples); precision_alpha1 = 1 / var(alpha1_samples)
    m_beta0  = mean(beta0_samples);  precision_beta0  = 1 / var(beta0_samples)
    m_betaL  = mean(betaL_samples);  precision_betaL  = 1 / var(betaL_samples)
    m_betaR  = mean(betaR_samples);  precision_betaR  = 1 / var(betaR_samples)
    m_beta2  = mean(beta2_samples);  precision_beta2  = 1 / var(beta2_samples)
    m_beta_peak = mean(beta_peak_samples)
    precision_beta_peak = 1 / var(beta_peak_samples)
    
    if (earlystop == 1L) break
    
    ## ------- Choose the starting dose for the next stage (key: allow new doses) -------
    
    if (k < num_stages) {
      
      Utility_mean <- colMeans(utility_samples)
      
      avail_curr <- available_set[[k]]
      avail_next <- available_set[[k + 1L]]
      
      # In stage k, admissible doses that have been used and not labeled 1 or 3
      admissible_curr <- which((elimi[k, ] == 0 | elimi[k, ] == 2) & n[k, ] > 0)
      admissible_curr <- intersect(admissible_curr, avail_curr)
      
      if (length(admissible_curr) > 0L) {
        ordered_curr <- order(Utility_mean[admissible_curr], decreasing = TRUE)
        n_sel_curr <- min(L, length(admissible_curr))
        best_from_curr <- admissible_curr[ordered_curr[1:n_sel_curr]]
      } else {
        best_from_curr <- integer(0L)
      }
      
      # Doses that newly appear in the next stage
      new_in_next <- setdiff(avail_next, avail_curr)
      
      # Candidate set: top L doses from current stage plus newly introduced doses
      candidate_next <- sort(unique(c(best_from_curr, new_in_next)))
      if (length(candidate_next) == 0L) {
        candidate_next <- avail_next
      }
      
      ordered_next <- order(Utility_mean[candidate_next], decreasing = TRUE)
      n_to_select  <- min(L, length(candidate_next))
      d <- sort(candidate_next[ordered_next[1:n_to_select]])
    }
    
  } # end for stages
  
  ## ---------- Final recommended dose ----------
  
  last_stage <- max(which(rowSums(n) > 0))
  n_final    <- n[last_stage, ]
  Utility    <- colMeans(utility_samples)
  
  admissible_dose <- which(n_final != 0)
  
  elmi_monitor_T = too_toxic_dose
  elmi_monitor_E = futile_dose
  
  if (length(elmi_monitor_T) != 0L) {
    admissible_dose <- admissible_dose[
      !(admissible_dose %in% (min(elmi_monitor_T):num_doses))
    ]
  }
  if (length(elmi_monitor_E) != 0L) {
    admissible_dose <- admissible_dose[
      !(admissible_dose %in% elmi_monitor_E)
    ]
  }
  
  if (length(admissible_dose) != 0L) {
    best_idx_global <- admissible_dose[which.max(Utility[admissible_dose])]
  } else {
    best_idx_global <- 0L
  }
  
  best_dose_value <- if (best_idx_global > 0L) dosages[best_idx_global] else NA_real_
  
  ## ---------- Output data.frame ----------
  
  df = data.frame(
    cohort      = rep(seq_len(num_stages), each = num_doses),
    dose_index  = rep(seq_len(num_doses), times = num_stages),
    dose        = rep(dosages, times = num_stages),
    y_T         = reshape_matrix(y_T),
    y_E         = reshape_matrix(y_E),
    n           = reshape_matrix(n),
    pi_T_hat    = reshape_matrix(pi_T_hat),
    pi_E_hat    = reshape_matrix(pi_E_hat),
    p_T_monitor = reshape_matrix(p_T_monitor),
    p_E_monitor = reshape_matrix(p_E_monitor),
    Utility     = as.numeric(Utility),
    elimi       = reshape_matrix(elimi)
  )
  
  if (is_list_dose) {
    # Map global dose value back to the local index of the last stage
    last_doses <- dosages_input[[last_stage]]
    best_dose_stage <- if (!is.na(best_dose_value)) {
      match(best_dose_value, last_doses)
    } else {
      NA_integer_
    }
    
    output = list(
      data            = df,
      best_dose_index = best_dose_stage,     # Index within the last stage (corresponding to dosages_list[[last_stage]])
      best_dose_value = best_dose_value      # Actual standardized dose value
    )
  } else {
    output = list(
      data      = df,
      best_dose = best_idx_global            # Same as original: index in the 'dosages' vector
    )
  }
  
  return(output)
}

dosages_list <- list(
  c(0.0, 0.2, 0.4, 0.6, 0.8     , 1.0),     # Stage 1: 6 doses
  c(0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0)      # Stage 2: 7 doses
)

p_T_true <- list(
  c(0.01, 0.02, 0.03, 0.03, 0.05      , 0.38),
  c(0.02, 0.03, 0.04, 0.05, 0.06, 0.11, 0.35)
)

p_E_true <- list(
  c(0.02, 0.03, 0.19, 0.25, 0.41      , 0.41),
  c(0.08, 0.12, 0.15, 0.25, 0.35, 0.65, 0.35)
)

out <- BAR12_design(
  dosages = dosages_list,
  p_T_sim = p_T_true,
  p_E_sim = p_E_true,
  num_stages = 2,
  ncohort = c(6, 8),
  cohortsize = c(3, 3),
  seed = 1,
  print.out = FALSE
)

out$best_dose_value                # Recommended dose value in the last stage
