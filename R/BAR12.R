# Stan Models
path_stage1 = system.file("stan/Stage1_prior.stan", package = "BAR12")
path_stagek = system.file("stan/Stagek_prior.stan", package = "BAR12")

stage1_prior = stan_model(file = path_stage1)
stagek_prior = stan_model(file = path_stagek)

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

# BAR12
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
  num_doses = length(dosages)
  npts = ncohort * cohortsize

  # Record data
  d_tol = y_T_tol = y_E_tol = NULL
  y_T = matrix(0, nrow = num_stages, ncol = num_doses)
  y_E = matrix(0, nrow = num_stages, ncol = num_doses)
  n   = matrix(0, nrow = num_stages, ncol = num_doses)
  earlystop = 0
  d = startdose
  elimi = matrix(0, nrow = num_stages, ncol = num_doses)

  # Estimation
  pi_T_hat = matrix(0, nrow = num_stages, ncol = num_doses)
  pi_E_hat = matrix(0, nrow = num_stages, ncol = num_doses)

  # Monitor
  p_T_monitor = matrix(0, nrow = num_stages, ncol = num_doses)
  p_E_monitor = matrix(0, nrow = num_stages, ncol = num_doses)

  p_T1_E1 = p_T1_E0 = p_T0_E1 = p_T0_E0 = matrix(0, nrow = num_stages, ncol = num_doses)
  too_toxic_dose = NULL
  futile_dose = NULL

  weights_samples = matrix(1, ncol = num_doses, nrow = num_stages) #rep(1, num_stages)
  t_rounds = 1

  for(k in 1:num_stages){

    ncohort_k = ncohort[k]
    too_toxic_dose = NULL
    futile_dose = NULL
    d_tol_k = y_T_tol_k = y_E_tol_k = NULL

    i = 0
    while (i < ncohort_k) {

      if(length(d) > 1){

        if ((sum(n[k,]) + cohortsize[k]*L) > npts[k]) {
          nremain = npts[k] - sum(n[k,])
          if(nremain > 0){
            d_sample = sample(d, size = nremain,
                              prob = Utility_mean[d]^kappa/sum(Utility_mean[d]^kappa),
                              replace = TRUE)
            d_sample = sort(d_sample)
          }else{
            break
          }
        }else{
          d_sample = sample(d, size = cohortsize[k]*L,
                            prob = Utility_mean[d]^kappa/sum(Utility_mean[d]^kappa),
                            replace = TRUE)
          d_sample = sort(d_sample)
        }

        for(idx in 1:length(d_sample)){
          Y_out = outcome(idx_dose = d_sample[idx], n_size = 1,
                          prob_T = p_T_sim[[k]],
                          prob_E = p_E_sim[[k]])

          d_tol_k = c(d_tol_k, d_sample[idx])
          y_T_tol_k = c(y_T_tol_k, Y_out$Y_T)
          y_E_tol_k = c(y_E_tol_k, Y_out$Y_E)

          y_T[k,d_sample[idx]] = y_T[k,d_sample[idx]] + sum(Y_out$Y_T)
          y_E[k,d_sample[idx]] = y_E[k,d_sample[idx]] + sum(Y_out$Y_E)
          n[k,d_sample[idx]] = n[k,d_sample[idx]] + 1
        }
        i = i + length(d_sample)/(cohortsize[k])

      }else{
        Y_out = outcome(idx_dose = d, n_size = cohortsize[k],
                        prob_T = p_T_sim[[k]],
                        prob_E = p_E_sim[[k]])

        d_tol_k = c(d_tol_k, rep(d, cohortsize[k]))
        y_T_tol_k = c(y_T_tol_k, Y_out$Y_T)
        y_E_tol_k = c(y_E_tol_k, Y_out$Y_E)

        if ((sum(n[k,]) + cohortsize[k]) > npts[k]) {
          nremain = npts[k] - sum(n[k,])
          y_T[k,d] = y_T[k,d] + sum(Y_out$Y_T[1:nremain])
          y_E[k,d] = y_E[k,d] + sum(Y_out$Y_E[1:nremain])
          n[k,d] = n[k,d] + nremain
          break
        } else {
          y_T[k,d] = y_T[k,d] + sum(Y_out$Y_T)
          y_E[k,d] = y_E[k,d] + sum(Y_out$Y_E)
          n[k,d] = n[k,d] + cohortsize[k]
        }

        i = i + 1
      }

      if(k == 1){
        # Estimate the posterior probability of toxicity and efficacy
        stan_data = list(
          N = sum(n[k,]),
          Y_T = y_T_tol_k,
          Y_E = y_E_tol_k,
          d = dosages[d_tol_k]
        )

        stan_fit = rstan::sampling(stage1_prior,
                                   data = stan_data,
                                   chains = 4,
                                   iter = 10000,
                                   warmup = 5000,
                                   thin = 2,
                                   seed = 123,
                                   verbose = FALSE,  # Suppress detailed messages
                                   refresh = 0)      # Disable progress updates)

      }else{

        stan_data = list(
          N = sum(n[k,]),
          Y_T = y_T_tol_k,
          Y_E = y_E_tol_k,
          d = dosages[d_tol_k],
          m_alpha0 = m_alpha0,
          precision_alpha0 = precision_alpha0,
          m_alpha1 = m_alpha1,
          precision_alpha1 = precision_alpha1,
          m_beta0 = m_beta0,
          precision_beta0 = precision_beta0,
          m_betaL = m_betaL,
          precision_betaL = precision_betaL,
          m_betaR = m_betaR,
          precision_betaR = precision_betaR,
          m_beta2 = m_beta2,
          precision_beta2 = precision_beta2,
          m_beta_peak = m_beta_peak,
          precision_beta_peak = precision_beta_peak
        )

        stan_fit = rstan::sampling(stagek_prior,
                                   data = stan_data,
                                   chains = 4,
                                   iter = 10000,
                                   warmup = 5000,
                                   thin = 2,
                                   seed = 123,
                                   verbose = FALSE,  # Suppress detailed messages
                                   refresh = 0)      # Disable progress updates
      }

      samples = rstan::extract(stan_fit)
      alpha0_samples = samples$alpha0
      alpha1_samples = samples$alpha1
      beta0_samples = samples$beta0
      betaL_samples = samples$betaL
      betaR_samples = samples$betaR
      beta2_samples = samples$beta2
      beta_peak_samples = samples$beta_peak
      sigma2_epsilon_samples = samples$sigma2_epsilon

      utility_samples = matrix(NA, nrow = length(alpha0_samples), ncol = num_doses)

      for (j in 1:num_doses) {

        pi_T_marginal_samples_j = numeric(length(alpha0_samples))
        pi_E_marginal_samples_j = numeric(length(alpha0_samples))

        for(s in 1:length(alpha0_samples)){

          epsilon = rnorm(n_mc_epsilon, 0, sd = sqrt(sigma2_epsilon_samples[s]))
          pi_T_j = invlogit(alpha0_samples[s] + exp(alpha1_samples)[s] * dosages[j] + epsilon)
          pi_E_j = invlogit(beta0_samples[s] +
                              exp(betaL_samples)[s] * dosages[j] +
                              betaR_samples[s] * pmax(0, dosages[j] - beta_peak_samples[s]) + beta2_samples[s] * dosages[j]^2 +
                              epsilon)

          pi_T_marginal_samples_j[s] = mean(pi_T_j)
          pi_E_marginal_samples_j[s] = mean(pi_E_j)

          #utility_samples[s,j] = mean(pi_T_j <= target_tox) * mean(pi_E_j > target_eff)
          utility_samples[s,j] = mean(w11 * pi_T_j * pi_E_j +
                                        w00 * (1 - pi_T_j) * (1 - pi_E_j) +
                                        w01 * (1 - pi_T_j) * pi_E_j +
                                        w10 * pi_T_j * (1 - pi_E_j))
        }

        pi_T_hat[k,j] = mean(pi_T_marginal_samples_j)
        pi_E_hat[k,j] = mean(pi_E_marginal_samples_j)

        p_T_monitor[k,j] = mean(pi_T_marginal_samples_j > target_tox)
        p_E_monitor[k,j] = mean(pi_E_marginal_samples_j < target_eff)
      }

      n_tmp = as.numeric(colSums(weights_samples * n))
      Utility_mean = colMeans(utility_samples)

      elmi_monitor_T = which(p_T_monitor[k,] > cutoff_tox)
      elmi_monitor_T = elmi_monitor_T[elmi_monitor_T %in% which(n[k,] > 0)]
      elmi_monitor_E = which(p_E_monitor[k,] > cutoff_eff)
      elmi_monitor_E = elmi_monitor_E[elmi_monitor_E %in% which(n[k,] > 0)]

      if (length(elmi_monitor_T) != 0) {
        min_tox_dose = min(elmi_monitor_T)
        too_toxic_dose = c(too_toxic_dose, min_tox_dose:num_doses)
        too_toxic_dose = sort(unique(too_toxic_dose))
        if (n[k,min_tox_dose] > 0) elimi[k,too_toxic_dose] = 1
      }
      if (length(elmi_monitor_E) != 0) {
        futile_dose = c(futile_dose, elmi_monitor_E)
        futile_dose = sort(unique(futile_dose))
        elimi[k,futile_dose] = 3
      }
      elimi[k, n[k,] >= max_allocate_dose] = 2
      if(all(elimi[k,] >= 1)) break

      if(k == 1){
        if (d < num_doses) {
          if (n[1,d + 1] == 0 & elimi[1,d + 1] == 0) {
            best_dose = d + 1
          } else {
            admissible_dose_set = which(elimi[k,] == 0)
            best_dose = admissible_dose_set[which.max(Utility_mean[admissible_dose_set])]
          }
        } else if (d == num_doses) {
          admissible_dose_set = which(elimi[k,] == 0)
          best_dose = admissible_dose_set[which.max(Utility_mean[admissible_dose_set])]
        }

        if (d != best_dose) {
          if (elimi[k, best_dose] == 0) {
            d = best_dose
          } else if (elimi[k, best_dose] >= 1 & elimi[k, d] >= 1) {
            admissible_dose_set = which(elimi[k,] == 0)
            if (length(admissible_dose_set) != 0) {
              d = admissible_dose_set[which.max(Utility_mean[admissible_dose_set])]
            } else {
              break
            }
          }
        }else if (d == best_dose) {
          if (elimi[k,d] == 0) {
            d = d
          } else if (elimi[k,d] >= 1) {
            admissible_dose_set = which(elimi[k,] == 0)
            if (length(admissible_dose_set) != 0) {
              d = admissible_dose_set[which.max(Utility_mean[admissible_dose_set])]
            } else {
              break
            }
          }
        }
      }else{

        admissible_dose_set = which((elimi[k,] == 0))

        ###############################################################
        d_n = which(n[k,] > 0)
        tmp_left = d_n - 1
        tmp_right= d_n + 1
        tmp_left = tmp_left[tmp_left >= 1 & tmp_left<=num_doses]
        tmp_right= tmp_right[tmp_right >= 1 & tmp_right<=num_doses]
        d_n = union(d_n, union(tmp_left, tmp_right))
        d_n = sort(d_n)
        if(length(d_n)!=0 & length(admissible_dose_set)!=0){
          admissible_dose_set = intersect(d_n, admissible_dose_set)
        }else{
          admissible_dose_set = integer(0)
        }
        ###############################################################

        if (length(admissible_dose_set) > 0) {  # Ensure there are admissible doses
          # Order admissible doses by Utility in descending order
          ordered_indices = order(Utility_mean[admissible_dose_set], decreasing = TRUE)
          # Select top L or all if fewer than L
          n_to_select = min(L, length(admissible_dose_set))
          d = admissible_dose_set[ordered_indices[1:n_to_select]]
          d = sort(d)
        } else {
          break
        }
      }

      if(print.out){
        print(paste0("Npts: ", paste(n[k,], collapse = ", ")))
        print(paste0("Utility: ", paste(round(as.numeric(Utility_mean), 3), collapse = ", ")))
        print(paste0("pi_T: ", paste(round(as.numeric(pi_T_hat[k,]), 3), collapse = ", ")))
        print(paste0("pi_E: ", paste(round(as.numeric(pi_E_hat[k,]), 3), collapse = ", ")))
        print(paste0("Monitor_Tox: ", paste(round(as.numeric(p_T_monitor[k,]), 3), collapse = ", ")))
        print(paste0("Monitor_Eff: ", paste(round(as.numeric(p_E_monitor[k,]), 3), collapse = ", ")))
      }

    }

    admissible_dose_set = which((elimi[k,] == 0 | elimi[k,] == 2) & n[k,] > 0)
    Utility_mean = colMeans(utility_samples)

    # Select top L doses (or all if fewer than L)
    if (length(admissible_dose_set) > 0) {  # Ensure there are admissible doses
      # Order admissible doses by Utility in descending order
      ordered_indices = order(Utility_mean[admissible_dose_set], decreasing = TRUE)
      # Select top L or all if fewer than L
      n_to_select = min(L, length(admissible_dose_set))
      d = admissible_dose_set[ordered_indices[1:n_to_select]]
      d = sort(d)
    } else {
      break
    }

    #
    d_tol = c(d_tol, d_tol_k)
    y_T_tol = c(y_T_tol, y_T_tol_k)
    y_E_tol = c(y_E_tol, y_E_tol_k)

    m_alpha0 = mean(alpha0_samples); precision_alpha0 = 1/var(alpha0_samples)
    m_alpha1 = mean(alpha1_samples); precision_alpha1 = 1/var(alpha1_samples)
    m_beta0  = mean(beta0_samples) ; precision_beta0  = 1/var(beta0_samples)
    m_betaL  = mean(betaL_samples) ; precision_betaL  = 1/var(betaL_samples)
    m_betaR  = mean(betaR_samples) ; precision_betaR  = 1/var(betaR_samples)
    m_beta2  = mean(beta2_samples) ; precision_beta2  = 1/var(beta2_samples)
    m_beta_peak = mean(beta_peak_samples);
    precision_beta_peak = 1/var(beta_peak_samples);

  }

  admissible_dose = which(n[num_stages,] != 0)
  Utility = colMeans(utility_samples)

  elmi_monitor_T = too_toxic_dose
  elmi_monitor_E = futile_dose

  if (length(elmi_monitor_T) != 0) {
    admissible_dose = admissible_dose[!(admissible_dose %in% (min(elmi_monitor_T):num_doses))]
  }
  if (length(elmi_monitor_E) != 0) {
    admissible_dose = admissible_dose[!(admissible_dose %in% elmi_monitor_E)]
  }
  if (length(admissible_dose) != 0) {
    best_dose = admissible_dose[which.max(Utility[admissible_dose])]
  } else {
    best_dose = 0
  }

  df = data.frame(
    cohort = rep(1:num_stages, each = num_doses),
    dose_index = rep(1:num_doses, times = num_stages),
    dose = rep(dosages, times = num_stages),
    y_T = reshape_matrix(y_T),
    y_E = reshape_matrix(y_E),
    n = reshape_matrix(n),
    pi_T_hat = reshape_matrix(pi_T_hat),
    pi_E_hat = reshape_matrix(pi_E_hat),
    p_T_monitor = reshape_matrix(p_T_monitor),
    p_E_monitor = reshape_matrix(p_E_monitor),
    Utility = as.numeric(Utility),
    elimi = reshape_matrix(elimi)
  )

  output = list(data = df,
                best_dose = best_dose)

  return(output)
}
