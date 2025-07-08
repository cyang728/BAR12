// Tox_Eff_v2_k.stan
data {
  int<lower=0> N;              // Number of observations
  array[N] int<lower=0, upper=1> Y_T;  // Toxicity outcomes (0 or 1)
  array[N] int<lower=0, upper=1> Y_E;  // Efficacy outcomes (0 or 1)
  vector[N] d;                 // Dosage covariate
  real m_alpha0;               // Spike mean for alpha0
  real<lower=0> precision_alpha0; // Spike precision for alpha0
  real m_alpha1;               // Spike mean for alpha1
  real<lower=0> precision_alpha1; // Spike precision for alpha1
  real m_beta0;                // Spike mean for beta0
  real<lower=0> precision_beta0;  // Spike precision for beta0
  real m_betaL;                // Spike mean for betaL
  real<lower=0> precision_betaL;  // Spike precision for betaL
  real m_betaR;                // Spike mean for betaR
  real<lower=0> precision_betaR;  // Spike precision for betaR
  real m_beta2;                // Spike mean for betaR
  real<lower=0> precision_beta2;  // Spike precision for betaR
  real m_beta_peak;
  real<lower=0> precision_beta_peak;  // Spike precision for beta_peak
}

parameters {
  real<lower=0> sigma2_epsilon;
  vector[N] epsilon;           // Random effects
  real<lower=0,upper=1> z_alpha0; // Mixture weight for alpha0
  real<lower=0,upper=1> z_beta0;  // Mixture weight for beta0
  real alpha0_spike;           // Spike component for alpha0
  real alpha0_slab;            // Slab component for alpha0
  real alpha1_spike;           // Spike component for alpha1
  real alpha1_slab;            // Slab component for alpha1
  real beta0_spike;            // Spike component for gamma0
  real beta0_slab;             // Slab component for gamma0
  real betaL_spike;            // Spike component for gamma1
  real betaL_slab;             // Slab component for gamma1
  real betaR_spike;            // Spike component for gamma2
  real betaR_slab;             // Slab component for gamma2
  real beta2_spike;            // Spike component for gamma2
  real beta2_slab;             // Slab component for gamma2
  real<lower=0.0, upper=1.5> beta_peak_spike;
  real<lower=0.0, upper=1.5> beta_peak_slab;
  //real<lower=0, upper=1> z_peak;
}

transformed parameters {
  real alpha0;                 // Intercept for toxicity
  real alpha1;                 // Slope for toxicity (log scale)
  real beta0;                  // Intercept for efficacy
  real betaL;                  // Linear term for efficacy
  real betaR;                  // Quadratic term for efficacy
  real beta2;
  real<lower=0.0, upper=1.5> beta_peak;
  vector[N] eta_T;             // Linear predictor for toxicity
  vector[N] eta_E;             // Linear predictor for efficacy
  vector[N] pi_T;              // Probability of toxicity
  vector[N] pi_E;              // Probability of efficacy

  // Mixture for coefficients
  alpha0 = z_alpha0 * alpha0_slab + (1 - z_alpha0) * alpha0_spike;
  alpha1 = z_alpha0 * alpha1_slab + (1 - z_alpha0) * alpha1_spike;
  beta0  = z_beta0 * beta0_slab + (1 - z_beta0) * beta0_spike;
  betaL  = z_beta0 * betaL_slab + (1 - z_beta0) * betaL_spike;
  betaR  = z_beta0 * betaR_slab + (1 - z_beta0) * betaR_spike;
  beta2  = z_beta0 * beta2_slab + (1 - z_beta0) * beta2_spike;
  beta_peak = z_beta0 * beta_peak_slab + (1 - z_beta0) * beta_peak_spike;

  // Linear predictors
  eta_T = alpha0 + exp(alpha1) * d + epsilon;
  eta_E = beta0 + exp(betaL) * d + betaR * fmax(0, d - beta_peak) + beta2 * d .* d + epsilon;

  // Inverse logit transformation
  pi_T = inv_logit(eta_T);
  pi_E = inv_logit(eta_E);
}

model {
  // Likelihood
  Y_T ~ bernoulli_logit(eta_T); // Using bernoulli_logit for numerical stability
  Y_E ~ bernoulli_logit(eta_E); // Using bernoulli_logit for numerical stability

  // Random effects
  sigma2_epsilon ~ cauchy(0, 1);
  epsilon ~ normal(0, sqrt(sigma2_epsilon));

  // Continuous mixture weights
  z_alpha0 ~ beta(1, 1);
  z_beta0 ~ beta(1, 1);

  // Spike and slab priors
  alpha0_spike ~ normal(m_alpha0, sqrt(1 / precision_alpha0));
  alpha0_slab  ~ normal(-3.0, sqrt(5));

  alpha1_spike ~ normal(m_alpha1, sqrt(1 / precision_alpha1));
  alpha1_slab  ~ normal( 1.1, sqrt(2));

  beta0_spike  ~ normal(m_beta0, sqrt(1 / precision_beta0));
  beta0_slab   ~ normal(-1.0, sqrt(6));

  betaL_spike  ~ normal(m_betaL, sqrt(1 / precision_betaL));
  betaL_slab   ~ normal( 1.1, sqrt(2));

  betaR_spike  ~ normal(m_betaR, sqrt(1 / precision_betaR));
  betaR_slab   ~ normal(-5.0, sqrt(3));

  beta2_spike  ~ normal(m_beta2, sqrt(1 / precision_beta2));
  beta2_slab   ~ normal(-0.1, sqrt(2));

  beta_peak_spike ~ normal(m_beta_peak, sqrt(0.1 / precision_beta_peak));
  beta_peak_slab  ~ normal(0.8, 1);
}
