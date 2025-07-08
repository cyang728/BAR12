data {
  int<lower=0> N;              // Number of observations
  array[N] int<lower=0, upper=1> Y_T;  // Toxicity outcomes (0 or 1)
  array[N] int<lower=0, upper=1> Y_E;  // Efficacy outcomes (0 or 1)
  vector[N] d;                 // Standardized Doses
}

parameters {
  // Toxicity parameters
  real alpha0;                 // Intercept for toxicity
  real alpha1;                 // Slope for toxicity (log scale)

  // Efficacy hinge model parameters
  real beta0;                        // Intercept
  real<lower=0.0, upper=1.5> beta_peak; // Peak location for efficacy
  real betaL;                       // Magnitude of increasing slope for efficacy
  real betaR;                       // Magnitude of decreasing slope for efficacy
  real beta2;

  // Random effects (assuming shared random effect for simplicity)
  real<lower=0> sigma2_epsilon;
  vector[N] epsilon;
}

transformed parameters {
  vector[N] eta_T;             // Linear predictor for toxicity
  vector[N] eta_E;             // Linear predictor for efficacy
  vector[N] pi_T;              // Probability of toxicity
  vector[N] pi_E;              // Probability of efficacy

  // Linear predictor for toxicity (linear model)
  eta_T = alpha0 + exp(alpha1) * d + epsilon;

  // Linear predictor for Efficacy (change-point model)
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

  // Priors for toxicity (Y_T) - Renamed alpha0, alpha1 to alpha0_T, alpha1_T for clarity
  alpha0 ~ normal(-3.0, sqrt(5));
  alpha1 ~ normal( 1.1, sqrt(2));

  // Priors for efficacy (Y_E) - change-point model
  beta0 ~ normal(-1.0, sqrt(6));
  betaL ~ normal( 1.1, sqrt(2));
  betaR ~ normal(-5.0, sqrt(3));
  beta2 ~ normal(-0.1, sqrt(2));
  beta_peak ~ normal(0.8, sqrt(1));
}
