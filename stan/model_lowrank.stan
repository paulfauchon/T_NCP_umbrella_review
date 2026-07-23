// ============================================================================
//  model_lowrank.stan  —  Reduced-rank (latent-factor) trait x NCP model
//
//  A model MATCHED to the data structure: it lets a trait have OPPOSITE effects
//  on different NCPs (trade-offs) while still pooling — via K latent "syndrome"
//  axes — instead of the additive model's one-direction-per-trait assumption.
//
//    alpha[k] = sum_f L[trait(k), f] * F[service(k), f]  +  sigma_alpha * z[k]
//
//  where L (K_trait x K) are trait loadings, F (K_service x K) are NCP loadings,
//  and the small residual z[k] catches pairs off the low-rank structure.
//
//  Identification (standard Bayesian factor analysis): F is lower-triangular with
//  a positive diagonal — this fixes rotation and sign. L is free. Requires
//  K <= K_service.
//
//  Everything else (review MVN with kappa*Sigma + diag(tau_j^2), quality->
//  precision, Dirichlet-multinomial, exact censoring, ordered thresholds) is
//  identical to model_traitNCP.stan, so LOO is directly comparable.
//  Requires Stan >= 2.32.
// ============================================================================

functions {
  real cprob(real thr, real eta, int link) {
    return link == 0 ? inv_logit(thr - eta) : Phi(thr - eta);
  }
}

data {
  int<lower=1> N;
  int<lower=1> K_ts;
  int<lower=1> K_trait;
  int<lower=1> K_service;
  int<lower=1> J;
  int<lower=1> K;                      // latent rank (e.g. 2); K <= K_service

  int<lower=0, upper=2> quality_mode;  // 2 = precision (primary)
  int<lower=0, upper=1> link;          // 0 logit
  int<lower=0, upper=1> use_dm;        // 1 Dirichlet-multinomial
  int<lower=0, upper=1> estimate_kappa;
  real<lower=0>         kappa_fixed;

  real<lower=0> ps_sigma;
  real<lower=0> ps_tau;
  real<lower=0> ps_gamma;
  real<lower=0> ps_kappa_ln;

  array[K_ts] int<lower=1, upper=K_trait>   trait_id;
  array[K_ts] int<lower=1, upper=K_service> service_id;

  array[N, 3] int          counts;
  array[N]    int<lower=1> n_trials;
  array[N]    int<lower=1, upper=K_ts> ts_id;
  array[N]    int<lower=1, upper=J>    j_id;
  vector<lower=0>[N]        study_w;

  int<lower=0> N_cens;
  array[N_cens] int<lower=1, upper=K_ts> ts_id_cens;
  array[N_cens] int<lower=1, upper=J>    j_id_cens;
  vector<lower=0>[N_cens]   study_w_cens;
  int<lower=1>              total_configs;
  array[total_configs, 3]  int all_configs;
  array[N_cens] int<lower=1> config_start;
  array[N_cens] int<lower=1> config_count;

  vector<lower=0, upper=1>[J] quality_j;
  cov_matrix[J]               Sigma;
}

transformed data {
  int use_gamma = (quality_mode == 2);
  int n_Ffree = 0;                       // number of free below-diagonal F entries
  for (j in 1:K) n_Ffree += (K_service - j);

  vector[N]      w_obs;
  vector[N_cens] w_cens;
  for (i in 1:N)      w_obs[i]  = study_w[i]      * ((quality_mode == 1) ? quality_j[j_id[i]]      : 1.0);
  for (i in 1:N_cens) w_cens[i] = study_w_cens[i] * ((quality_mode == 1) ? quality_j[j_id_cens[i]] : 1.0);
}

parameters {
  matrix[K_trait, K]   L;              // trait loadings (free)
  vector<lower=0>[K]   F_diag;         // positive diagonal of F
  vector[n_Ffree]      F_free;         // free below-diagonal entries of F

  vector[K_ts]  z_alpha;               // small idiosyncratic pair residual
  real<lower=0> sigma_alpha;

  vector[J]     z_u;
  real<lower=0> tau;
  real<lower=0> gamma;
  real<lower=0> kappa_raw;

  ordered[2]              c;
  real<lower=-6, upper=9> log_phi;
}

transformed parameters {
  real phi   = exp(log_phi);
  real kappa = estimate_kappa ? kappa_raw : kappa_fixed;

  // Lower-triangular NCP loadings with positive diagonal (fixes rotation/sign)
  matrix[K_service, K] F = rep_matrix(0, K_service, K);
  {
    int pos = 1;
    for (j in 1:K) {
      F[j, j] = F_diag[j];
      for (i in (j + 1):K_service) { F[i, j] = F_free[pos]; pos += 1; }
    }
  }

  // Quality-dependent review SDs
  vector<lower=0>[J] tau_j;
  for (j in 1:J) tau_j[j] = use_gamma ? tau * exp(gamma * (1 - quality_j[j])) : tau;

  // Pair effects: low-rank trait x NCP interaction + small residual
  vector[K_ts] alpha;
  for (k in 1:K_ts)
    alpha[k] = dot_product(L[trait_id[k]], F[service_id[k]])
             + sigma_alpha * z_alpha[k];

  // Review effects
  vector[J] u;
  {
    matrix[J, J] Cov = kappa * Sigma;
    for (j in 1:J) Cov[j, j] += square(tau_j[j]) + 1e-6;   // PD nugget (negligible)
    u = cholesky_decompose(Cov) * z_u;
  }
}

model {
  to_vector(L) ~ std_normal();
  F_diag       ~ normal(0, 1);
  F_free       ~ std_normal();

  z_alpha      ~ std_normal();
  sigma_alpha  ~ normal(0, ps_sigma * 0.5);   // tight: factors do the work
  z_u          ~ std_normal();
  tau          ~ normal(0, ps_tau);
  gamma        ~ normal(0, ps_gamma);
  kappa_raw    ~ lognormal(0, ps_kappa_ln);

  c       ~ normal(0, 2);
  log_phi ~ normal(log(10), 1.5);

  for (i in 1:N) {
    real eta = alpha[ts_id[i]] + u[j_id[i]];
    real a1 = cprob(c[1], eta, link);
    real a2 = cprob(c[2], eta, link);
    vector[3] p = [a1, a2 - a1, 1 - a2]';
    real ll = use_dm ? dirichlet_multinomial_lpmf(counts[i] | phi * p + 1e-6)
                     : multinomial_lpmf(counts[i] | p);
    target += w_obs[i] * ll;
  }
  for (i in 1:N_cens) {
    real eta = alpha[ts_id_cens[i]] + u[j_id_cens[i]];
    real a1 = cprob(c[1], eta, link);
    real a2 = cprob(c[2], eta, link);
    vector[3] p = [a1, a2 - a1, 1 - a2]';
    int s = config_start[i]; int m = config_count[i];
    vector[m] lp;
    for (g in 1:m)
      lp[g] = use_dm ? dirichlet_multinomial_lpmf(all_configs[s + g - 1] | phi * p + 1e-6)
                     : multinomial_lpmf(all_configs[s + g - 1] | p);
    target += w_cens[i] * log_sum_exp(lp);
  }
}

generated quantities {
  vector[N + N_cens] log_lik;
  for (i in 1:N) {
    real eta = alpha[ts_id[i]] + u[j_id[i]];
    real a1 = cprob(c[1], eta, link); real a2 = cprob(c[2], eta, link);
    vector[3] p = [a1, a2 - a1, 1 - a2]';
    log_lik[i] = use_dm ? dirichlet_multinomial_lpmf(counts[i] | phi * p + 1e-6)
                        : multinomial_lpmf(counts[i] | p);
  }
  for (i in 1:N_cens) {
    real eta = alpha[ts_id_cens[i]] + u[j_id_cens[i]];
    real a1 = cprob(c[1], eta, link); real a2 = cprob(c[2], eta, link);
    vector[3] p = [a1, a2 - a1, 1 - a2]';
    int s = config_start[i]; int m = config_count[i];
    vector[m] lp;
    for (g in 1:m)
      lp[g] = use_dm ? dirichlet_multinomial_lpmf(all_configs[s + g - 1] | phi * p + 1e-6)
                     : multinomial_lpmf(all_configs[s + g - 1] | p);
    log_lik[N + i] = log_sum_exp(lp);
  }
}
