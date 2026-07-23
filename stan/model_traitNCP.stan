// ============================================================================
//  model_traitNCP.stan   (v3 — post external+internal review)
//  Hierarchical Bayesian latent-ordinal model for a second-order (umbrella)
//  synthesis of plant trait -> Nature's Contributions to People (NCP) evidence.
//
//  ONE model file; every model + every robustness variant is a set of data flags
//  so the entire sensitivity battery runs without touching the Stan code.
//
//  ---- structural flags ------------------------------------------------------
//   trait_structure : 0 pair | 1 service | 2 trait+service | 3 organ+trait+service
//   quality_mode    : 0 none (nq) | 1 likelihood-weight | 2 precision (qp, PRIMARY)
//   link            : 0 logit | 1 probit                      (sensitivity)
//   use_dm          : 1 Dirichlet-multinomial | 0 multinomial (sensitivity; phi off)
//   estimate_kappa  : 1 estimate Sigma scale | 0 fix to kappa_fixed (PRIMARY: fix=1)
//
//  ---- identifiability decisions (from statistical review) -------------------
//   * kappa is FIXED to 1 in the primary model. On this ultra-sparse design
//     (median n=1) kappa, tau and the fixed diag(Sigma) share one ridge; freeing
//     kappa worsens it. estimate_kappa=1 + a tight prior is a DIAGNOSTIC only.
//   * mu_trait / mu_service / mu_cat are SUM-TO-ZERO so the overall location sits
//     cleanly in the ordered thresholds (no global intercept; avoids a soft ridge
//     between the thresholds and the collective mean of the random effects).
//   * All prior scales are passed as DATA so prior-sensitivity is a data change.
//   * study_w are OPTIONAL fractional weights (1/multiplicity) to blunt the
//     double-counting of primary studies shared across reviews (interim fix; the
//     proper fix is a crossed study-level random effect once study IDs exist).
//
//  Review effects:  u = X_rev*beta_rev + chol(kappa*Sigma + diag(tau_j^2))*z_u.
//  Censored records handled by exact marginalisation. Requires Stan >= 2.32.
// ============================================================================

functions {
  // cumulative probability up to a threshold, logit (link=0) or probit (link=1)
  real cprob(real thr, real eta, int link) {
    return link == 0 ? inv_logit(thr - eta) : Phi(thr - eta);
  }
}

data {
  int<lower=1> N;
  int<lower=1> K_ts;
  int<lower=1> K_trait;
  int<lower=1> K_service;
  int<lower=1> K_cat;
  int<lower=1> J;
  int<lower=0> P;

  int<lower=0, upper=3> trait_structure;
  int<lower=0, upper=2> quality_mode;
  int<lower=0, upper=1> link;
  int<lower=0, upper=1> use_dm;
  int<lower=0, upper=1> estimate_kappa;
  real<lower=0>         kappa_fixed;

  // Prior scales (data-driven, for prior-sensitivity)
  real<lower=0> ps_sigma;     // half-normal SD for sigma_trait/service/cat/alpha
  real<lower=0> ps_tau;       // half-normal SD for tau
  real<lower=0> ps_gamma;     // half-normal SD for gamma
  real<lower=0> ps_kappa_ln;  // lognormal SD for kappa (when estimated)

  array[K_ts] int<lower=1, upper=K_trait>   trait_id;
  array[K_ts] int<lower=1, upper=K_service> service_id;
  array[K_ts] int<lower=1, upper=K_cat>     cat_id;

  array[N, 3] int          counts;
  array[N]    int<lower=1> n_trials;
  array[N]    int<lower=1, upper=K_ts> ts_id;
  array[N]    int<lower=1, upper=J>    j_id;
  vector<lower=0>[N]        study_w;         // fractional study weights (default 1)

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
  matrix[J, P]                X_rev;
}

transformed data {
  int use_service = (trait_structure >= 1);
  int use_trait   = (trait_structure >= 2);
  int use_cat     = (trait_structure >= 3);
  int use_gamma   = (quality_mode == 2);

  vector[N]      w_obs;
  vector[N_cens] w_cens;
  for (i in 1:N)      w_obs[i]  = study_w[i]      * ((quality_mode == 1) ? quality_j[j_id[i]]      : 1.0);
  for (i in 1:N_cens) w_cens[i] = study_w_cens[i] * ((quality_mode == 1) ? quality_j[j_id_cens[i]] : 1.0);
}

parameters {
  vector[K_ts]  z_alpha;
  real<lower=0> sigma_alpha;

  vector[K_trait]   z_trait;
  vector[K_service] z_service;
  vector[K_cat]     z_cat;
  real<lower=0>     sigma_trait;
  real<lower=0>     sigma_service;
  real<lower=0>     sigma_cat;

  vector[J]     z_u;
  real<lower=0> tau;
  real<lower=0> gamma;
  vector[P]     beta_rev;
  real<lower=0> kappa_raw;

  ordered[2] c;
  // Bounded to keep phi = exp(log_phi) finite: phi in [~0.0025, ~8100]. This
  // stops transient +Inf concentrations during warmup (harmless rejections, but
  // frequent because phi is weakly identified at n=1). The bounds sit far in the
  // tails of the normal(log 10, 1.5) prior, so the posterior is unaffected.
  real<lower=-6, upper=9> log_phi;
}

transformed parameters {
  real phi   = exp(log_phi);
  real kappa = estimate_kappa ? kappa_raw : kappa_fixed;

  // Sum-to-zero level means (centred so the thresholds carry the location)
  vector[K_trait]   mu_trait;
  vector[K_service] mu_service;
  vector[K_cat]     mu_cat;
  {
    vector[K_trait]   rt = use_trait   ? sigma_trait   * z_trait   : rep_vector(0, K_trait);
    vector[K_service] rs = use_service ? sigma_service * z_service : rep_vector(0, K_service);
    vector[K_cat]     rc = use_cat     ? sigma_cat     * z_cat     : rep_vector(0, K_cat);
    mu_trait   = rt - mean(rt);
    mu_service = rs - mean(rs);
    mu_cat     = rc - mean(rc);
  }

  vector<lower=0>[J] tau_j;
  for (j in 1:J)
    tau_j[j] = use_gamma ? tau * exp(gamma * (1 - quality_j[j])) : tau;

  vector[K_ts] alpha;
  for (k in 1:K_ts)
    alpha[k] = mu_cat[cat_id[k]] + mu_trait[trait_id[k]]
             + mu_service[service_id[k]] + sigma_alpha * z_alpha[k];

  vector[J] u;
  {
    matrix[J, J] Cov = kappa * Sigma;
    // + 1e-6 nugget guarantees positive-definiteness (removes transient warmup
    // Cholesky rejections when tau_j -> 0); negligible vs tau_j^2 and diag(Sigma).
    for (j in 1:J) Cov[j, j] += square(tau_j[j]) + 1e-6;
    u = (P > 0 ? X_rev * beta_rev : rep_vector(0, J)) + cholesky_decompose(Cov) * z_u;
  }
}

model {
  z_alpha   ~ std_normal();
  z_trait   ~ std_normal();
  z_service ~ std_normal();
  z_cat     ~ std_normal();
  z_u       ~ std_normal();

  sigma_alpha   ~ normal(0, ps_sigma);
  sigma_trait   ~ normal(0, ps_sigma);
  sigma_service ~ normal(0, ps_sigma);
  sigma_cat     ~ normal(0, ps_sigma);
  tau           ~ normal(0, ps_tau);
  gamma         ~ normal(0, ps_gamma);
  beta_rev      ~ normal(0, 1);
  kappa_raw     ~ lognormal(0, ps_kappa_ln);

  c       ~ normal(0, 2);
  log_phi ~ normal(log(10), 1.5);

  for (i in 1:N) {
    real eta = alpha[ts_id[i]] + u[j_id[i]];
    real a1  = cprob(c[1], eta, link);
    real a2  = cprob(c[2], eta, link);
    vector[3] p = [a1, a2 - a1, 1 - a2]';
    real ll = use_dm
      ? dirichlet_multinomial_lpmf(counts[i] | phi * p + 1e-6)
      : multinomial_lpmf(counts[i] | p);
    target += w_obs[i] * ll;
  }

  for (i in 1:N_cens) {
    real eta = alpha[ts_id_cens[i]] + u[j_id_cens[i]];
    real a1  = cprob(c[1], eta, link);
    real a2  = cprob(c[2], eta, link);
    vector[3] p = [a1, a2 - a1, 1 - a2]';
    int s = config_start[i];
    int m = config_count[i];
    vector[m] lp;
    for (g in 1:m)
      lp[g] = use_dm
        ? dirichlet_multinomial_lpmf(all_configs[s + g - 1] | phi * p + 1e-6)
        : multinomial_lpmf(all_configs[s + g - 1] | p);
    target += w_cens[i] * log_sum_exp(lp);
  }
}

generated quantities {
  array[N, 3] int y_rep;
  vector[N + N_cens] log_lik;          // unweighted, for honest PSIS-LOO

  for (i in 1:N) {
    real eta = alpha[ts_id[i]] + u[j_id[i]];
    real a1  = cprob(c[1], eta, link);
    real a2  = cprob(c[2], eta, link);
    vector[3] p = [a1, a2 - a1, 1 - a2]';
    if (use_dm) {
      vector[3] p_dm = dirichlet_rng(phi * p + 1e-6);
      y_rep[i]   = multinomial_rng(p_dm, n_trials[i]);
      log_lik[i] = dirichlet_multinomial_lpmf(counts[i] | phi * p + 1e-6);
    } else {
      y_rep[i]   = multinomial_rng(p, n_trials[i]);
      log_lik[i] = multinomial_lpmf(counts[i] | p);
    }
  }
  for (i in 1:N_cens) {
    real eta = alpha[ts_id_cens[i]] + u[j_id_cens[i]];
    real a1  = cprob(c[1], eta, link);
    real a2  = cprob(c[2], eta, link);
    vector[3] p = [a1, a2 - a1, 1 - a2]';
    int s = config_start[i];
    int m = config_count[i];
    vector[m] lp;
    for (g in 1:m)
      lp[g] = use_dm
        ? dirichlet_multinomial_lpmf(all_configs[s + g - 1] | phi * p + 1e-6)
        : multinomial_lpmf(all_configs[s + g - 1] | p);
    log_lik[N + i] = log_sum_exp(lp);
  }
}
