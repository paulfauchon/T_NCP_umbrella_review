# Trait → NCP second-order synthesis — modelling pipeline

Clean, reproducible Bayesian pipeline for the umbrella-review of plant
trait → Nature's Contributions to People (NCP) relationships. One Stan model,
selected by flags; numbered R scripts; every headline number traceable to a
script and a saved output.

## Quick start (self-contained)

Everything needed is in this folder — no external files.

1. Open `clean_pipeline.Rproj` in RStudio (anchors all paths).
2. `source("run_all.R")`.

`run_all.R` installs the R dependencies and CmdStan on first run, then executes
the pipeline. It starts in **QUICK mode** (`QUICK <- TRUE`: ~300 iterations,
primary model only — a few minutes) so you can confirm it runs end-to-end; set
`QUICK <- FALSE` for the real analysis, and `RUN_HEAVY <- TRUE` to add the
recovery simulation and full sensitivity battery. `data/sigma.csv` ships
prebuilt; set `BUILD_SIGMA <- TRUE` to regenerate it from the Excel.

## Directory layout

```
.
├── clean_pipeline.Rproj      # open this first (path anchor)
├── .here                     # path anchor for plain Rscript
├── run_all.R                 # one-click end-to-end runner
├── data/                     # ALL inputs included
│   ├── data_model.csv        # one row = one trait-NCP record within a review
│   ├── sigma.csv             # review dependence structure (prebuilt)
│   ├── study_scores.csv      # per-review methodological-quality scores
│   └── 2nd order SR Traits - Services.xlsx   # source for 00b_build_sigma.R
├── stan/
│   └── model_traitNCP.stan   # the single model (all variants via flags)
├── R/
│   ├── 00_config.R           # paths, seed, sampler, MODEL_SPECS
│   ├── _helpers.R            # censoring enumeration + variance decomposition
│   ├── 01_prepare_data.R     # build stan_data (+ integrity checks)
│   ├── 02_fit_models.R       # fit every model spec
│   ├── 03_diagnostics.R      # Rhat/ESS/divergences, PPC, LOO comparison
│   ├── 04_variance_decomposition.R   # THE central result (corrected)
│   ├── 05_grouping_sensitivity.R     # trait-grouping / exchangeability test
│   ├── 06_recovery_simulation.R      # simulation-based recovery + coverage
│   ├── 07_results_pairs.R    # per-pair posterior support vs vote counting
│   └── 08_identification_checks.R    # prior/posterior + Sigma-scale sensitivity
└── outputs/                  # created at runtime (fits + csv results)
```

## Run order

```r
source("R/00b_build_sigma.R")   # -> data/sigma.csv (from the primary-paper lists)
source("R/01_prepare_data.R")   # -> outputs/prepared_data.rds
source("R/02_fit_models.R")     # -> outputs/fits/fit_<name>.rds
source("R/03_diagnostics.R")    # -> outputs/model_comparison.csv
source("R/04_variance_decomposition.R")
source("R/05_grouping_sensitivity.R")
source("R/07_results_pairs.R")
source("R/08_sensitivity_battery.R")   # kappa/DM/link/priors/Sigma battery + corr(tau,kappa,gamma)
source("R/06_recovery_simulation.R")   # heavy (N_REP fits); run last
```

## The model (one file, two flags)

Latent 3-category cumulative-logit → Dirichlet-multinomial over
(Nneg, Nnull, Npos). For record *i* (pair *k*, review *j*):

```
eta_i   = mu_cat[c(k)] + mu_trait[t(k)] + mu_service[s(k)] + sigma_alpha·z_k + u_j
u       = X_rev·beta_rev + chol(kappa·Sigma + diag(tau_j²))·z_u
tau_j   = tau·exp(gamma·(1 − quality_j))          # quality → precision (qp)
```

No global intercept: the ordered thresholds carry the location (an intercept
would be non-identified against them — this is what made an earlier `mu_global`
variant ill-posed).

| flag | values |
|------|--------|
| `trait_structure` | 0 pair · 1 service · 2 trait+service · 3 organ+trait+service |
| `quality_mode`    | 0 none (nq) · 1 likelihood-weight (sensitivity only) · 2 precision (**qp, primary**) |
| `link`            | 0 logit (**primary**) · 1 probit (sensitivity) |
| `use_dm`          | 1 Dirichlet-multinomial (**primary**) · 0 multinomial / no overdispersion (sensitivity) |
| `estimate_kappa`  | 1 estimate the logit-scale per-study variance (**primary**, identified by the n_j spread) · 0 fix to `kappa_fixed` (used in the 08 robustness scan) |

`MODEL_SPECS` (in `00_config.R`) lists the six models fitted; **`qp` is primary.**
Prior scales (`ps_sigma`, `ps_tau`, `ps_gamma`, `ps_kappa_ln`) and fractional
study weights (`study_w`) are passed as **data**, so the whole robustness battery
runs by changing inputs, never the Stan code.

`mu_trait`, `mu_service`, `mu_cat` are **sum-to-zero** (location lives in the
ordered thresholds).

**Sigma and kappa.** `Sigma` is built by `00b_build_sigma.R` from the trait×service
primary-paper lists (tabs E1/E2/E3) as a UNITLESS dependence structure:
`Sigma[j,j]=1/n_j`, `Sigma[j,k]=shared(j,k)/(n_j n_k)`, PSD by construction. The
logit-scale per-study sampling variance is therefore the estimated parameter
`kappa` (`Var(u_j)=kappa/n_j + tau_j^2`). Because `n_j` spans ~12–213 (an 18×
spread), the heteroscedastic `kappa/n_j` term is distinguishable from the flat
`tau`, so **kappa is estimated in the primary model** (the earlier fix-to-1 was
only needed when `Sigma` carried an arbitrary scale). `R/08` reports
`corr(tau,kappa,gamma)` and the fixed-kappa robustness scan.

## Methods ↔ code correspondence

| Manuscript claim | Where |
|---|---|
| Latent ordinal DM likelihood, censored marginalisation | `stan/model_traitNCP.stan` |
| Quality as precision (not likelihood weight) | `quality_mode = 2`, `tau_j` |
| Review dependence via fixed Sigma + estimated scale kappa | `transformed parameters` (Cov) |
| Variance decomposition (headline) | `R/04` + `variance_decomposition()` |
| "5 models converge" robustness | `R/03` (LOO) |
| Trait-grouping does not drive the result | `R/05` |
| Result is not a shrinkage artefact | `R/06` (recovery + coverage) |
| Split is identified / robust to Sigma | `R/08` |

## Corrections applied after statistical review

Three independent reviewers (internal) audited the model and the decomposition.
Fixed here:

- **Variance decomposition put on one footing.** All components are now
  per-effect latent variances; the logistic link variance (π²/3) is in the full
  denominator; the DM overdispersion φ is reported separately. The previous
  `n_j/N` weighting made the sampling term cancel to `J/N` — removed.
- **Sampling floor made data-informed.** A scale `kappa` on `Sigma` is estimated
  (identified through the off-diagonals); `R/08` runs the fixed-`kappa`
  sensitivity so the headline is shown to be robust (or flagged as an assumption).
- **Recovery study fixed and strengthened.** Truth and estimator now share the
  per-effect basis; many replicates report bias *and* 90 % coverage; a dedicated
  scenario exercises the censored marginalisation path.
- **PPC uses the Dirichlet-multinomial** (was multinomial → underdispersed).
- **`gamma` prior tightened** to `half-N(0,0.5)` (weak identification / funnel).
- **Review order aligned once** across `j_id`, `quality_j`, `Sigma` (the original
  scripts built these from different orderings — a silent-mismatch risk).
- Sigma symmetrised + PD-checked; `rowSums(counts)==ntot` asserted; empty-window
  censoring guarded.

## Second review round (external: ChatGPT + Gemini) — actioned

- **Sigma rebuilt as a unitless structure** from the trait×service primary-paper
  lists (`00b_build_sigma.R`); `kappa` now carries the logit-scale per-study
  variance and is estimated (identified by the n_j spread 12–213). This resolves
  both the "Sigma on the wrong (count) scale" and the "kappa not identified"
  objections at once.
- **Full sensitivity battery** (`R/08`): kappa {0.5,1,2}, DM vs multinomial,
  logit vs probit, wide priors, Sigma diagonal +/-20% — with the genuine-
  heterogeneity share reported for each, plus posterior `corr(tau,kappa,gamma)`.
- **Sum-to-zero** on the level means; **prior scales and study weights are data**.
- **phi stays OUT of the latent denominator** (both external reviewers concur);
  reported separately, with a multinomial (no-overdispersion) sensitivity fit and
  the count of records that can inform phi (n>1).

## Known limitations to state in the manuscript (not fully fixable in code)

1. **[RESOLVED] Sigma scale.** Previously `diag(Sigma)=1/n_j` baked a count-scale
   variance into a log-odds effect. Now `Sigma` is a unitless *structure*
   (`00b_build_sigma.R`) and the logit-scale per-study variance is the estimated
   `kappa`. The construction script ships with the pipeline. Remaining caveat:
   the trait×service primary-paper lists (E1/E2/E3) were an "essayé" extraction —
   spot-check the reviews whose n_j moved most (e.g. 382, 54) before final submission.
2. **Double-counting of shared primary studies.** Overlapping studies enter
   several records' counts; `Sigma` correlates the review *effects* but the DM
   likelihood still treats overlapping *counts* as independent → precision is
   overstated and PSIS-LOO is optimistic. Proper fix needs primary-study IDs:
   deduplicate to unique studies with a crossed study×review effect. Report the
   number of unique primary studies vs the 1043 records.
3. **Review effect confounds synthesis with ecological context.** Reviews are
   ecosystem/climate-specific, so `u_j` mixes method and context. Plug review
   covariates into `X_rev` (`01_prepare_data.R`) and report residual
   heterogeneity after adjustment before calling it "synthesis noise".
4. **Direction-only.** The 3-category scheme discards effect magnitude; state
   as a limitation and, where magnitudes exist, run a subset sensitivity.
5. **Sigma is itself estimated, treated as known.** Its uncertainty is not
   propagated; the diagonal +/-20% battery in `R/08` bounds but does not remove
   this. Ideally bootstrap the overlap/sampling estimates and refit.
6. **Informative censoring (MNAR / publication bias).** Exact marginalisation
   assumes the censoring (null/mixed reporting) is ignorable given observed data.
   If negative effects are preferentially reported as null, mu and tau are biased.
   Justify the assumption; a selection-model tilt on the censored-config weights
   is the sensitivity to add.
```
