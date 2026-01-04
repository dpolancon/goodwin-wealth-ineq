# ============================================================
# 02_stage2_finance.R  (AUDITED: SOFT vs FATAL GATING + PHI1 LOCK)
# Stage 2: Finance/discipline objects + phi1 calibration (fixed)
# Reads Stage 1 survivors from CSV
# Outputs -> outputs/wealth_goodwin/stage2_finance/stage2_finance_scan.csv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(tibble)
  library(stringr)
})

# Stage 00 globals (paths, params, cfg, functions)
if (!exists("dirs") || !exists("par_base")) source("00_globals_utils.R", local = FALSE)

dir_stage1 <- dirs$stage1_backbone
dir_stage2 <- dirs$stage2_finance
dir.create(dir_stage2, showWarnings = FALSE, recursive = TRUE)

cat("\n====================\nStage 2: Finance/Discipline + phi1 calibration (soft/fatal gating)\n====================\n")

stage1_file <- file.path(dir_stage1, "stage1_backbone.csv")
stopifnot(file.exists(stage1_file))

stage1_res <- readr::read_csv(stage1_file, show_col_types = FALSE)
stage1_ok  <- stage1_res %>% filter(backbone_ok)

# ----------------------------
# Output schema
# ----------------------------
stage2_schema <- tibble(
  try_id = integer(),
  cand_id = integer(),
  sigma = double(), g_n = double(), i = double(), delta = double(),
  kappa_min = double(), kappa_max = double(),
  kappa_star = double(), r_star = double(), d_star = double(), omega_star = double(),
  rF = double(), psi = double(), phi2 = double(),
  lambda_star = double(), f_star = double(), Z_star = double(),
  
  phi1_fixed = double(),
  e_implied  = double(),
  phi1_calib_ok = logical(),
  phi1_calib_reason = character(),
  
  econ_ok = logical(),
  reason = character()
)

# ----------------------------
# Phi1 calibration settings
# ----------------------------
phi1_prior_min   <- par_base$phi1_min %||% 0.10
phi1_prior_max   <- par_base$phi1_max %||% 5.00
phi1_draws       <- if (exists("phi1_draws_stage2")) phi1_draws_stage2 else 400L

e_center <- cfg$targets$e_target %||% 0.94
e_sd     <- if (exists("e_target_sd")) e_target_sd else 0.05
e_min_ok <- if (exists("e_min_ok")) e_min_ok else 0.50
e_max_ok <- if (exists("e_max_ok")) e_max_ok else 1.20

w_min_ok <- if (exists("phi1_weight_min_ok")) phi1_weight_min_ok else 1e-6

safe_num <- function(x) suppressWarnings(as.numeric(x))

calibrate_phi1_once <- function(A, phi1_min, phi1_max, e_target, e_sd, e_min_ok, e_max_ok, draws=400L) {
  if (!is.finite(A)) return(list(phi1_fixed=NA_real_, e_implied=NA_real_, ok=FALSE, reason="bad_A"))
  if (!is.finite(phi1_min) || !is.finite(phi1_max) || phi1_min <= 0 || phi1_max <= phi1_min) {
    return(list(phi1_fixed=NA_real_, e_implied=NA_real_, ok=FALSE, reason="bad_phi1_bounds"))
  }
  
  phi1_draw <- runif(draws, min=phi1_min, max=phi1_max)
  e_hat <- A / phi1_draw
  
  ok_e <- is.finite(e_hat) & (e_hat >= e_min_ok) & (e_hat <= e_max_ok)
  
  w <- rep(0, length(phi1_draw))
  w[ok_e] <- exp(-0.5 * ((e_hat[ok_e] - e_target) / e_sd)^2)
  
  if (all(w <= 0) || !any(is.finite(w))) {
    if (any(is.finite(e_hat))) {
      j <- which.min(abs(e_hat - e_target))
      return(list(phi1_fixed=phi1_draw[j], e_implied=e_hat[j], ok=FALSE, reason="no_mass_fallback_minabs"))
    }
    return(list(phi1_fixed=NA_real_, e_implied=NA_real_, ok=FALSE, reason="no_mass"))
  }
  
  probs <- w / sum(w, na.rm=TRUE)
  j <- sample.int(length(phi1_draw), size=1, prob=probs)
  
  ok <- is.finite(w[j]) && (w[j] >= w_min_ok)
  list(phi1_fixed=phi1_draw[j], e_implied=e_hat[j], ok=ok, reason=if (ok) "ok" else "low_weight_selected")
}

# ============================================================
# Early exit
# ============================================================
if (nrow(stage1_ok) == 0) {
  cat("Stage 2: no Stage 1 survivors. Writing empty Stage 2 outputs.\n")
  write_csv(stage2_schema, file.path(dir_stage2, "stage2_finance_scan.csv"))
  write_csv(tibble(reason="no_stage1_survivors", n=0), file.path(dir_stage2, "failure_reasons.csv"))
  quit(save="no")
}

# ============================================================
# GRID
# ============================================================
stage2_grid <- tidyr::crossing(
  cand_id = seq_len(nrow(stage1_ok)),
  rF = rF_grid,
  psi = psi_grid,
  phi2 = phi2_grid
)

# ============================================================
# MAIN
# ============================================================
stage2_res <- purrr::pmap_dfr(stage2_grid, function(cand_id, rF, psi, phi2) {
  
  row <- stage1_ok[cand_id, ]
  
  p <- par_base
  p$kappa_max <- safe_num(row$kappa_max)
  
  fatal <- character()
  soft  <- character()
  
  r_star <- safe_num(row$r_star)
  d_star <- safe_num(row$d_star)
  g_n    <- safe_num(row$g_n)
  
  # --- lambda ---
  lambda_star <- lambda_fun(r_star, rF=rF, psi=psi)
  if (!is.finite(lambda_star)) fatal <- c(fatal, "nonfinite_lambda_star")
  
  # saturation is usually a warning, not a death sentence
  lam_eps_soft <- 1e-8
  if (is.finite(lambda_star) && (lambda_star < lam_eps_soft || lambda_star > 1 - lam_eps_soft)) {
    soft <- c(soft, "lambda_saturation")
  }
  
  # --- f* ---
  f_star <- NA_real_
  if (is.finite(lambda_star) && is.finite(r_star) && is.finite(g_n) && g_n > 0) {
    if (lambda_star >= 1) {
      fatal <- c(fatal, "lambda_ge_1")
    } else if (lambda_star <= 0) {
      fatal <- c(fatal, "lambda_le_0")
    } else {
      iotaF_star <- r_star * (lambda_star / (1 - lambda_star))
      f_star <- iotaF_star / g_n
    }
  } else {
    fatal <- c(fatal, "bad_inputs_for_f_star")
  }
  
  if (!is.finite(f_star)) fatal <- c(fatal, "nonfinite_f_star")
  if (is.finite(f_star) && f_star < 0) soft <- c(soft, "f_star_negative")  # often indicates region to avoid, but keep as soft
  
  # --- Z* ---
  Z_star <- NA_real_
  if (is.finite(d_star) && is.finite(f_star)) {
    Z_star <- Z_fun(d_star, f_star, p)
  } else {
    fatal <- c(fatal, "bad_inputs_for_Z")
  }
  if (!is.finite(Z_star)) fatal <- c(fatal, "nonfinite_Z_star")
  
  # --- phi1 calibration (fixed) ---
  phi1_fixed <- NA_real_
  e_implied  <- NA_real_
  phi1_calib_ok <- FALSE
  phi1_calib_reason <- "not_run"
  
  if (is.finite(Z_star) && is.finite(phi2)) {
    A <- (p$alpha - p$phi0 + phi2 * Z_star)
    cal <- calibrate_phi1_once(A, phi1_prior_min, phi1_prior_max, e_center, e_sd, e_min_ok, e_max_ok, draws=phi1_draws)
    
    phi1_fixed <- cal$phi1_fixed
    e_implied  <- cal$e_implied
    phi1_calib_ok <- isTRUE(cal$ok)
    phi1_calib_reason <- as.character(cal$reason)
    
    if (!is.finite(phi1_fixed)) fatal <- c(fatal, "nonfinite_phi1_fixed")
    if (!is.finite(e_implied))  fatal <- c(fatal, "nonfinite_e_implied")
    if (!phi1_calib_ok) soft <- c(soft, paste0("phi1_calib_", phi1_calib_reason))
  } else {
    fatal <- c(fatal, "bad_inputs_for_phi1_calib")
    phi1_calib_reason <- "bad_inputs_for_phi1_calib"
  }
  
  # econ_ok based ONLY on fatal
  econ_ok <- (length(fatal) == 0)
  
  reasons <- c(fatal, paste0("soft:", soft))
  reasons <- reasons[!is.na(reasons) & nzchar(reasons)]
  reason  <- if (length(reasons) == 0) "ok" else paste(unique(reasons), collapse="|")
  
  tibble(
    try_id = row$try_id,
    cand_id = cand_id,
    sigma = row$sigma, g_n = row$g_n, i = row$i, delta = row$delta,
    kappa_min = row$kappa_min, kappa_max = row$kappa_max,
    kappa_star = row$kappa_star, r_star = row$r_star, d_star = row$d_star, omega_star = row$omega_star,
    
    rF = rF, psi = psi, phi2 = phi2,
    lambda_star = lambda_star, f_star = f_star, Z_star = Z_star,
    
    phi1_fixed = phi1_fixed,
    e_implied  = e_implied,
    phi1_calib_ok = phi1_calib_ok,
    phi1_calib_reason = phi1_calib_reason,
    
    econ_ok = econ_ok,
    reason = reason
  )
}) %>% bind_rows(stage2_schema)

cat("Stage 2 counts:\n")
cat("  n_try     =", nrow(stage2_res), "\n")
cat("  n_econ_ok =", sum(stage2_res$econ_ok, na.rm = TRUE), "\n")
cat("  n_phi1_calib_ok =", sum(stage2_res$phi1_calib_ok, na.rm = TRUE), "\n")

write_csv(stage2_res, file.path(dir_stage2, "stage2_finance_scan.csv"))

stage2_reason_counts <- stage2_res %>%
  mutate(reason = ifelse(is.na(reason) | reason=="", "ok", reason)) %>%
  separate_rows(reason, sep="\\|") %>%
  count(reason, sort=TRUE)
write_csv(stage2_reason_counts, file.path(dir_stage2, "failure_reasons.csv"))

# quick plots (optional)
p1 <- ggplot(stage2_res, aes(x=lambda_star, fill=econ_ok)) +
  geom_histogram(bins=40, alpha=0.7) +
  facet_wrap(~kappa_max) +
  theme_minimal() + labs(title="Stage 2: lambda* distribution", x="lambda*", y="count")
ggsave(file.path(dir_stage2, "lambda_dist.png"), p1, width=10, height=6, dpi=160)
