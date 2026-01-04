# ============================================================
# 03_stage3_stability_hopf_FAST.R  (STANDALONE + COARSE-TO-FINE)
# Stage 3: analytic stability + Hopf root scan (fast)
#
# Fixes:
#  - Robust path resolution for Stage 2 input (works even if wd is "wrong")
#  - Forces analytic core to load (NO passthrough; tripwire hard-stops)
#  - Maps phi1_fixed -> phi1 (PHI1 LOCK consistent with your World 1 rule)
#  - Baseline gate BEFORE hopf scan
#  - Coarse-to-fine hopf scanning (big speed-up)
#
# Outputs (canonical):
#   outputs/wealth_goodwin/stage3_stability/stage3_candidates.csv
#   outputs/wealth_goodwin/stage3_stability/manifest_stage3.csv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(tibble)
  library(stringr)
})

cat("\n====================\n")
cat("Stage 3: Stability + Hopf scan (FAST)\n")
cat("====================\n")

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

# ----------------------------
# Load globals/dirs if available (standalone-safe)
# ----------------------------
if (!exists("dirs", inherits = TRUE)) {
  if (file.exists("00_globals_utils.R")) {
    source("00_globals_utils.R", local = FALSE)
  } else {
    # Minimal fallback if you run from a random working directory
    # (because humans love chaos)
    dirs <- list(
      out_root = "outputs/wealth_goodwin",
      stage2_finance = file.path("outputs/wealth_goodwin", "stage2_finance"),
      stage2 = file.path("outputs/wealth_goodwin", "stage2"),
      stage3_stability = file.path("outputs/wealth_goodwin", "stage3_stability")
    )
  }
}

# ----------------------------
# Source analytic core (HARD TRIPWIRE: must exist)
# ----------------------------
core_candidates <- c(
  "00b_wealthIneq_analytic_core.R",
  file.path("R", "00b_wealthIneq_analytic_core.R"),
  file.path("R", "hooks", "00b_wealthIneq_analytic_core.R"),
  file.path("codes", "goodwin_wealth", "00b_wealthIneq_analytic_core.R")
)
core_file <- core_candidates[file.exists(core_candidates)][1]

if (is.na(core_file) || !nzchar(core_file)) {
  stop("Stage 3: cannot find analytic core 00b_wealthIneq_analytic_core.R.\nLooked for:\n  - ",
       paste(core_candidates, collapse = "\n  - "))
}

source(core_file, local = FALSE)

# Tripwire: refuse to continue in passthrough land
stopifnot(exists("compute_at_rF", mode = "function", inherits = TRUE))

cat("Stage 3: using analytic core:\n  ", core_file, "\n", sep = "")

# ----------------------------
# Robust Stage 2 input resolution
# ----------------------------
stage2_candidates <- c(
  # canonical
  file.path(dirs$stage2_finance %||% file.path(dirs$out_root, "stage2_finance"), "stage2_candidates.csv"),
  file.path(dirs$stage2 %||% file.path(dirs$out_root, "stage2"), "stage2_candidates.csv"),
  
  # common relative
  "outputs/wealth_goodwin/stage2_finance/stage2_candidates.csv",
  "outputs/wealth_goodwin/stage2/stage2_candidates.csv",
  
  # last resort: search under outputs/wealth_goodwin
  list.files("outputs/wealth_goodwin", pattern = "stage2_candidates\\.csv$", recursive = TRUE, full.names = TRUE)
) %>% unlist() %>% unique()

in2 <- stage2_candidates[file.exists(stage2_candidates)][1]
if (is.na(in2) || !nzchar(in2)) {
  stop("Stage 3: cannot find Stage 2 candidates.\nLooked for:\n  - ",
       paste(stage2_candidates, collapse = "\n  - "))
}
cat("Stage 3: using Stage 2 input:\n  ", in2, "\n", sep = "")

df2 <- readr::read_csv(in2, show_col_types = FALSE)

if (nrow(df2) == 0) stop("Stage 3: Stage 2 input is empty: ", in2)

# Optional Stage2 gating if econ_ok exists
if ("econ_ok" %in% names(df2)) {
  before <- nrow(df2)
  df2 <- df2 %>% filter(.data$econ_ok %in% c(TRUE, 1))
  cat("Stage 3: econ_ok filter kept ", nrow(df2), " / ", before, " rows\n", sep = "")
}

# ----------------------------
# Output dir
# ----------------------------
out_dir <- dirs$stage3_stability %||% file.path(dirs$out_root, "stage3_stability") %||% "outputs/wealth_goodwin/stage3_stability"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_candidates <- file.path(out_dir, "stage3_candidates.csv")
out_manifest   <- file.path(out_dir, "manifest_stage3.csv")

# ----------------------------
# Paper-grade admissibility bands (same idea as Stage 4)
# ----------------------------
B <- list(
  e_min = 0.20, e_max = 0.98,
  w_min = 0.40, w_max = 0.95,
  d_min = 0.00, d_max = 3.00
)

# ----------------------------
# Hopf scan settings (coarse-to-fine)
# ----------------------------
hopf_grid_coarse    <- seq(0.005, 0.20, by = 0.005)
hopf_refine_span    <- 0.01
hopf_grid_fine_step <- 0.001

# ----------------------------
# Helper: linearized root from sign change
# ----------------------------
roots_from_scan <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 2) return(numeric(0))
  
  s <- y[-length(y)] * y[-1]
  idx <- which(s <= 0)
  
  if (length(idx) == 0) return(numeric(0))
  
  map_dbl(idx, function(i) {
    x1 <- x[i]; x2 <- x[i+1]
    y1 <- y[i]; y2 <- y[i+1]
    if (!is.finite(y1) || !is.finite(y2)) return(NA_real_)
    if (abs(y2 - y1) < 1e-12) return(mean(c(x1, x2)))
    x1 - y1 * (x2 - x1) / (y2 - y1)
  }) %>% .[is.finite(.)] %>% sort()
}

# ----------------------------
# Core per-row evaluator
# ----------------------------
eval_one <- function(df_row) {
  row <- as.list(df_row)
  
  # ----- schema harmonization -----
  # PHI1 LOCK: analytic core expects row$phi1
  if (is.null(row$phi1) || !is.finite(as.numeric(row$phi1))) {
    if (!is.null(row$phi1_fixed) && is.finite(as.numeric(row$phi1_fixed))) {
      row$phi1 <- as.numeric(row$phi1_fixed)
    }
  }
  
  # baseline rF for evaluation
  rF0 <- as.numeric(row$rF %||% NA_real_)
  
  if (!is.finite(rF0)) {
    return(tibble(
      try_id = row$try_id %||% NA_integer_,
      cand_id = row$cand_id %||% NA_integer_,
      ok = FALSE,
      reason = "missing_rF",
      n_roots = 0L, rF_root_med = NA_real_, rF_root_p25 = NA_real_, rF_root_p75 = NA_real_,
      has_hopf = FALSE
    ))
  }
  
  # ----- baseline analytic evaluation at rF0 -----
  base <- tryCatch(compute_at_rF(row, rFj = rF0),
                   error = function(e) list(ok = FALSE, reason = paste0("compute_error:", conditionMessage(e))))
  
  if (!isTRUE(base$ok)) {
    return(tibble(
      try_id = row$try_id %||% NA_integer_,
      cand_id = row$cand_id %||% NA_integer_,
      sigma = row$sigma %||% NA_real_,
      g_n = row$g_n %||% NA_real_,
      i = row$i %||% NA_real_,
      delta = row$delta %||% NA_real_,
      kappa_min = row$kappa_min %||% NA_real_,
      kappa_max = row$kappa_max %||% NA_real_,
      
      r_star = row$r_star %||% NA_real_,
      d_star = row$d_star %||% NA_real_,
      omega_star = row$omega_star %||% NA_real_,
      
      rF = rF0, psi = row$psi %||% NA_real_, phi2 = row$phi2 %||% NA_real_,
      phi1_fixed = row$phi1 %||% row$phi1_fixed %||% NA_real_,
      e_implied = NA_real_,
      
      stable = NA, RH_ok = NA, has_complex = NA, hopf_val = NA_real_,
      ok = FALSE,
      reason = base$reason %||% "compute_failed",
      
      n_roots = 0L, rF_root_med = NA_real_, rF_root_p25 = NA_real_, rF_root_p75 = NA_real_,
      has_hopf = FALSE
    ))
  }
  
  # baseline gate before scan
  gate_interior0 <-
    is.finite(base$e_star) && base$e_star >= B$e_min && base$e_star <= B$e_max &&
    is.finite(base$omega_star) && base$omega_star >= B$w_min && base$omega_star <= B$w_max &&
    is.finite(base$d_star) && base$d_star >= B$d_min && base$d_star <= B$d_max
  
  if (!isTRUE(gate_interior0)) {
    return(tibble(
      try_id = row$try_id %||% NA_integer_,
      cand_id = row$cand_id %||% NA_integer_,
      sigma = row$sigma %||% NA_real_,
      g_n = row$g_n %||% NA_real_,
      i = row$i %||% NA_real_,
      delta = row$delta %||% NA_real_,
      kappa_min = row$kappa_min %||% NA_real_,
      kappa_max = row$kappa_max %||% NA_real_,
      
      r_star = base$r %||% row$r_star %||% NA_real_,
      d_star = base$d_star,
      omega_star = base$omega_star,
      
      rF = rF0, psi = row$psi %||% NA_real_, phi2 = row$phi2 %||% NA_real_,
      phi1_fixed = base$phi1 %||% row$phi1_fixed %||% NA_real_,
      e_implied = base$e_star,
      
      stable = isTRUE(base$stable),
      RH_ok = isTRUE(base$RH_ok),
      has_complex = isTRUE(base$has_complex),
      hopf_val = base$hopf_val,
      
      ok = TRUE,
      reason = "ok_baseline_but_reject_interior_bounds",
      
      n_roots = 0L, rF_root_med = NA_real_, rF_root_p25 = NA_real_, rF_root_p75 = NA_real_,
      has_hopf = FALSE
    ))
  }
  
  # ----- coarse scan -----
  scanC <- purrr::map_dfr(hopf_grid_coarse, function(rFj) {
    out <- tryCatch(compute_at_rF(row, rFj = rFj),
                    error = function(e) list(ok = FALSE))
    tibble(rFj = rFj, hopf_val = if (isTRUE(out$ok)) out$hopf_val else NA_real_)
  })
  
  roots <- numeric(0)
  hv <- scanC$hopf_val
  xx <- scanC$rFj
  okv <- is.finite(hv) & is.finite(xx)
  hv <- hv[okv]; xx <- xx[okv]
  
  brackets <- integer(0)
  if (length(xx) >= 2) {
    s <- hv[-length(hv)] * hv[-1]
    brackets <- which(s <= 0)
  }
  
  # ----- refine only around brackets -----
  if (length(brackets) > 0) {
    for (j in brackets) {
      mid <- mean(c(xx[j], xx[j+1]))
      fine <- seq(max(0.0, mid - hopf_refine_span),
                  min(1.0, mid + hopf_refine_span),
                  by = hopf_grid_fine_step)
      
      scanF <- purrr::map_dfr(fine, function(rFj) {
        out <- tryCatch(compute_at_rF(row, rFj = rFj),
                        error = function(e) list(ok = FALSE))
        tibble(rFj = rFj, hopf_val = if (isTRUE(out$ok)) out$hopf_val else NA_real_)
      })
      
      roots <- c(roots, roots_from_scan(scanF$rFj, scanF$hopf_val))
    }
  }
  
  roots <- sort(unique(round(roots, 6)))
  n_roots <- length(roots)
  
  tibble(
    try_id = row$try_id %||% NA_integer_,
    cand_id = row$cand_id %||% NA_integer_,
    
    sigma = row$sigma %||% NA_real_,
    g_n = row$g_n %||% NA_real_,
    i = row$i %||% NA_real_,
    delta = row$delta %||% NA_real_,
    kappa_min = row$kappa_min %||% NA_real_,
    kappa_max = row$kappa_max %||% NA_real_,
    
    r_star = base$r %||% row$r_star %||% NA_real_,
    d_star = base$d_star,
    omega_star = base$omega_star,
    
    rF = rF0, psi = row$psi %||% NA_real_, phi2 = row$phi2 %||% NA_real_,
    
    lambda_star = base$lam %||% NA_real_,
    f_star = base$f %||% NA_real_,
    Z_star = base$Z %||% NA_real_,
    
    phi1_fixed = base$phi1 %||% row$phi1_fixed %||% NA_real_,
    e_implied = base$e_star,
    
    stable = isTRUE(base$stable),
    RH_ok = isTRUE(base$RH_ok),
    has_complex = isTRUE(base$has_complex),
    hopf_val = base$hopf_val %||% NA_real_,
    
    ok = TRUE,
    reason = "ok",
    
    n_roots = as.integer(n_roots),
    rF_root_med = if (n_roots > 0) median(roots) else NA_real_,
    rF_root_p25 = if (n_roots > 0) as.numeric(stats::quantile(roots, 0.25)) else NA_real_,
    rF_root_p75 = if (n_roots > 0) as.numeric(stats::quantile(roots, 0.75)) else NA_real_,
    has_hopf = (n_roots > 0)
  )
}

# ----------------------------
# Run Stage 3
# ----------------------------
row_lists <- split(df2, seq_len(nrow(df2)))

cat("Stage 3: evaluating ", length(row_lists), " candidates...\n", sep = "")

stage3 <- purrr::map_dfr(row_lists, eval_one)

readr::write_csv(stage3, out_candidates)

manifest <- tibble(
  file = c(out_candidates),
  exists = file.exists(c(out_candidates)),
  n_rows = c(nrow(stage3))
)
readr::write_csv(manifest, out_manifest)

cat("\nStage 3 wrote:\n  - ", out_candidates, "\n", sep = "")
cat("Stage 3 done.\n")
