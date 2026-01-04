# ============================================================
# 01_stage1_backbone.R  (V2: ROBUST + CONSISTENT WITH STAGE 3–5)
# Stage 1: Backbone feasibility (with kappa_max scan)
# Outputs -> outputs/wealth_goodwin/stage1_backbone/
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(tibble)
})

# ----------------------------
# Load globals + paths (idempotent, standalone-safe)
# ----------------------------
if (!exists("dirs", inherits = TRUE) || !exists("par_base", inherits = TRUE) ||
    !exists("sigma_grid", inherits = TRUE) || !exists("gn_grid", inherits = TRUE) ||
    !exists("i_grid", inherits = TRUE) || !exists("delta_grid", inherits = TRUE) ||
    !exists("kappa_max_grid", inherits = TRUE)) {
  
  get_script_dir_local <- function() {
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
      p <- rstudioapi::getActiveDocumentContext()$path
      if (!is.null(p) && nzchar(p)) return(normalizePath(dirname(p), winslash = "/", mustWork = FALSE))
    }
    args <- commandArgs(trailingOnly = FALSE)
    hit <- grep("^--file=", args, value = TRUE)
    if (length(hit) > 0) {
      p <- sub("^--file=", "", hit[1])
      return(normalizePath(dirname(p), winslash = "/", mustWork = FALSE))
    }
    normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  }
  
  if (file.exists("00_globals_utils.R")) {
    source("00_globals_utils.R", local = FALSE)
  } else {
    code_dir <- get_script_dir_local()
    globals_path <- file.path(code_dir, "00_globals_utils.R")
    stopifnot(file.exists(globals_path))
    source(globals_path, local = FALSE)
  }
}

# Fallback helper if not defined in globals
if (!exists("collapse_reasons", inherits = TRUE)) {
  collapse_reasons <- function(x) {
    x <- x[is.finite(match(x, x))] # keeps order, drops NA-ish
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x) == 0) return("ok")
    paste(unique(x), collapse = "|")
  }
}

# Stage 1 output dir (canonical)
dir_stage1 <- dirs$stage1_backbone %||% file.path(dirs$out_root %||% "outputs/wealth_goodwin", "stage1_backbone")
dir.create(dir_stage1, showWarnings = FALSE, recursive = TRUE)

cat("\n====================\nStage 1: Backbone (with kappa_max scan) [V2]\n====================\n")
cat("Output dir:", dir_stage1, "\n")

# ----------------------------
# Grid
# ----------------------------
stage1_grid <- tidyr::crossing(
  sigma     = sigma_grid,
  g_n       = gn_grid,
  i         = i_grid,
  delta     = delta_grid,
  kappa_max = kappa_max_grid
) %>%
  mutate(try_id = row_number())

# Empty schema for type stability
stage1_schema <- tibble(
  try_id = integer(),
  sigma = double(), g_n = double(), i = double(), delta = double(),
  kappa_min = double(), kappa_max = double(),
  kappa_star = double(), r_star = double(),
  d_star = double(), omega_star = double(),
  backbone_ok = logical(),
  reason = character()
)

# ----------------------------
# Compute feasibility
# ----------------------------
kappa_bound_eps <- if (exists("cfg", inherits = TRUE) && !is.null(cfg$clamps$kappa_eps)) cfg$clamps$kappa_eps else 1e-10

stage1_res <- purrr::pmap_dfr(stage1_grid, function(sigma, g_n, i, delta, kappa_max, try_id) {
  
  p <- par_base
  p$kappa_max <- kappa_max
  
  reasons <- character()
  
  # κ* = σ (g_n + δ)
  kappa_star <- sigma * (g_n + delta)
  if (!is.finite(kappa_star)) reasons <- c(reasons, "nonfinite_kappa_star")
  
  # Must lie strictly inside logistic bounds (avoid inversion instability)
  if (is.finite(kappa_star) &&
      (kappa_star <= p$kappa_min + kappa_bound_eps || kappa_star >= p$kappa_max - kappa_bound_eps)) {
    reasons <- c(reasons, "kappa_star_out_of_bounds")
  }
  
  r_star <- NA_real_
  d_star <- NA_real_
  omega_star <- NA_real_
  
  if (!("kappa_star_out_of_bounds" %in% reasons) && is.finite(kappa_star)) {
    
    # r* from inverse kappa mapping
    r_star <- kappa_inv(kappa_star, p)
    if (!is.finite(r_star)) reasons <- c(reasons, "nonfinite_r_star")
    
    # d* = (κ* - σ r*) / g_n
    if (is.finite(r_star) && is.finite(g_n) && g_n > 0) {
      d_star <- (kappa_star - sigma * r_star) / g_n
    } else {
      reasons <- c(reasons, "bad_gn_for_d_star")
    }
    
    if (!is.finite(d_star)) reasons <- c(reasons, "nonfinite_d_star")
    if (is.finite(d_star) && d_star < 0) reasons <- c(reasons, "d_star_negative")
    
    # ω* = 1 - i d* - σ r*
    if (is.finite(d_star) && is.finite(r_star)) {
      omega_star <- 1 - i * d_star - sigma * r_star
    }
    if (!is.finite(omega_star)) reasons <- c(reasons, "nonfinite_omega_star")
    
    # Feasible ω* range
    if (is.finite(omega_star) && (omega_star <= 0 || omega_star >= 1)) reasons <- c(reasons, "omega_outside_0_1")
    if (exists("omega_band", inherits = TRUE) && is.numeric(omega_band) && length(omega_band) == 2) {
      if (is.finite(omega_star) && !(omega_star >= omega_band[1] && omega_star <= omega_band[2])) {
        reasons <- c(reasons, "omega_outside_target_band")
      }
    }
  }
  
  backbone_ok <- (length(reasons) == 0)
  
  tibble(
    try_id = try_id,
    sigma = sigma, g_n = g_n, i = i, delta = delta,
    kappa_min = p$kappa_min, kappa_max = p$kappa_max,
    kappa_star = kappa_star,
    r_star = r_star,
    d_star = d_star,
    omega_star = omega_star,
    backbone_ok = backbone_ok,
    reason = collapse_reasons(reasons)
  )
}) %>%
  bind_rows(stage1_schema)

cat("Stage 1 counts:\n")
cat("  n_try         =", nrow(stage1_res), "\n")
cat("  n_backbone_ok =", sum(stage1_res$backbone_ok, na.rm = TRUE), "\n")

# ----------------------------
# Diagnostics (optional)
# ----------------------------
RUN_STAGE1_DIAG <- TRUE
if (isTRUE(RUN_STAGE1_DIAG) && exists("run_stage1_diagnostics", mode = "function", inherits = TRUE)) {
  diag <- run_stage1_diagnostics(stage1_res)
  readr::write_csv(diag$tab_condition, file.path(dir_stage1, "diag_condition_counts.csv"))
  readr::write_csv(diag$gap_stats,     file.path(dir_stage1, "diag_gap_stats.csv"))
  readr::write_csv(diag$term_stats,    file.path(dir_stage1, "diag_term_stats.csv"))
  readr::write_csv(diag$omega_dir,     file.path(dir_stage1, "diag_omega_dir.csv"))
  readr::write_csv(diag$reasons_all,   file.path(dir_stage1, "diag_reasons_all.csv"))
  readr::write_csv(diag$reasons_second_order, file.path(dir_stage1, "diag_reasons_second_order.csv"))
}

# ----------------------------
# Exports
# ----------------------------
readr::write_csv(stage1_res, file.path(dir_stage1, "stage1_backbone.csv"))

stage1_reason_counts <- stage1_res %>%
  tidyr::separate_rows(reason, sep = "\\|") %>%
  count(reason, sort = TRUE)

readr::write_csv(stage1_reason_counts, file.path(dir_stage1, "failure_reasons.csv"))

stage1_reason_by_kmax <- stage1_res %>%
  tidyr::separate_rows(reason, sep = "\\|") %>%
  count(kappa_max, reason, sort = TRUE)

readr::write_csv(stage1_reason_by_kmax, file.path(dir_stage1, "failure_reasons_by_kappa_max.csv"))

stage1_backbone_rate <- stage1_res %>%
  group_by(kappa_max) %>%
  summarise(
    n = n(),
    n_ok = sum(backbone_ok, na.rm = TRUE),
    ok_rate = n_ok / n,
    .groups = "drop"
  )

readr::write_csv(stage1_backbone_rate, file.path(dir_stage1, "backbone_rate_by_kappa_max.csv"))

# ----------------------------
# Plots (keep, but cheap)
# ----------------------------
p1 <- ggplot(stage1_res, aes(x = r_star, y = omega_star, color = backbone_ok)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~kappa_max) +
  labs(title = "Stage 1: omega* vs r* (by kappa_max)", x = "r*", y = "omega*") +
  theme_minimal()

p2 <- ggplot(stage1_res, aes(x = r_star, y = d_star, color = backbone_ok)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~kappa_max) +
  labs(title = "Stage 1: d* vs r* (by kappa_max)", x = "r*", y = "d*") +
  theme_minimal()

p3 <- ggplot(stage1_res, aes(x = omega_star, fill = backbone_ok)) +
  geom_histogram(bins = 40, alpha = 0.7) +
  { if (exists("omega_band", inherits = TRUE) && is.numeric(omega_band)) geom_vline(xintercept = omega_band, linetype = "dashed") } +
  facet_wrap(~kappa_max) +
  labs(title = "Stage 1: omega* histogram (by kappa_max)", x = "omega*", y = "count") +
  theme_minimal()

p4 <- ggplot(stage1_res, aes(x = d_star, fill = backbone_ok)) +
  geom_histogram(bins = 40, alpha = 0.7) +
  facet_wrap(~kappa_max) +
  labs(title = "Stage 1: d* histogram (by kappa_max)", x = "d*", y = "count") +
  theme_minimal()

p4b <- ggplot(stage1_backbone_rate, aes(x = factor(kappa_max), y = ok_rate)) +
  geom_col() +
  labs(title = "Stage 1: backbone OK rate by kappa_max", x = "kappa_max", y = "OK rate") +
  theme_minimal()

ggsave(file.path(dir_stage1, "omega_vs_r.png"), p1, width = 10, height = 6, dpi = 160)
ggsave(file.path(dir_stage1, "d_vs_r.png"),     p2, width = 10, height = 6, dpi = 160)
ggsave(file.path(dir_stage1, "hist_omega.png"), p3, width = 10, height = 6, dpi = 160)
ggsave(file.path(dir_stage1, "hist_d.png"),     p4, width = 10, height = 6, dpi = 160)
ggsave(file.path(dir_stage1, "backbone_rate_by_kappa_max.png"), p4b, width = 8, height = 5, dpi = 160)

cat("\nStage 1 top binding constraints:\n")
print(stage1_reason_counts %>% slice_head(n = 8))
