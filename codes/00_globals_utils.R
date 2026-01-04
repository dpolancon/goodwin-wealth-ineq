# ============================
# 00_globals_utils.R (PATHS + GLOBALS)  [AUDITED + CONSOLIDATED]
# Repo layout:
#   .../goodwin_model/codes/goodwin_wealth
# Outputs:
#   .../goodwin_model/outputs/wealth_goodwin/{stage...}
# ============================

options(stringsAsFactors = FALSE)

# ----------------------------
# Null-coalescing helper (used everywhere)
# ----------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ----------------------------
# Path helpers
# ----------------------------
norm_path <- function(p) normalizePath(p, winslash = "/", mustWork = FALSE)

dir_create_safe <- function(p) {
  dir.create(p, showWarnings = FALSE, recursive = TRUE)
  invisible(p)
}

# Robust: find directory of the current script (RStudio or Rscript)
get_script_dir <- function() {
  # RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- rstudioapi::getActiveDocumentContext()$path
    if (!is.null(p) && nzchar(p)) return(norm_path(dirname(p)))
  }
  # Rscript --file=...
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) > 0) {
    p <- sub("^--file=", "", hit[1])
    return(norm_path(dirname(p)))
  }
  # Fallback
  norm_path(getwd())
}

# ----------------------------
# Root paths
# ----------------------------
CODE_ROOT <- get_script_dir()  # .../codes/goodwin_wealth
REPO_ROOT <- norm_path(file.path(CODE_ROOT, "..", ".."))  # .../goodwin_model
OUT_ROOT  <- norm_path(file.path(REPO_ROOT, "outputs", "wealth_goodwin"))

# ----------------------------
# Output folders by stage
# (Includes aliases to avoid downstream breakage)
# ----------------------------
dirs <- list(
  code_root = CODE_ROOT,
  repo_root = REPO_ROOT,
  out_root  = OUT_ROOT,
  
  # canonical stage folders
  stage1_backbone  = file.path(OUT_ROOT, "stage1_backbone"),
  stage2_finance   = file.path(OUT_ROOT, "stage2_finance"),
  stage3_stability = file.path(OUT_ROOT, "stage3_stability"),
  stage4_scoring   = file.path(OUT_ROOT, "stage4_scoring"),
  stage5_simulation= file.path(OUT_ROOT, "stage5_simulation"),
  
  # short aliases (compat)
  stage1 = file.path(OUT_ROOT, "stage1_backbone"),
  stage2 = file.path(OUT_ROOT, "stage2_finance"),
  stage3 = file.path(OUT_ROOT, "stage3_stability"),
  stage4 = file.path(OUT_ROOT, "stage4_scoring"),
  stage5 = file.path(OUT_ROOT, "stage5_simulation"),
  
  # Stage 5 subfolders
  stage5A = file.path(OUT_ROOT, "stage5_simulation", "5A"),
  stage5B = file.path(OUT_ROOT, "stage5_simulation", "5B"),
  stage5C = file.path(OUT_ROOT, "stage5_simulation", "5C"),
  stage5D = file.path(OUT_ROOT, "stage5_simulation", "5D"),
  stage5E = file.path(OUT_ROOT, "stage5_simulation", "5E"),
  
  # analysis sink (optional)
  analysis = file.path(OUT_ROOT, "analysis")
)

# Create output folders (safe)
invisible(lapply(
  dirs[names(dirs) %in% c(
    "stage1_backbone","stage2_finance","stage3_stability","stage4_scoring","stage5_simulation",
    "stage1","stage2","stage3","stage4","stage5",
    "stage5A","stage5B","stage5C","stage5D","stage5E","analysis"
  )],
  dir_create_safe
))

cat("\n[PATHS]\n",
    "CODE_ROOT: ", dirs$code_root, "\n",
    "REPO_ROOT: ", dirs$repo_root, "\n",
    "OUT_ROOT : ", dirs$out_root,  "\n", sep = "")

# ============================================================
# GLOBAL TARGETS / GRIDS / BASE PARAMS
# Keep as globals for compatibility + also store in cfg
# ============================================================

# ----------------------------
# Targets
# ----------------------------
omega_target <- 0.65
omega_band   <- c(0.55, 0.70)

e_target <- 0.94
e_band   <- c(0.88, 0.95)

# ----------------------------
# Minimal grids (Stage 1 backbone)
# ----------------------------
sigma_grid <- c(2.25, 2.35, 2.50)
gn_grid    <- seq(0.05, 0.08, length.out = 4)
i_grid     <- seq(0.025, 0.040, length.out = 4)
delta_grid <- seq(0.03, 0.06,  length.out = 4)

# scan kappa_max to avoid hard-choking feasibility
kappa_max_grid <- c(0.25, 0.35, 0.45)

# ----------------------------
# Stage 2 scans
# ----------------------------
rF_grid   <- seq(0.02, 0.12, length.out = 9)
psi_grid  <- c(5, 10, 15, 20)
phi2_grid <- seq(0.0, 3.0, length.out = 7)

# Hopf scan along rF around each candidate (Stage 3)
hopf_rF_span <- 0.05
hopf_n_grid  <- 31

# ----------------------------
# Base parameters
# IMPORTANT: only kappa_max is scanned; everything else stays the same
# ----------------------------
par_base <- list(
  # kappa(r): logistic bounds + center/steepness
  kappa_min = 0.02,
  kappa_max = 0.25,   # baseline (overridden in scan)
  kappa0    = 0.10,
  kappa1    = 30.0,
  
  # Z(d,f)
  phi3 = 8.0,
  phi4 = 1.0,
  
  # omega dynamics
  phi0  = -0.02,
  alpha = 0.02,
  
  # phi1 bounds (Stage 2 calibration; fixed for Stage 3–5)
  phi1_min = 0.10,
  phi1_max = 5.00
)

# ============================================================
# NUMERIC SAFETY / COMMON HELPERS
# ============================================================

# Stable logistic (prevents overflow; faster than homebrew exp)
logistic <- function(x) plogis(x)

# Make sure we can safely coerce scalars
as_num1 <- function(x, default = NA_real_) {
  v <- suppressWarnings(as.numeric(x))
  if (length(v) == 0 || is.na(v[1]) || !is.finite(v[1])) return(default)
  v[1]
}

# ------------------------------------------------------------
# Core model functions
# ------------------------------------------------------------
kappa_fun <- function(r, p) {
  p$kappa_min + (p$kappa_max - p$kappa_min) * logistic(p$kappa1 * (r - p$kappa0))
}

kappa_prime <- function(r, p) {
  s <- logistic(p$kappa1 * (r - p$kappa0))
  (p$kappa_max - p$kappa_min) * p$kappa1 * s * (1 - s)
}

kappa_inv <- function(kappa_star, p) {
  y <- (kappa_star - p$kappa_min) / (p$kappa_max - p$kappa_min)
  y <- pmin(pmax(y, 1e-12), 1 - 1e-12)
  p$kappa0 + (1 / p$kappa1) * log(y / (1 - y))
}

Z_fun <- function(d, f, p) {
  logistic(p$phi3 * ((d - 1) + p$phi4 * (f - 1)))
}

lambda_fun <- function(r, rF, psi) logistic(psi * (r - rF))

collapse_reasons <- function(reasons) {
  reasons <- unique(reasons[!is.na(reasons) & reasons != ""])
  if (length(reasons) == 0) "ok" else paste(reasons, collapse = "|")
}

# ============================================================
# Stage 1 diagnostics helper (kept as-is, but stable)
# ============================================================
run_stage1_diagnostics <- function(stage1_res, omega_band = NULL) {
  suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(tibble)
  })
  
  if (is.null(omega_band)) {
    omega_band <- get0("omega_band", envir = .GlobalEnv, inherits = FALSE)
    if (is.null(omega_band)) stop("Provide omega_band or define numeric omega_band in .GlobalEnv.")
  }
  
  omega_band <- as.numeric(omega_band)
  if (length(omega_band) != 2 || any(!is.finite(omega_band))) {
    stop("omega_band must be numeric length-2 with finite values.")
  }
  
  tab_condition <- stage1_res %>%
    mutate(gnd = g_n + delta,
           condition = ifelse(is.finite(r_star) & is.finite(gnd), r_star <= gnd, NA)) %>%
    count(condition, backbone_ok, name = "n") %>%
    arrange(desc(n))
  
  gap_stats <- stage1_res %>%
    mutate(gnd = g_n + delta, gap = r_star - gnd) %>%
    summarise(
      n_finite = sum(is.finite(gap)),
      gap_p50 = median(gap, na.rm = TRUE),
      gap_p90 = quantile(gap, 0.90, na.rm = TRUE, names = FALSE),
      gap_p95 = quantile(gap, 0.95, na.rm = TRUE, names = FALSE)
    )
  
  term_stats <- stage1_res %>%
    filter(is.finite(omega_star), is.finite(d_star), is.finite(r_star)) %>%
    mutate(term_id = i * d_star,
           term_sr = sigma * r_star) %>%
    summarise(
      n = n(),
      med_id = median(term_id, na.rm = TRUE),
      med_sr = median(term_sr, na.rm = TRUE),
      med_omega = median(omega_star, na.rm = TRUE)
    )
  
  omega_dir <- stage1_res %>%
    mutate(omega_dir = case_when(
      !is.finite(omega_star) ~ "nonfinite",
      omega_star < omega_band[1] ~ "too_low",
      omega_star > omega_band[2] ~ "too_high",
      TRUE ~ "in_band"
    )) %>%
    count(omega_dir, backbone_ok, name = "n") %>%
    arrange(desc(n))
  
  reasons_all <- stage1_res %>%
    mutate(reason = ifelse(is.na(reason) | reason == "", "ok", reason)) %>%
    separate_rows(reason, sep = "\\|") %>%
    mutate(reason = ifelse(is.na(reason) | reason == "", "ok", reason)) %>%
    count(reason, sort = TRUE, name = "n")
  
  reasons_second_order <- stage1_res %>%
    mutate(gnd = g_n + delta,
           condition = ifelse(is.finite(r_star) & is.finite(gnd), r_star <= gnd, NA)) %>%
    filter(condition %in% TRUE, backbone_ok %in% FALSE) %>%
    mutate(reason = ifelse(is.na(reason) | reason == "", "ok", reason)) %>%
    separate_rows(reason, sep = "\\|") %>%
    mutate(reason = ifelse(is.na(reason) | reason == "", "ok", reason)) %>%
    count(reason, sort = TRUE, name = "n")
  
  list(
    tab_condition = tab_condition,
    gap_stats = gap_stats,
    term_stats = term_stats,
    omega_dir = omega_dir,
    reasons_all = reasons_all,
    reasons_second_order = reasons_second_order
  )
}

# ============================================================
# Stage 3 LaTeX helpers
# ============================================================
latex_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([%&#_{}$])", "\\\\\\1", x, perl = TRUE)
  x
}

df_to_latex_tabular <- function(df, caption = NULL, label = NULL, align = NULL) {
  stopifnot(is.data.frame(df))
  if (is.null(align)) align <- paste0(rep("l", ncol(df)), collapse = "")
  header <- paste(latex_escape(names(df)), collapse = " & ")
  body <- apply(df, 1, function(row) paste(latex_escape(row), collapse = " & "))
  body <- paste0(body, " \\\\")
  out <- c(
    "\\begin{table}[H]",
    "\\centering",
    if (!is.null(caption)) paste0("\\caption{", caption, "}"),
    if (!is.null(label))   paste0("\\label{", label, "}"),
    paste0("\\begin{tabular}{", align, "}"),
    "\\hline",
    paste0(header, " \\\\"),
    "\\hline",
    body,
    "\\hline",
    "\\end{tabular}",
    "\\end{table}"
  )
  paste(out, collapse = "\n")
}

# ============================================================
# cfg object (used by Stage 2–5 + analytic core)
# This is where your earlier file was missing critical fields.
# ============================================================
cfg <- list(
  paths = dirs,
  targets = list(
    omega_target = omega_target,
    omega_band   = omega_band,
    e_target     = e_target,
    e_band       = e_band
  ),
  grids = list(
    sigma_grid     = sigma_grid,
    gn_grid        = gn_grid,
    i_grid         = i_grid,
    delta_grid     = delta_grid,
    kappa_max_grid = kappa_max_grid,
    rF_grid        = rF_grid,
    psi_grid       = psi_grid,
    phi2_grid      = phi2_grid
  ),
  hopf = list(
    hopf_rF_span = hopf_rF_span,
    hopf_n_grid  = hopf_n_grid
  ),
  clamps = list(
    lambda_eps = 1e-6,  # used in Stage 2 + simulate_model clamp
    kappa_eps  = 1e-10
  ),
  seeds = list(
    stage2 = 123
  ),
  stage2 = list(
    top_per_backbone = 30L
  ),
  tolerances = list(
    eig_im_tol = 1e-8,
    lambda_fd_eps = 1e-6
  ),
  par_base = par_base
)

# ============================================================
# Source analytic core ONCE (and only once), AFTER cfg exists
# ============================================================
core_path <- file.path(CODE_ROOT, "00b_wealthIneq_analytic_core.R")
if (file.exists(core_path)) {
  source(core_path, local = FALSE)
} else {
  warning("Analytic core not found: ", core_path, " (Stage 3 will fallback/passthrough)")
}

if (exists("compute_at_rF", mode = "function")) {
  message("[OK] compute_at_rF is available (Stage 3 analytic mode enabled).")
} else {
  message("[WARN] compute_at_rF missing (Stage 3 will fallback).")
}
