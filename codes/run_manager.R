# ============================================================
# 99_run_manager_all_stages.R  (AUDITED + CONSISTENT)
# Run manager: stages 01 -> 05
#
# Assumes scripts live in same folder:
#   00_globals_utils.R   (ALREADY sources analytic core ONCE)
#   01_stage1_backbone.R
#   02_stage2_finance.R
#   03_stage3_stability_hopf.R
#   04_stage4_scoring_shortlist.R
#   05_stage5_pathwayA_light.R
# ============================================================

# ----------------------------
# Robust: find directory of the current script (RStudio or Rscript)
# ----------------------------
get_script_dir <- function() {
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

CODE_DIR <- get_script_dir()
oldwd <- getwd()
setwd(CODE_DIR)
on.exit(setwd(oldwd), add = TRUE)

cat("\n[RUN MANAGER]\nCODE_DIR =", CODE_DIR, "\n")

# ----------------------------
# Packages (optionally auto-install)
# ----------------------------
auto_install <- FALSE

req_pkgs <- c("dplyr","tidyr","purrr","readr","ggplot2","stringr","tibble","rstudioapi")
missing <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  cat("\nMissing packages:\n  ", paste(missing, collapse = ", "), "\n", sep = "")
  if (isTRUE(auto_install)) {
    install.packages(missing, dependencies = TRUE)
  } else {
    stop("Missing required packages. Install them or set auto_install=TRUE.")
  }
}

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
# Run controls
# ----------------------------
RUN_STAGE1 <- TRUE
RUN_STAGE2 <- TRUE
RUN_STAGE3 <- TRUE
RUN_STAGE4 <- TRUE
RUN_STAGE5 <- TRUE

# If TRUE, skip running a stage if its key outputs already exist
SKIP_IF_EXISTS <- FALSE

# If TRUE, Stage 3 hard-fails when compute_at_rF is missing
REQUIRE_ANALYTIC_HOPF <- TRUE

# ----------------------------
# Source Stage 00 (paths + globals + cfg + analytic core ONCE)
# ----------------------------
globals_path <- file.path(CODE_DIR, "00_globals_utils.R")
stopifnot(file.exists(globals_path))
source(globals_path, local = FALSE)

# After sourcing globals, these should exist:
stopifnot(exists("dirs"), exists("cfg"), exists("par_base"))

# Hard requirement (recommended) for analytic Hopf
if (isTRUE(REQUIRE_ANALYTIC_HOPF)) {
  if (!exists("compute_at_rF", mode = "function")) {
    stop("compute_at_rF not found. Fix: ensure 00b_wealthIneq_analytic_core.R exists and is sourced by 00_globals_utils.R.")
  }
} else {
  if (!exists("compute_at_rF", mode = "function")) {
    message("[INFO] compute_at_rF not found: Stage 3 may run in passthrough/proxy mode (if your Stage 3 script supports it).")
  }
}

# Canonical stage dirs
dir_stage1 <- dirs$stage1_backbone
dir_stage2 <- dirs$stage2_finance
dir_stage3 <- dirs$stage3_stability
dir_stage4 <- dirs$stage4_scoring
dir_stage5 <- dirs$stage5_simulation

# ----------------------------
# Helper: source with logging + timing
# ----------------------------
source_stage <- function(fname) {
  f <- file.path(CODE_DIR, fname)
  stopifnot(file.exists(f))
  cat("\n----------------------------------------\n")
  cat("SOURCING:", fname, "\n")
  cat("----------------------------------------\n")
  t0 <- proc.time()[3]
  source(f, local = FALSE)
  t1 <- proc.time()[3]
  cat(sprintf("[DONE] %s in %.2f sec\n", fname, (t1 - t0)))
  invisible(TRUE)
}

assert_exists <- function(paths, label = "") {
  ok <- file.exists(paths)
  if (!all(ok)) {
    cat("\n[MISSING OUTPUTS]", label, "\n")
    print(tibble(path = paths, exists = ok))
    stop("Some expected outputs are missing. Stage failed or paths mismatch.")
  }
  invisible(TRUE)
}

warn_if_missing <- function(paths, label = "") {
  ok <- file.exists(paths)
  if (!all(ok)) {
    cat("\n[WARN] Some optional outputs missing:", label, "\n")
    print(tibble(path = paths, exists = ok))
  }
  invisible(TRUE)
}

# ----------------------------
# Stage filenames
# ----------------------------
STAGE_FILES <- list(
  s1 = "01_stage1_backbone.R",
  s2 = "02_stage2_finance.R",
  s3 = "03_stage3_stability_hopf.R",
  s4 = "04_stage4_scoring_shortlist.R",
  s5 = "05_stage5_pathwayA_light.R"
)

# ----------------------------
# Key outputs (authoritative names)
# ----------------------------
OUT1 <- file.path(dir_stage1, "stage1_backbone.csv")
OUT2 <- file.path(dir_stage2, "stage2_finance_scan.csv")
OUT3 <- file.path(dir_stage3, c("stage3_candidates.csv","hopf_roots.csv"))
OUT4 <- file.path(dir_stage4, "stage4_scored_candidates.csv")

# Stage 5 expected outputs
OUT5 <- c(
  file.path(dirs$stage5A, "stage5A_candidates.csv"),
  file.path(dirs$stage5B, "stage5B_flip_metrics.csv"),
  file.path(dirs$stage5B, "stage5B2_expand_local_sweep.csv"),
  file.path(dirs$stage5C, "stage5C_robustness.csv"),
  file.path(dirs$stage5D, "stage5D_ranked.csv"),
  file.path(dirs$stage5E, "stage5E_pivot_report.csv"),
  file.path(dirs$stage5E, "stage5E_pivot_recommendation.txt"),
  file.path(dirs$stage5E, "stage5E_failure_mode_shares.csv")
)

# ============================================================
# STAGE 1
# ============================================================
if (isTRUE(RUN_STAGE1)) {
  if (isTRUE(SKIP_IF_EXISTS) && file.exists(OUT1)) {
    cat("\n[SKIP] Stage 1 outputs already exist:", OUT1, "\n")
  } else {
    source_stage(STAGE_FILES$s1)
    assert_exists(OUT1, "Stage 1")
    warn_if_missing(file.path(dir_stage1, "failure_reasons.csv"), "Stage 1")
  }
} else {
  cat("\n[SKIP] RUN_STAGE1=FALSE\n")
}

# ============================================================
# STAGE 2
# ============================================================
if (isTRUE(RUN_STAGE2)) {
  if (isTRUE(SKIP_IF_EXISTS) && file.exists(OUT2)) {
    cat("\n[SKIP] Stage 2 outputs already exist:", OUT2, "\n")
  } else {
    source_stage(STAGE_FILES$s2)
    assert_exists(OUT2, "Stage 2")
    warn_if_missing(file.path(dir_stage2, "failure_reasons.csv"), "Stage 2")
  }
} else {
  cat("\n[SKIP] RUN_STAGE2=FALSE\n")
}

# ============================================================
# STAGE 3 PREFLIGHT (cheap audit before expensive run)
# ============================================================
preflight_stage2 <- function(path_out2) {
  if (!file.exists(path_out2)) stop("Stage 2 file missing: ", path_out2)
  df2 <- readr::read_csv(path_out2, show_col_types = FALSE)
  
  # minimum columns Stage 3 will need
  need_cols <- c("sigma","g_n","i","delta","kappa_max","r_star","d_star","omega_star","rF","psi","phi2")
  miss <- setdiff(need_cols, names(df2))
  if (length(miss) > 0) stop("Stage 2 output missing required columns: ", paste(miss, collapse = ", "))
  
  # survivors: prefer econ_ok TRUE if present
  if ("econ_ok" %in% names(df2)) {
    n_ok <- sum(isTRUE(df2$econ_ok), na.rm = TRUE)
    cat("\n[PREFLIGHT] Stage 2 econ_ok survivors:", n_ok, "\n")
    if (n_ok == 0) stop("No econ_ok survivors in Stage 2. Do NOT run Stage 3. Fix Stage 1–2 grid/filters first.")
  } else {
    cat("\n[PREFLIGHT] Stage 2 has no econ_ok column. Proceeding, but you’re flying without instruments.\n")
  }
  
  # phi1 lock check (high value sanity)
  if (any(c("phi1_fixed","phi1") %in% names(df2))) {
    col <- if ("phi1_fixed" %in% names(df2)) "phi1_fixed" else "phi1"
    n_phi <- sum(is.finite(df2[[col]]))
    cat("[PREFLIGHT] Stage 2 finite ", col, ": ", n_phi, "\n", sep = "")
    if (n_phi == 0 && isTRUE(REQUIRE_ANALYTIC_HOPF)) {
      stop("No finite phi1 values in Stage 2. But analytic Hopf requires phi1 locked. Fix Stage 2 calibration.")
    }
  }
  
  invisible(TRUE)
}

# ============================================================
# STAGE 3
# ============================================================
if (isTRUE(RUN_STAGE3)) {
  preflight_stage2(OUT2)
  
  if (isTRUE(SKIP_IF_EXISTS) && all(file.exists(OUT3))) {
    cat("\n[SKIP] Stage 3 outputs already exist:\n  ", paste(OUT3, collapse = "\n  "), "\n", sep = "")
  } else {
    source_stage(STAGE_FILES$s3)
    assert_exists(OUT3, "Stage 3")
    warn_if_missing(file.path(dir_stage3, "failure_reasons.csv"), "Stage 3")
  }
} else {
  cat("\n[SKIP] RUN_STAGE3=FALSE\n")
}

# ============================================================
# STAGE 4
# ============================================================
if (isTRUE(RUN_STAGE4)) {
  if (isTRUE(SKIP_IF_EXISTS) && file.exists(OUT4)) {
    cat("\n[SKIP] Stage 4 outputs already exist:", OUT4, "\n")
  } else {
    source_stage(STAGE_FILES$s4)
    assert_exists(OUT4, "Stage 4")
    warn_if_missing(
      file.path(dir_stage4, c("manifest_stage4.csv","stage4_shortlist_top50_GATED.csv")),
      "Stage 4 (optional)"
    )
  }
} else {
  cat("\n[SKIP] RUN_STAGE4=FALSE\n")
}

# ============================================================
# STAGE 5
# ============================================================
if (isTRUE(RUN_STAGE5)) {
  source_stage(STAGE_FILES$s5)
  
  if (!exists("run_stage5_pathwayA_light")) {
    stop("Stage 5 script sourced, but run_stage5_pathwayA_light() not found.")
  }
  
  # Stage 5 requires a simulator hook
  if (!exists("simulate_model", mode = "function")) {
    stop("simulate_model() not found. It should be provided by 00b_wealthIneq_analytic_core.R via 00_globals_utils.R.")
  }
  
  cat("\n----------------------------------------\n")
  cat("RUNNING: Stage 5 (run_stage5_pathwayA_light)\n")
  cat("----------------------------------------\n")
  out5 <- run_stage5_pathwayA_light(cfg = cfg, dirs = dirs)
  
  assert_exists(OUT5, "Stage 5")
  
} else {
  cat("\n[SKIP] RUN_STAGE5=FALSE\n")
}

# ============================================================
# MANIFEST: plots produced
# ============================================================
PLOT_EXT <- c("png","pdf","jpg","jpeg","svg")
OUT_ROOT <- dirs$out_root

match_dir <- function(path, dir) {
  if (is.null(dir) || length(dir) == 0) return(FALSE)
  if (!is.character(dir) || !nzchar(dir[1])) return(FALSE)
  stringr::str_detect(tolower(path), stringr::fixed(tolower(dir[1])))
}

all_files <- list.files(OUT_ROOT, recursive = TRUE, full.names = TRUE)
manifest_plots <- tibble(
  path = all_files,
  ext  = tolower(tools::file_ext(all_files))
) %>%
  filter(ext %in% PLOT_EXT) %>%
  mutate(
    bytes = file.size(path),
    mtime = as.POSIXct(file.info(path)$mtime, tz = "")
  ) %>%
  arrange(path) %>%
  mutate(
    stage = case_when(
      match_dir(path, dirs$stage1_backbone) ~ "stage1_backbone",
      match_dir(path, dirs$stage2_finance)  ~ "stage2_finance",
      match_dir(path, dirs$stage3_stability)~ "stage3_stability",
      match_dir(path, dirs$stage4_scoring)  ~ "stage4_scoring",
      match_dir(path, dirs$stage5A) ~ "stage5A",
      match_dir(path, dirs$stage5B) ~ "stage5B",
      match_dir(path, dirs$stage5C) ~ "stage5C",
      match_dir(path, dirs$stage5D) ~ "stage5D",
      match_dir(path, dirs$stage5E) ~ "stage5E",
      match_dir(path, dirs$stage5_simulation) ~ "stage5_simulation",
      TRUE ~ "other"
    )
  )

dir.create(dirs$analysis, showWarnings = FALSE, recursive = TRUE)
manifest_path <- file.path(dirs$analysis, "manifest_plots.csv")
readr::write_csv(manifest_plots, manifest_path)

cat("\n====================\nPLOT MANIFEST SAVED\n====================\n")
cat("CSV:", manifest_path, "\n")

cat("\n====================\nPLOT SUMMARY BY STAGE/EXT\n====================\n")
print(manifest_plots %>% count(stage, ext, sort = TRUE))

cat("\n[RUN MANAGER COMPLETE]\n")
