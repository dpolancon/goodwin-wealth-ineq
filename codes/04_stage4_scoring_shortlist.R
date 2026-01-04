# ============================================================
# 04_stage4_scoring_shortlist.R  (CONSOLIDATED + WORKFLOW-SAFE)
# Stage 4: gating + scoring + shortlist
#
# - Robustly finds Stage 3 input
# - Normalizes schema: e_star created from e_implied if needed
# - Normalizes ok/reason fields
# - Optional: auto-relax hopf_required if Stage 3 has no Hopf roots
#
# Outputs -> dirs$stage4_scoring:
#   - stage4_scored_candidates.csv
#   - stage4_ranked.csv
#   - stage4_shortlist_top50_GATED.csv
#   - stage4_shortlist.csv
#   - stage4_gate_counts.csv
#   - manifest_stage4.csv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

cat("\n====================\n")
cat("Stage 4: Scoring + Shortlist (standalone-safe)\n")
cat("====================\n")

# ----------------------------
# Load globals + dirs if needed
# ----------------------------
if (!exists("dirs")) {
  if (file.exists("00_globals_utils.R")) {
    source("00_globals_utils.R", local = FALSE)
  } else {
    stop("Stage 4: dirs not found and 00_globals_utils.R not in working directory.")
  }
}

# ----------------------------
# Find Stage 3 input robustly
# ----------------------------
stage3_candidates <- c(
  file.path(dirs$stage3_stability, "stage3_candidates.csv"),
  file.path(dirs$out_root, "stage3_stability", "stage3_candidates.csv"),
  file.path(dirs$out_root, "stage3", "stage3_candidates.csv"),
  "outputs/wealth_goodwin/stage3_stability/stage3_candidates.csv",
  "outputs/wealth_goodwin/stage3/stage3_candidates.csv"
)

in_file <- stage3_candidates[file.exists(stage3_candidates)][1]
if (is.na(in_file) || !nzchar(in_file)) {
  stop("Stage 4: cannot find Stage 3 input. Looked for:\n  - ",
       paste(stage3_candidates, collapse = "\n  - "))
}
cat("Stage 4: using input:\n  ", in_file, "\n", sep = "")

# ----------------------------
# Outputs
# ----------------------------
out_dir <- dirs$stage4_scoring
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_scored   <- file.path(out_dir, "stage4_scored_candidates.csv")
out_ranked   <- file.path(out_dir, "stage4_ranked.csv")
out_short50  <- file.path(out_dir, "stage4_shortlist_top50_GATED.csv")
out_short    <- file.path(out_dir, "stage4_shortlist.csv")
out_gates    <- file.path(out_dir, "stage4_gate_counts.csv")
out_manifest <- file.path(out_dir, "manifest_stage4.csv")

df <- readr::read_csv(in_file, show_col_types = FALSE)

# ----------------------------
# Bands (EDIT HERE)
# ----------------------------
B <- list(
  e_min = 0.20, e_max = 0.98,
  w_min = 0.40, w_max = 0.95,
  d_min = 0.00, d_max = 3.00,
  
  hopf_required = TRUE,
  auto_relax_hopf_if_none = TRUE   # lifesaver when Stage 3 has no roots
)

# ----------------------------
# Normalize schema
# ----------------------------
if (!("e_star" %in% names(df))) {
  if ("e_implied" %in% names(df)) {
    df <- df %>% mutate(e_star = .data$e_implied)
  } else if (exists("cfg", inherits = TRUE) && !is.null(get("cfg", inherits = TRUE)$targets$e_target)) {
    df <- df %>% mutate(e_star = get("cfg", inherits = TRUE)$targets$e_target)
  } else {
    df <- df %>% mutate(e_star = NA_real_)
  }
}

if (!("ok" %in% names(df))) {
  df$ok <- if ("sim_ok" %in% names(df)) df$sim_ok else TRUE
}

# prefer reason, else reason.y/x, else sim_reason, else "ok"
if (!("reason" %in% names(df))) {
  if ("reason.y" %in% names(df))      df$reason <- df$reason.y
  else if ("reason.x" %in% names(df)) df$reason <- df$reason.x
  else if ("sim_reason" %in% names(df)) df$reason <- df$sim_reason
  else df$reason <- "ok"
} else {
  # fill blank/NA reason
  df <- df %>% mutate(
    reason = ifelse(is.na(.data$reason) | .data$reason == "",
                    dplyr::coalesce(.data$reason.y, .data$reason.x, .data$sim_reason, "ok"),
                    .data$reason)
  )
}

# Ensure stability fields exist (may be NA)
if (!("stable" %in% names(df)))      df$stable <- NA
if (!("RH_ok" %in% names(df)))       df$RH_ok <- NA
if (!("has_complex" %in% names(df))) df$has_complex <- NA

# Ensure hopf fields exist
for (nm in c("has_hopf","n_roots","rF_root_med")) {
  if (!(nm %in% names(df))) df[[nm]] <- NA
}

# Ensure core state vars exist
if (!("omega_star" %in% names(df))) df$omega_star <- NA_real_
if (!("d_star" %in% names(df)))     df$d_star <- NA_real_

# Auto-relax Hopf requirement if none exist (optional)
hopf_required <- isTRUE(B$hopf_required)
if (hopf_required && isTRUE(B$auto_relax_hopf_if_none)) {
  has_any_hopf <- any(df$has_hopf %in% TRUE, na.rm = TRUE)
  if (!has_any_hopf) {
    hopf_required <- FALSE
    cat("[Stage 4] No Hopf roots in Stage 3 -> auto_relax_hopf_if_none = TRUE, disabling hopf_required.\n")
  }
}

# ----------------------------
# Gates
# ----------------------------
df2 <- df %>%
  mutate(
    ok_flag = dplyr::coalesce(as.logical(.data$ok), TRUE),
    
    gate_interior =
      is.finite(.data$e_star) & .data$e_star >= B$e_min & .data$e_star <= B$e_max &
      is.finite(.data$omega_star) & .data$omega_star >= B$w_min & .data$omega_star <= B$w_max &
      is.finite(.data$d_star) & .data$d_star >= B$d_min & .data$d_star <= B$d_max,
    
    gate_hopf = if (hopf_required) {
      (.data$has_hopf %in% TRUE) &
        is.finite(.data$n_roots) & .data$n_roots > 0 &
        is.finite(.data$rF_root_med)
    } else TRUE,
    
    gate_ok = .data$ok_flag & .data$gate_interior & .data$gate_hopf,
    
    gate_reason = case_when(
      !.data$ok_flag ~ paste0("analytic_fail:", .data$reason),
      !.data$gate_interior ~ "reject_interior_bounds",
      hopf_required & !(.data$has_hopf %in% TRUE) ~ "reject_no_hopf_root",
      TRUE ~ "pass"
    )
  )

gate_counts <- df2 %>%
  count(.data$gate_reason, sort = TRUE) %>%
  mutate(share = .data$n / sum(.data$n))
readr::write_csv(gate_counts, out_gates)

# ----------------------------
# Scoring (only for gate_ok)
# ----------------------------
mid <- list(e = (B$e_min + B$e_max)/2, w = (B$w_min + B$w_max)/2, d = 0.75)

score_one <- function(e, w, d, stable, RH_ok, has_complex) {
  if (!is.finite(e) || !is.finite(w) || !is.finite(d)) return(NA_real_)
  pe <- abs(e - mid$e)
  pw <- abs(w - mid$w)
  pd <- abs(d - mid$d) / 2
  
  bonus <- 0
  if (isTRUE(RH_ok)) bonus <- bonus + 2
  if (isTRUE(stable)) bonus <- bonus + 1
  if (isTRUE(has_complex)) bonus <- bonus + 0.5
  
  bonus - (2*pe + 2*pw + pd)
}

df3 <- df2 %>%
  mutate(
    score_stage4 = ifelse(
      .data$gate_ok,
      mapply(score_one, .data$e_star, .data$omega_star, .data$d_star,
             .data$stable, .data$RH_ok, .data$has_complex),
      NA_real_
    )
  ) %>%
  arrange(desc(.data$gate_ok), desc(.data$score_stage4)) %>%
  mutate(rank_stage4 = row_number())

TOP_N <- 50
short <- df3 %>% filter(.data$gate_ok) %>% slice_head(n = TOP_N)

# ----------------------------
# Write outputs
# ----------------------------
readr::write_csv(df3, out_scored)
readr::write_csv(df3, out_ranked)
readr::write_csv(short, out_short50)
readr::write_csv(short, out_short)

manifest <- tibble(
  file = c(out_scored, out_ranked, out_short50, out_short, out_gates),
  exists = file.exists(c(out_scored, out_ranked, out_short50, out_short, out_gates))
)
readr::write_csv(manifest, out_manifest)

cat("\nStage 4 outputs written:\n")
cat("  - ", out_scored, "\n", sep = "")
cat("  - ", out_ranked, "\n", sep = "")
cat("  - ", out_short50, "\n", sep = "")
cat("  - ", out_short, "\n", sep = "")
cat("  - ", out_gates, "\n", sep = "")
cat("Stage 4 done.\n")
