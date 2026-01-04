# ============================================================
# 05_stage5_pathwayA_light_REDO_CONSOLIDATED.R
# Stage 5: Hopf neighborhood sims (paper-first, workflow-safe)
#
# Fixes:
#  - NO group_modify() (avoid dplyr grouping gotchas)
#  - bad_x0 gets 2 fallbacks: unperturbed x0, then x0=NULL
#  - cycle_metrics() guards against all-NA vectors (no max/min warnings)
#  - bad_x0 counted as numeric taint
#
# Outputs -> dirs$stage5_simulation (or out_root fallback):
#   5A/stage5A_candidates.csv
#   5B/stage5B_flip_metrics.csv
#   5B/stage5B2_expand_local_sweep.csv
#   5C/stage5C_robustness.csv
#   5D/stage5D_ranked.csv
#   5E/stage5E_pivot_report.csv
#   5E/stage5E_failure_mode_shares.csv
# ============================================================

stage5_run <- function() {
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(purrr)
    library(tibble)
  })
  
  cat("\n====================\n")
  cat("Stage 5: Pathway A (REDO) — Hopf neighborhood sims\n")
  cat("====================\n")
  
  `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
  
  empty_with_cols <- function(cols) {
    out <- tibble::tibble()
    for (nm in names(cols)) out[[nm]] <- cols[[nm]][0]
    out
  }
  
  # ----------------------------
  # Load globals + dirs if needed
  # ----------------------------
  if (!exists("dirs")) {
    if (file.exists("00_globals_utils.R")) {
      source("00_globals_utils.R", local = FALSE)
    } else {
      stop("Stage 5: dirs not found and 00_globals_utils.R not in working directory.")
    }
  }
  
  if (!exists("simulate_model", mode = "function")) {
    stop("Stage 5: simulate_model() missing. Source 00b_wealthIneq_analytic_core.R first.")
  }
  
  # ----------------------------
  # Find Stage 4 shortlist robustly
  # ----------------------------
  short_candidates <- c(
    file.path(dirs$stage4_scoring, "stage4_shortlist.csv"),
    file.path(dirs$stage4_scoring, "stage4_shortlist_top50_GATED.csv"),
    file.path(dirs$out_root, "stage4_scoring", "stage4_shortlist.csv"),
    file.path(dirs$out_root, "stage4_scoring", "stage4_shortlist_top50_GATED.csv"),
    "outputs/wealth_goodwin/stage4_scoring/stage4_shortlist.csv",
    "outputs/wealth_goodwin/stage4_scoring/stage4_shortlist_top50_GATED.csv",
    "outputs/wealth_goodwin/stage4/stage4_shortlist.csv"
  )
  
  in_short <- short_candidates[file.exists(short_candidates)][1]
  if (is.na(in_short) || !nzchar(in_short)) {
    stop("Stage 5: cannot find Stage 4 shortlist. Looked for:\n  - ",
         paste(short_candidates, collapse = "\n  - "))
  }
  cat("Stage 5: using shortlist:\n  ", in_short, "\n", sep = "")
  
  # ----------------------------
  # Output dirs
  # ----------------------------
  base_dir <- dirs$stage5_simulation %||% file.path(dirs$out_root, "stage5_simulation")
  dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  
  dirA <- file.path(base_dir, "5A"); dir.create(dirA, recursive = TRUE, showWarnings = FALSE)
  dirB <- file.path(base_dir, "5B"); dir.create(dirB, recursive = TRUE, showWarnings = FALSE)
  dirC <- file.path(base_dir, "5C"); dir.create(dirC, recursive = TRUE, showWarnings = FALSE)
  dirD <- file.path(base_dir, "5D"); dir.create(dirD, recursive = TRUE, showWarnings = FALSE)
  dirE <- file.path(base_dir, "5E"); dir.create(dirE, recursive = TRUE, showWarnings = FALSE)
  
  outA  <- file.path(dirA, "stage5A_candidates.csv")
  outB  <- file.path(dirB, "stage5B_flip_metrics.csv")
  outB2 <- file.path(dirB, "stage5B2_expand_local_sweep.csv")
  outC  <- file.path(dirC, "stage5C_robustness.csv")
  outD  <- file.path(dirD, "stage5D_ranked.csv")
  outE  <- file.path(dirE, "stage5E_pivot_report.csv")
  outE2 <- file.path(dirE, "stage5E_failure_mode_shares.csv")
  
  # ----------------------------
  # Tunables (EDIT HERE)
  # ----------------------------
  t_end <- 600
  dt    <- 0.05
  
  hopf_eps <- 0.002
  perturb  <- c(e = 0.002, omega = -0.002, d = 0.02)
  
  bounds <- list(
    e     = c(0.0, 1.2),
    omega = c(0.0, 0.999),
    d     = c(0.0, 10.0)
  )
  
  var_min <- 1e-5
  amp_min <- 0.02
  burn_in <- 200
  
  ROBUST_N <- 12
  set.seed(123)
  
  # ----------------------------
  # Helpers
  # ----------------------------
  as_rowlist <- function(df_row) as.list(df_row)
  
  get_stop_reason <- function(sim) {
    sr <- attr(sim, "stop_reason")
    if (!is.null(sr) && length(sr) > 0) return(as.character(sr)[1])
    "ok"
  }
  
  # Treat bad_x0 as taint (it is: you have no trajectory, just a refusal)
  is_tainted <- function(stop_reason) {
    if (is.na(stop_reason) || !nzchar(stop_reason)) return(FALSE)
    grepl("bad_x0|event_|rk|nan|inf|step|overflow|underflow|singular|ode_error|simulate_error",
          stop_reason, ignore.case = TRUE)
  }
  
  check_bounded <- function(sim) {
    if (!all(c("e","omega","d") %in% names(sim))) return(list(ok=FALSE, reason="missing_state_cols"))
    e_ok <- all(sim$e >= bounds$e[1] & sim$e <= bounds$e[2], na.rm = TRUE)
    w_ok <- all(sim$omega >= bounds$omega[1] & sim$omega <= bounds$omega[2], na.rm = TRUE)
    d_ok <- all(sim$d >= bounds$d[1] & sim$d <= bounds$d[2], na.rm = TRUE)
    ok <- e_ok & w_ok & d_ok
    reason <- if (ok) "ok" else paste(
      c(if(!e_ok) "e_out_of_bounds" else NULL,
        if(!w_ok) "omega_out_of_bounds" else NULL,
        if(!d_ok) "d_out_of_bounds" else NULL),
      collapse=";"
    )
    list(ok = ok, reason = reason)
  }
  
  cycle_metrics <- function(sim) {
    if (!("time" %in% names(sim))) {
      return(tibble(has_cycle=FALSE, reason="missing_time", amp=NA_real_, period=NA_real_, omega_mean=NA_real_))
    }
    sim2 <- as_tibble(sim) %>% filter(.data$time >= burn_in)
    if (!("omega" %in% names(sim2))) {
      return(tibble(has_cycle=FALSE, reason="missing_omega", amp=NA_real_, period=NA_real_, omega_mean=NA_real_))
    }
    
    x <- sim2$omega
    t <- sim2$time
    
    # critical guard: no finite points -> no warnings -> no garbage
    if (length(x) < 5 || all(!is.finite(x))) {
      return(tibble(has_cycle=FALSE, reason="all_na_or_too_short", amp=NA_real_, period=NA_real_, omega_mean=NA_real_))
    }
    
    v <- var(x, na.rm = TRUE)
    a <- (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    m <- mean(x, na.rm = TRUE)
    
    if (!is.finite(v) || !is.finite(a)) {
      return(tibble(has_cycle=FALSE, reason="na_metrics", amp=NA_real_, period=NA_real_, omega_mean=m))
    }
    if (v < var_min || a < amp_min) {
      return(tibble(has_cycle=FALSE, reason="low_var_or_amp", amp=a, period=NA_real_, omega_mean=m))
    }
    
    dx <- diff(x)
    sgn <- sign(dx)
    turns <- which(diff(sgn) == -2)  # local maxima heuristic
    if (length(turns) < 3) {
      return(tibble(has_cycle=TRUE, reason="cycle_detected_weak_period", amp=a, period=NA_real_, omega_mean=m))
    }
    
    tt <- t[turns]
    per <- suppressWarnings(median(diff(tt), na.rm = TRUE))
    
    tibble(has_cycle=TRUE, reason="cycle_detected", amp=a, period=per, omega_mean=m)
  }
  
  # main simulator wrapper with fallbacks for bad_x0
  safe_simulate <- function(row, rF_sim, t_end, dt, x0_try, x0_fallback, allow_null = TRUE) {
    
    run_once <- function(x0) {
      tryCatch(
        simulate_model(row, rF_sim = rF_sim, t_end = t_end, dt = dt, x0 = x0),
        error = function(e) {
          sim <- tibble(time = NA_real_, e = NA_real_, omega = NA_real_, d = NA_real_)
          attr(sim, "stop_reason") <- paste0("simulate_error:", conditionMessage(e))
          sim
        }
      )
    }
    
    sim1 <- run_once(x0_try)
    sr1  <- get_stop_reason(sim1)
    if (!identical(sr1, "bad_x0")) return(sim1)
    
    sim2 <- run_once(x0_fallback)
    sr2  <- get_stop_reason(sim2)
    if (!identical(sr2, "bad_x0")) {
      attr(sim2, "stop_reason") <- paste0(sr2, "|fallback_unperturbed")
      return(sim2)
    }
    
    if (isTRUE(allow_null)) {
      sim3 <- run_once(NULL)
      sr3  <- get_stop_reason(sim3)
      attr(sim3, "stop_reason") <- paste0(sr3, "|fallback_null")
      return(sim3)
    }
    
    attr(sim2, "stop_reason") <- paste0(sr2, "|still_bad_x0")
    sim2
  }
  
  # ----------------------------
  # Load shortlist and normalize schema
  # ----------------------------
  short <- readr::read_csv(in_short, show_col_types = FALSE)
  
  if (!("candidate_id" %in% names(short))) short <- short %>% mutate(candidate_id = row_number())
  
  if (!("e_star" %in% names(short))) {
    if ("e_implied" %in% names(short)) short <- short %>% mutate(e_star = .data$e_implied)
    else short <- short %>% mutate(e_star = NA_real_)
  }
  if (!("omega_star" %in% names(short))) short <- short %>% mutate(omega_star = NA_real_)
  if (!("d_star" %in% names(short)))     short <- short %>% mutate(d_star = NA_real_)
  
  for (nm in c("has_hopf","n_roots","rF_root_med","rF")) {
    if (!(nm %in% names(short))) short[[nm]] <- NA
  }
  
  # ----------------------------
  # If shortlist empty: write schema-safe outputs and exit
  # ----------------------------
  if (nrow(short) == 0) {
    warning("Stage 5: shortlist is empty -> writing empty outputs (schema-safe) and exiting.")
    
    stage5A_empty <- empty_with_cols(list(
      candidate_id = integer(),
      rF_anchor = numeric(),
      anchor_source = character()
    ))
    
    flip_empty <- empty_with_cols(list(
      candidate_id = integer(),
      sim_ok = logical(),
      sim_reason = character(),
      rF_below = numeric(),
      rF_above = numeric(),
      stop_reason_below = character(),
      stop_reason_above = character(),
      tainted_below = logical(),
      tainted_above = logical(),
      bounded_below = logical(),
      bounded_reason_below = character(),
      bounded_above = logical(),
      bounded_reason_above = character(),
      cycle_below = logical(),
      cycle_reason_below = character(),
      cycle_above = logical(),
      cycle_reason_above = character(),
      amp_below = numeric(),
      amp_above = numeric(),
      period_below = numeric(),
      period_above = numeric(),
      omega_mean_below = numeric(),
      omega_mean_above = numeric(),
      flip_ok = logical(),
      class_5B = character()
    ))
    
    expand_empty <- empty_with_cols(list(
      candidate_id = integer(),
      rF_center = numeric(),
      rF_sim = numeric(),
      bounded = logical(),
      has_cycle = logical(),
      tainted = logical(),
      stop_reason = character(),
      reason = character()
    ))
    
    robust_empty <- empty_with_cols(list(
      candidate_id = integer(),
      survival_rate = numeric(),
      bounded_rate = numeric(),
      cycle_rate = numeric(),
      amp_med = numeric(),
      period_med = numeric(),
      taint_share = numeric()
    ))
    
    ranked_empty <- empty_with_cols(list(
      candidate_id = integer(),
      score_stage5 = numeric(),
      class_5D = character()
    ))
    
    pivot_empty <- empty_with_cols(list(
      n_total_sim_ok = integer(),
      n_flip_ok = integer(),
      flip_rate = numeric(),
      robust_survival_med = numeric(),
      share_unbounded = numeric(),
      share_tainted = numeric(),
      recommended_mode = character()
    ))
    
    failshares_empty <- empty_with_cols(list(
      class_5D = character(),
      n = integer(),
      share = numeric()
    ))
    
    readr::write_csv(stage5A_empty, outA)
    readr::write_csv(flip_empty, outB)
    readr::write_csv(expand_empty, outB2)
    readr::write_csv(robust_empty, outC)
    readr::write_csv(ranked_empty, outD)
    readr::write_csv(pivot_empty, outE)
    readr::write_csv(failshares_empty, outE2)
    
    message("Stage 5: done (empty shortlist). Outputs written to: ", base_dir)
    return(invisible(NULL))
  }
  
  # ----------------------------
  # Stage 5A candidates (anchor)
  # ----------------------------
  stage5A <- short %>%
    mutate(
      rF_anchor = ifelse(is.finite(.data$rF_root_med), .data$rF_root_med, .data$rF),
      anchor_source = ifelse(is.finite(.data$rF_root_med), "hopf_root_med", "baseline_rF")
    )
  
  readr::write_csv(stage5A, outA)
  
  # ----------------------------
  # Stage 5B flip test
  # ----------------------------
  row_lists <- split(stage5A, seq_len(nrow(stage5A)))
  
  flip_res <- purrr::map_dfr(row_lists, function(df_row) {
    row <- as_rowlist(df_row)
    
    has_hopf_ok <- (row$has_hopf %in% TRUE) &&
      is.finite(row$rF_root_med) &&
      is.finite(row$n_roots) && row$n_roots > 0
    
    if (!has_hopf_ok) {
      return(tibble(
        candidate_id = row$candidate_id,
        sim_ok = TRUE, sim_reason = "skip_no_hopf",
        rF_below = NA_real_, rF_above = NA_real_,
        stop_reason_below = NA_character_, stop_reason_above = NA_character_,
        tainted_below = NA, tainted_above = NA,
        bounded_below = NA, bounded_reason_below = NA_character_,
        bounded_above = NA, bounded_reason_above = NA_character_,
        cycle_below = NA, cycle_reason_below = NA_character_,
        cycle_above = NA, cycle_reason_above = NA_character_,
        amp_below = NA_real_, amp_above = NA_real_,
        period_below = NA_real_, period_above = NA_real_,
        omega_mean_below = NA_real_, omega_mean_above = NA_real_,
        flip_ok = FALSE,
        class_5B = "SKIP_NO_HOPF"
      ))
    }
    
    rF_below <- row$rF_root_med - hopf_eps
    rF_above <- row$rF_root_med + hopf_eps
    
    x0_try <- c(
      e     = as.numeric(row$e_star)     + perturb["e"],
      omega = as.numeric(row$omega_star) + perturb["omega"],
      d     = as.numeric(row$d_star)     + perturb["d"]
    )
    
    x0_fallback <- c(
      e     = as.numeric(row$e_star),
      omega = as.numeric(row$omega_star),
      d     = as.numeric(row$d_star)
    )
    
    sim_below <- safe_simulate(row, rF_below, t_end, dt, x0_try, x0_fallback, allow_null = TRUE)
    sim_above <- safe_simulate(row, rF_above, t_end, dt, x0_try, x0_fallback, allow_null = TRUE)
    
    sr_below <- get_stop_reason(sim_below)
    sr_above <- get_stop_reason(sim_above)
    
    ta_below <- is_tainted(sr_below)
    ta_above <- is_tainted(sr_above)
    
    b_below <- check_bounded(sim_below)
    b_above <- check_bounded(sim_above)
    
    cy_below <- cycle_metrics(sim_below)
    cy_above <- cycle_metrics(sim_above)
    
    flip_ok <- isTRUE(b_below$ok) && isTRUE(b_above$ok) &&
      !ta_below && !ta_above &&
      xor(isTRUE(cy_below$has_cycle), isTRUE(cy_above$has_cycle))
    
    class <- dplyr::case_when(
      ta_below || ta_above ~ "FAIL_NUMERIC_TAINT",
      !b_below$ok || !b_above$ok ~ "FAIL_UNBOUNDED",
      !xor(isTRUE(cy_below$has_cycle), isTRUE(cy_above$has_cycle)) ~ "NO_FLIP",
      TRUE ~ "FLIP_OK"
    )
    
    tibble(
      candidate_id = row$candidate_id,
      sim_ok = TRUE, sim_reason = "ok",
      rF_below = rF_below, rF_above = rF_above,
      stop_reason_below = sr_below,
      stop_reason_above = sr_above,
      tainted_below = ta_below,
      tainted_above = ta_above,
      bounded_below = b_below$ok, bounded_reason_below = b_below$reason,
      bounded_above = b_above$ok, bounded_reason_above = b_above$reason,
      cycle_below = cy_below$has_cycle, cycle_reason_below = cy_below$reason,
      cycle_above = cy_above$has_cycle, cycle_reason_above = cy_above$reason,
      amp_below = cy_below$amp, amp_above = cy_above$amp,
      period_below = cy_below$period, period_above = cy_above$period,
      omega_mean_below = cy_below$omega_mean, omega_mean_above = cy_above$omega_mean,
      flip_ok = flip_ok,
      class_5B = class
    )
  })
  
  if (!("flip_ok" %in% names(flip_res))) flip_res$flip_ok <- logical(nrow(flip_res))
  readr::write_csv(flip_res, outB)
  
  # ----------------------------
  # Stage 5B2: expand local sweep ONLY when flip_ok
  # ----------------------------
  if (nrow(flip_res) == 0 || !any(flip_res$flip_ok %in% TRUE, na.rm = TRUE)) {
    
    expand_res <- empty_with_cols(list(
      candidate_id = integer(),
      rF_center = numeric(),
      rF_sim = numeric(),
      bounded = logical(),
      has_cycle = logical(),
      tainted = logical(),
      stop_reason = character(),
      reason = character()
    ))
    
  } else {
    
    flip_ids <- flip_res %>% filter(.data$flip_ok %in% TRUE) %>% pull(.data$candidate_id)
    
    expand_res <- purrr::map_dfr(flip_ids, function(cid) {
      row_tbl <- stage5A %>% filter(.data$candidate_id == cid) %>% slice(1)
      row <- as_rowlist(row_tbl)
      
      if (!is.finite(row$rF_root_med)) return(tibble())
      
      grid <- seq(row$rF_root_med - 5*hopf_eps, row$rF_root_med + 5*hopf_eps, by = hopf_eps)
      
      x0_try <- c(
        e     = as.numeric(row$e_star)     + perturb["e"],
        omega = as.numeric(row$omega_star) + perturb["omega"],
        d     = as.numeric(row$d_star)     + perturb["d"]
      )
      
      x0_fallback <- c(
        e     = as.numeric(row$e_star),
        omega = as.numeric(row$omega_star),
        d     = as.numeric(row$d_star)
      )
      
      purrr::map_dfr(grid, function(rF_sim) {
        sim <- safe_simulate(row, rF_sim, t_end, dt, x0_try, x0_fallback, allow_null = TRUE)
        sr  <- get_stop_reason(sim)
        ta  <- is_tainted(sr)
        b   <- check_bounded(sim)
        cy  <- cycle_metrics(sim)
        
        tibble(
          candidate_id = cid,
          rF_center = row$rF_root_med,
          rF_sim = rF_sim,
          bounded = b$ok,
          has_cycle = cy$has_cycle,
          tainted = ta,
          stop_reason = sr,
          reason = if (ta) "tainted" else if (!b$ok) "unbounded" else if (!cy$has_cycle) "no_cycle" else "cycle"
        )
      })
    })
  }
  
  readr::write_csv(expand_res, outB2)
  
  # ----------------------------
  # Stage 5C: robustness ring (no group_modify)
  # ----------------------------
  do_robust_ids <- flip_res %>%
    filter(.data$class_5B %in% c("FLIP_OK","NO_FLIP")) %>%
    pull(.data$candidate_id)
  
  robust <- purrr::map_dfr(stage5A$candidate_id, function(cid) {
    
    if (!(cid %in% do_robust_ids)) {
      return(tibble(
        candidate_id = cid,
        survival_rate = 0, bounded_rate = 0, cycle_rate = 0,
        amp_med = NA_real_, period_med = NA_real_,
        taint_share = 1
      ))
    }
    
    row_tbl <- stage5A %>% filter(.data$candidate_id == cid) %>% slice(1)
    row <- as_rowlist(row_tbl)
    
    set.seed(123 + cid)
    
    sims <- purrr::map_dfr(seq_len(ROBUST_N), function(j) {
      
      jit <- c(
        e     = as.numeric(row$e_star)     + rnorm(1, 0, perturb["e"]),
        omega = as.numeric(row$omega_star) + rnorm(1, 0, abs(perturb["omega"])),
        d     = as.numeric(row$d_star)     + rnorm(1, 0, perturb["d"])
      )
      
      x0_fallback <- c(
        e     = as.numeric(row$e_star),
        omega = as.numeric(row$omega_star),
        d     = as.numeric(row$d_star)
      )
      
      rF_sim <- as.numeric(row$rF_anchor) + rnorm(1, 0, hopf_eps)
      
      sim <- safe_simulate(row, rF_sim, t_end, dt, jit, x0_fallback, allow_null = TRUE)
      sr  <- get_stop_reason(sim)
      ta  <- is_tainted(sr)
      b   <- check_bounded(sim)
      cy  <- cycle_metrics(sim)
      
      tibble(
        ok = !ta,
        bounded = b$ok,
        has_cycle = cy$has_cycle,
        amp = cy$amp,
        period = cy$period,
        tainted = ta
      )
    })
    
    tibble(
      candidate_id = cid,
      survival_rate = mean(sims$ok, na.rm=TRUE),
      bounded_rate  = mean(sims$bounded, na.rm=TRUE),
      cycle_rate    = mean(sims$has_cycle, na.rm=TRUE),
      amp_med       = suppressWarnings(median(sims$amp, na.rm=TRUE)),
      period_med    = suppressWarnings(median(sims$period, na.rm=TRUE)),
      taint_share   = mean(sims$tainted, na.rm=TRUE)
    )
  })
  
  readr::write_csv(robust, outC)
  
  # ----------------------------
  # Stage 5D: merge + rank
  # ----------------------------
  ranked <- stage5A %>%
    left_join(flip_res, by="candidate_id") %>%
    left_join(robust, by="candidate_id") %>%
    mutate(
      class_5D = case_when(
        .data$class_5B == "SKIP_NO_HOPF" ~ "SKIP_NO_HOPF",
        .data$class_5B == "FAIL_NUMERIC_TAINT" ~ "FAIL_NUMERIC_TAINT",
        .data$class_5B == "FAIL_UNBOUNDED" ~ "FAIL_UNBOUNDED",
        .data$flip_ok %in% TRUE ~ "FLIP_OK",
        TRUE ~ "NO_FLIP"
      ),
      score_stage5 = case_when(
        .data$class_5D %in% c("FAIL_NUMERIC_TAINT","FAIL_UNBOUNDED","SKIP_NO_HOPF") ~ 0,
        TRUE ~ 10*survival_rate + 5*bounded_rate + 2*cycle_rate + ifelse(.data$flip_ok %in% TRUE, 2, 0)
      )
    ) %>%
    arrange(desc(.data$score_stage5))
  
  readr::write_csv(ranked, outD)
  
  # ----------------------------
  # Stage 5E: pivot report
  # ----------------------------
  pivot <- ranked %>%
    summarise(
      n_total_sim_ok = n(),
      n_flip_ok = sum(.data$flip_ok %in% TRUE, na.rm=TRUE),
      flip_rate = ifelse(n_total_sim_ok > 0, n_flip_ok / n_total_sim_ok, NA_real_),
      robust_survival_med = suppressWarnings(median(.data$survival_rate, na.rm=TRUE)),
      share_unbounded = mean(.data$class_5D == "FAIL_UNBOUNDED", na.rm=TRUE),
      share_tainted = mean(.data$class_5D == "FAIL_NUMERIC_TAINT", na.rm=TRUE)
    ) %>%
    mutate(
      recommended_mode = case_when(
        share_tainted > 0.50 ~ "repair_numeric_then_retry",
        is.na(flip_rate) || flip_rate < 0.05 ~ "expand_hopf_scan_in_stage3_or_relax_stage4_hopf_gate",
        TRUE ~ "proceed_pathwayA"
      )
    )
  
  fail_shares <- ranked %>%
    count(.data$class_5D, name="n") %>%
    mutate(share = .data$n / sum(.data$n))
  
  readr::write_csv(pivot, outE)
  readr::write_csv(fail_shares, outE2)
  
  message("Stage 5: wrote outputs to ", base_dir)
  message("Stage 5: pivot recommendation = ", pivot$recommended_mode[1])
  
  invisible(list(base_dir = base_dir, pivot = pivot))
}

# Run Stage 5
stage5_run()
