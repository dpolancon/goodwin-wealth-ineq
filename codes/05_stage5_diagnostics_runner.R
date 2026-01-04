# ============================================================
# 05_stage5_diagnostics_runner.R
# Stage 5 diagnostics: flip existence, robustness join, clamp prevalence
# Writes:
#   stage5E/diag_stage5_report.txt
#   stage5E/diag_flip_presence_by_candidate.csv
#   stage5E/stage5D_ranked_augmented.csv   (if join possible)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(purrr)
})

# ----------------------------
# 0) Paths (edit if your base differs)
# ----------------------------
base_dir <- "C:/ReposGitHub/goodwin_model/outputs/wealth_goodwin/stage5_simulation"

dirs <- list(
  stage5B  = file.path(base_dir, "5B"),
  stage5B2 = file.path(base_dir, "5B2"),
  stage5C  = file.path(base_dir, "5C"),
  stage5D  = file.path(base_dir, "5D"),
  stage5E  = file.path(base_dir, "5E")
)

dir.create(dirs$stage5E, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# 1) Helpers
# ----------------------------
exists_csv <- function(path) file.exists(path) && grepl("\\.csv$", path)

read_csv_if_exists <- function(path) {
  if (!file.exists(path)) return(NULL)
  readr::read_csv(path, show_col_types = FALSE)
}

has_col <- function(df, nm) !is.null(df) && (nm %in% names(df))

# Find best-guess column among options
pick_col <- function(df, options) {
  if (is.null(df)) return(NULL)
  hit <- options[options %in% names(df)]
  if (length(hit) == 0) return(NULL)
  hit[[1]]
}

# robust "candidate key" columns commonly used
candidate_keys <- function(df) {
  if (is.null(df)) return(character(0))
  keys <- c("try_id", "cand_id", "id", "candidate_id")
  keys[keys %in% names(df)]
}

# detect a numeric rF column
pick_rF_col <- function(df) {
  pick_col(df, c("rF", "rF_sim", "rFj", "r_F", "rF_value"))
}

# detect a "Hopf / stability sign metric" column
# preference: something like max real eigenvalue, trace, hopf_metric, etc.
pick_hopf_metric_col <- function(df) {
  pick_col(df, c(
    "eig_max_re", "max_re_lambda", "lambda_max_re", "maxRe", "max_real_eig",
    "trace", "trJ", "jac_trace",
    "hopf_metric", "stability_metric", "re_max"
  ))
}

# Given a vector of metric values, compute sign with tolerance
sign_tol <- function(x, tol = 1e-10) {
  ifelse(is.na(x), NA_integer_,
         ifelse(x > tol,  1L,
                ifelse(x < -tol, -1L, 0L)))
}

# detect sign changes in an ordered series (ignoring zeros by carrying last nonzero)
has_sign_crossing <- function(sgn) {
  s <- sgn
  # carry forward last nonzero to avoid "0" hiding crossings
  last <- NA_integer_
  for (i in seq_along(s)) {
    if (is.na(s[i])) next
    if (s[i] == 0L) s[i] <- last else last <- s[i]
  }
  s2 <- s[!is.na(s)]
  if (length(s2) < 2) return(FALSE)
  any(diff(s2) != 0 & !is.na(diff(s2)))
}

# ----------------------------
# 2) Load Stage5 files (if present)
# ----------------------------
paths <- list(
  B  = file.path(dirs$stage5B,  "stage5B_flip_metrics.csv"),
  B2 = file.path(dirs$stage5B2, "stage5B2_expand_local_sweep.csv"),
  C  = file.path(dirs$stage5C,  "stage5C_robustness.csv"),
  D  = file.path(dirs$stage5D,  "stage5D_ranked.csv")
)

B  <- read_csv_if_exists(paths$B)
B2 <- read_csv_if_exists(paths$B2)
C  <- read_csv_if_exists(paths$C)
D  <- read_csv_if_exists(paths$D)

# ----------------------------
# 3) Diagnose flip existence (Stage5-native labels)
# ----------------------------
flip_diag <- NULL
flip_notes <- c()

if (is.null(B)) {
  flip_notes <- c(flip_notes, "No stage5B_flip_metrics.csv found. Cannot summarize flip_ok/class_5B.")
} else {
  flip_diag <- B %>%
    mutate(candidate_id = as.character(candidate_id)) %>%
    group_by(candidate_id) %>%
    summarise(
      rF_below = first(rF_below),
      rF_above = first(rF_above),
      flip_ok  = any(flip_ok %in% TRUE, na.rm = TRUE),
      class_5B = first(class_5B),
      tainted_any = any(tainted_below %in% TRUE | tainted_above %in% TRUE, na.rm = TRUE),
      unbounded_any = any(bounded_below %in% FALSE | bounded_above %in% FALSE, na.rm = TRUE),
      .groups = "drop"
    )
  
  flip_notes <- c(flip_notes, "Flip summary uses Stage5B fields: flip_ok + class_5B (not eigen crossings).")
}

# ----------------------------
# 3B) B2-based diagnostics: cycle onset + bracket miss rate
# ----------------------------
b2_diag <- NULL
b2_notes <- c()

if (is.null(B2)) {
  b2_notes <- c(b2_notes, "No Stage5B2 file found; cannot compute cycle onset diagnostics.")
} else {
  
  # Column expectations from your Stage5B2 schema
  need_cols <- c("candidate_id","rF_center","rF_sim","bounded","has_cycle","tainted")
  missing_cols <- setdiff(need_cols, names(B2))
  if (length(missing_cols) > 0) {
    b2_notes <- c(b2_notes, paste0("Stage5B2 missing cols: ", paste(missing_cols, collapse=", ")))
  } else {
    
    # Use Stage5B bracket (rF_below/rF_above) if present, to diagnose "eps_root missed onset"
    bracket <- NULL
    if (!is.null(B) && all(c("candidate_id","rF_below","rF_above") %in% names(B))) {
      bracket <- B %>%
        transmute(candidate_id = as.character(candidate_id),
                  rF_below = as.numeric(rF_below),
                  rF_above = as.numeric(rF_above))
    }
    
    b2_work <- B2 %>%
      mutate(
        candidate_id = as.character(candidate_id),
        rF_center = as.numeric(rF_center),
        rF_sim    = as.numeric(rF_sim),
        bounded   = as.logical(bounded),
        has_cycle = as.logical(has_cycle),
        tainted   = as.logical(tainted)
      )
    
    # Onset definition: first rF_sim (closest to center in absolute distance, OR smallest rF_sim) where
    # bounded==TRUE, tainted==FALSE, has_cycle==TRUE.
    # We'll do it properly: order by |rF_sim - rF_center| (your Stage5 builds that ordering internally),
    # but also report min and max rF in sweep for context.
    b2_diag <- b2_work %>%
      group_by(candidate_id) %>%
      summarise(
        rF_center = dplyr::first(rF_center[!is.na(rF_center)]),
        n_grid = n(),
        n_bounded = sum(bounded %in% TRUE, na.rm = TRUE),
        n_cycle_ok = sum(bounded %in% TRUE & has_cycle %in% TRUE & tainted %in% FALSE, na.rm = TRUE),
        onset_rF = {
          ok <- which(bounded %in% TRUE & has_cycle %in% TRUE & tainted %in% FALSE)
          if (length(ok) == 0) NA_real_ else {
            # choose onset as the closest-to-center cycle point
            d <- abs(rF_sim[ok] - rF_center[ok])
            rF_sim[ ok[ which.min(d) ] ]
          }
        },
        onset_dist = ifelse(is.na(onset_rF) | is.na(rF_center), NA_real_, abs(onset_rF - rF_center)),
        rF_min = suppressWarnings(min(rF_sim, na.rm = TRUE)),
        rF_max = suppressWarnings(max(rF_sim, na.rm = TRUE)),
        .groups = "drop"
      )
    
    # Join bracket info if available
    if (!is.null(bracket)) {
      b2_diag <- b2_diag %>%
        left_join(bracket, by = "candidate_id") %>%
        mutate(
          onset_in_bracket = ifelse(
            is.na(onset_rF) | is.na(rF_below) | is.na(rF_above),
            NA,
            onset_rF >= pmin(rF_below, rF_above) & onset_rF <= pmax(rF_below, rF_above)
          ),
          bracket_span = ifelse(is.na(rF_below) | is.na(rF_above), NA_real_, abs(rF_above - rF_below))
        )
      b2_notes <- c(b2_notes, "Computed onset + compared to Stage5B bracket (rF_below/rF_above).")
    } else {
      b2_notes <- c(b2_notes, "Computed onset, but Stage5B bracket not found (no rF_below/rF_above).")
    }
    
    # Write CSV
    out_b2 <- file.path(dirs$stage5E, "diag_b2_onset_by_candidate.csv")
    write_csv(b2_diag, out_b2)
    b2_notes <- c(b2_notes, paste0("Wrote: ", out_b2))
  }
}


# ----------------------------
# 4) Diagnose robustness presence and (optionally) auto-join into ranked
# ----------------------------
rob_notes <- c()
augmented_written <- FALSE

if (is.null(D)) {
  rob_notes <- c(rob_notes, "No Stage5D_ranked.csv found. Cannot diagnose/join robustness.")
} else {
  # detect robustness in D
  rob_col_D <- pick_col(D, c("robust_survival", "robust_survival_rate", "robust_survival_med"))
  if (!is.null(rob_col_D)) {
    rob_notes <- c(rob_notes, paste0("Ranked already contains robustness column: ", rob_col_D))
  } else {
    rob_notes <- c(rob_notes, "Ranked DOES NOT contain robustness columns (robust_survival / robust_survival_rate).")
    if (is.null(C)) {
      rob_notes <- c(rob_notes, "But stage5C_robustness.csv not found, so join not possible.")
    } else {
      keysD <- candidate_keys(D)
      keysC <- candidate_keys(C)
      join_keys <- intersect(keysD, keysC)
      
      # detect a robustness column in C
      rob_col_C <- pick_col(C, c("robust_survival", "robust_survival_rate", "survival_rate", "survival"))
      
      if (is.null(rob_col_C)) {
        rob_notes <- c(rob_notes, "Could not detect a robustness column inside Stage5C.")
      } else if (length(join_keys) == 0) {
        rob_notes <- c(rob_notes, "No common join keys between Stage5D and Stage5C (need try_id/cand_id/id/candidate_id).")
      } else {
        # ---- FIX: cast join keys to character in BOTH tables to avoid type mismatch ----
        D_join <- D %>%
          mutate(across(all_of(join_keys), ~ as.character(.x)))
        C_join <- C %>%
          mutate(across(all_of(join_keys), ~ as.character(.x)))
        
        D_aug <- D_join %>%
          left_join(
            C_join %>% select(all_of(join_keys), robust_survival_rate = all_of(rob_col_C)),
            by = join_keys
          )
        
        out_aug <- file.path(dirs$stage5E, "stage5D_ranked_augmented.csv")
        write_csv(D_aug, out_aug)
        augmented_written <- TRUE
        
        rob_notes <- c(rob_notes, paste0(
          "Wrote augmented ranked with robustness joined: ", out_aug,
          " | join_keys=", paste(join_keys, collapse = ", "),
          " | cast_keys_to=character",
          " | C_rob_col=", rob_col_C
        ))
      }
    }
  }
}

# ----------------------------
# 5) Clamp prevalence quick check (from ranked if present)
# ----------------------------
clamp_notes <- c()
if (!is.null(D)) {
  clamp_col <- pick_col(D, c("clamp_hit", "clamp_hits"))
  if (is.null(clamp_col)) {
    clamp_notes <- c(clamp_notes, "Ranked has no clamp column (clamp_hit/clamp_hits).")
  } else {
    ch <- D[[clamp_col]]
    # normalize if logical vs integer
    clamp_rate <- if (is.logical(ch)) mean(ch, na.rm = TRUE) else mean(ch > 0, na.rm = TRUE)
    clamp_notes <- c(clamp_notes, paste0("Clamp prevalence (share with clamp>0): ", signif(clamp_rate, 4), " [col=", clamp_col, "]"))
  }
} else {
  clamp_notes <- c(clamp_notes, "No ranked file to evaluate clamps.")
}

# ----------------------------
# 6) Write outputs
# ----------------------------
report_path <- file.path(dirs$stage5E, "diag_stage5_report.txt")
flip_out    <- file.path(dirs$stage5E, "diag_flip_presence_by_candidate.csv")

lines <- c(
  "STAGE 5 DIAGNOSTICS REPORT",
  paste0("base_dir: ", base_dir),
  "",
  "FILES FOUND:",
  paste0("  5B:  ", ifelse(file.exists(paths$B),  "YES", "NO"), "  -> ", paths$B),
  paste0("  5B2: ", ifelse(file.exists(paths$B2), "YES", "NO"), "  -> ", paths$B2),
  paste0("  5C:  ", ifelse(file.exists(paths$C),  "YES", "NO"), "  -> ", paths$C),
  paste0("  5D:  ", ifelse(file.exists(paths$D),  "YES", "NO"), "  -> ", paths$D),
  "",
  "B2 (CYCLE-ONSET) DIAGNOSTICS:",
  paste0("  notes: ", if (length(b2_notes) == 0) "none" else paste(b2_notes, collapse = " | ")),
  "",
  "FLIP / HOPF CROSSING DIAGNOSTICS:",
  paste0("  using: ", flip_source),
  paste0("  notes: ", if (length(flip_notes) == 0) "none" else paste(flip_notes, collapse = " | ")),
  "",
  "ROBUSTNESS JOIN DIAGNOSTICS:",
  paste0("  notes: ", if (length(rob_notes) == 0) "none" else paste(rob_notes, collapse = " | ")),
  "",
  "CLAMP DIAGNOSTICS:",
  paste0("  notes: ", if (length(clamp_notes) == 0) "none" else paste(clamp_notes, collapse = " | "))
)

writeLines(lines, report_path)

if (!is.null(flip_diag)) {
  write_csv(flip_diag, flip_out)
} else {
  # still write a placeholder so downstream doesn't break
  write_csv(tibble(note = "flip diagnostics unavailable (missing files or columns)"), flip_out)
}

cat("\n[OK] Wrote diagnostics report:\n  ", report_path, "\n", sep = "")
cat("[OK] Wrote flip summary:\n  ", flip_out, "\n", sep = "")
if (augmented_written) cat("[OK] Wrote augmented ranked:\n  ", file.path(dirs$stage5E, "stage5D_ranked_augmented.csv"), "\n", sep = "")

# ----------------------------
# 7) Console headline summary (fast eyeballing)
# ----------------------------
if (!is.null(flip_diag)) {
  cat("\n--- Flip summary headline ---\n")
  cat("Candidates/groups checked:", nrow(flip_diag), "\n")
  cat("Has sign crossing (any):", sum(flip_diag$has_crossing, na.rm = TRUE), "\n")
  cat("Share with crossing:", signif(mean(flip_diag$has_crossing, na.rm = TRUE), 4), "\n")
}

# ----------------------------
# 7B) Console summary for B2 onset (the evaluator you actually need)
# ----------------------------
if (!is.null(b2_diag) && nrow(b2_diag) > 0) {
  cat("\n--- B2 onset headline ---\n")
  cat("Candidates in B2:", nrow(b2_diag), "\n")
  cat("With any bounded+untainted cycle in sweep:", sum(!is.na(b2_diag$onset_rF)), "\n")
  cat("Share with onset:", signif(mean(!is.na(b2_diag$onset_rF)), 4), "\n")
  
  if ("onset_in_bracket" %in% names(b2_diag)) {
    ok <- !is.na(b2_diag$onset_in_bracket) & !is.na(b2_diag$onset_rF)
    if (any(ok)) {
      cat("Onset inside Stage5B bracket (conditional on onset):",
          sum(b2_diag$onset_in_bracket[ok]), "/", sum(ok),
          "=", signif(mean(b2_diag$onset_in_bracket[ok]), 4), "\n")
    }
  }
  
  # Quick plots (base R)
  try({
    hist(b2_diag$onset_dist, main = "B2 onset distance |onset_rF - center|",
         xlab = "distance", ylab = "count")
  }, silent = TRUE)
  
  if ("bracket_span" %in% names(b2_diag)) {
    try({
      plot(b2_diag$bracket_span, b2_diag$onset_dist,
           main = "Bracket span vs onset distance",
           xlab = "bracket_span = |rF_above - rF_below|",
           ylab = "onset_dist",
           pch = ifelse(isTRUE(b2_diag$onset_in_bracket), 19, 1))
      abline(0,1,lty=2)
    }, silent = TRUE)
  }
}

# ============================================================
# 8) Read + print flip diagnostics from CSV (console) + quick viz
# ============================================================

flip_csv <- file.path(dirs$stage5E, "diag_flip_presence_by_candidate.csv")
if (!file.exists(flip_csv)) {
  cat("\n[WARN] Missing flip diagnostics CSV:\n  ", flip_csv, "\n", sep = "")
} else {
  flip_tbl <- readr::read_csv(flip_csv, show_col_types = FALSE)
  
  cat("\n==============================\n")
  cat("FLIP DIAGNOSTICS (from CSV)\n")
  cat("==============================\n")
  cat("File:\n  ", flip_csv, "\n\n", sep = "")
  
  # If it's a placeholder file, just print it
  if ("note" %in% names(flip_tbl)) {
    print(flip_tbl)
  } else {
    # Basic checks
    if (!("has_crossing" %in% names(flip_tbl))) {
      cat("[WARN] CSV has no 'has_crossing' column. Columns are:\n")
      print(names(flip_tbl))
    } else {
      n_groups <- nrow(flip_tbl)
      n_cross  <- sum(flip_tbl$has_crossing, na.rm = TRUE)
      share    <- if (n_groups > 0) n_cross / n_groups else NA_real_
      
      cat("Groups/candidates in CSV:", n_groups, "\n")
      cat("With sign-crossing:", n_cross, "\n")
      cat("Share with crossing:", signif(share, 4), "\n\n")
      
      # Show top rows with crossing and decent coverage
      # (n_points descending first)
      if ("n_points" %in% names(flip_tbl)) {
        top_cross <- flip_tbl %>%
          filter(has_crossing) %>%
          arrange(desc(n_points)) %>%
          head(12)
      } else {
        top_cross <- flip_tbl %>%
          filter(has_crossing) %>%
          head(12)
      }
      
      if (nrow(top_cross) == 0) {
        cat("No candidates with detected crossing in this sweep.\n")
        cat("Interpretation: either (i) bracket misses boundary, or (ii) metric column isn't the right one.\n\n")
      } else {
        cat("Top candidates WITH crossing (up to 12):\n")
        print(top_cross)
        cat("\n")
      }
      
      # ---- Quick plots (base R) ----
      # Only plot when key columns exist
      if (all(c("met_min", "met_max") %in% names(flip_tbl))) {
        # Range of metric: max - min
        met_range <- flip_tbl$met_max - flip_tbl$met_min
        
        cat("Plotting metric ranges (met_max - met_min) ...\n")
        
        # 1) Histogram of metric ranges
        try({
          hist(met_range, main = "Metric range across rF sweep (met_max - met_min)",
               xlab = "Metric range", ylab = "Count")
        }, silent = TRUE)
        
        # 2) Scatter met_min vs met_max with crossing highlighted by point character
        try({
          pch_vec <- ifelse(isTRUE(flip_tbl$has_crossing), 19, 1)
          plot(flip_tbl$met_min, flip_tbl$met_max,
               main = "met_min vs met_max (crossing highlighted)",
               xlab = "met_min", ylab = "met_max", pch = pch_vec)
          abline(h = 0, v = 0, lty = 2)
        }, silent = TRUE)
      }
      
      if (all(c("rF_min", "rF_max") %in% names(flip_tbl))) {
        rF_span <- flip_tbl$rF_max - flip_tbl$rF_min
        cat("Plotting rF spans (rF_max - rF_min) ...\n")
        try({
          hist(rF_span, main = "rF span per candidate/group",
               xlab = "rF_max - rF_min", ylab = "Count")
        }, silent = TRUE)
      }
    }
  }
}

