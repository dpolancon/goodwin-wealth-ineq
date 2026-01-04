# ============================================================
# 00b_wealthIneq_analytic_core.R
# Analytic core for wealth-inequality reduced 3D system
# Provides:
#   - compute_at_rF(row, rFj)   for Stage 3 analytic mode
#   - simulate_model(...)       for Stage 5 simulations
#
# HARD RULE: PHI1 LOCK (World 1)
#   - phi1 must be taken from the candidate (row$phi1_fixed OR row$phi1)
#   - do NOT endogenize phi1 inside sweeps
# ============================================================

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

as_num1 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (is.null(x) || length(x) == 0 || !is.finite(x[1])) return(NA_real_)
  x[1]
}

# Default tolerances if missing
if (!exists("eig_im_tol", inherits = TRUE)) eig_im_tol <- 1e-8

# ------------------------------------------------------------
# Default global objects (only if missing)
# ------------------------------------------------------------
if (!exists("par_base", inherits = TRUE)) {
  par_base <- list(
    kappa_min = 0.02,
    kappa_max = 0.45,
    kappa0    = 0.10,
    kappa1    = 30.0,
    phi3      = 8.0,
    phi4      = 1.0,
    phi0      = -0.02,
    alpha     = 0.02
  )
}

if (!exists("cfg", inherits = TRUE)) {
  cfg <- list(
    targets = list(e_target = 0.94),
    clamps  = list(lambda_eps = 1e-10)
  )
}

logistic <- function(x) 1 / (1 + exp(-x))

# If these are already defined in 00_globals_utils.R, we do not override.
if (!exists("kappa_fun", inherits = TRUE)) {
  kappa_fun <- function(r, p) {
    p$kappa_min + (p$kappa_max - p$kappa_min) * logistic(p$kappa1 * (r - p$kappa0))
  }
}

if (!exists("kappa_prime", inherits = TRUE)) {
  kappa_prime <- function(r, p) {
    A <- (p$kappa_max - p$kappa_min)
    z <- p$kappa1 * (r - p$kappa0)
    L <- logistic(z)
    A * p$kappa1 * L * (1 - L)
  }
}

if (!exists("Z_fun", inherits = TRUE)) {
  Z_fun <- function(d, f, p) logistic(p$phi3 * ((d - 1) + p$phi4 * (f - 1)))
}

if (!exists("lambda_fun", inherits = TRUE)) {
  lambda_fun <- function(r, rF, psi) logistic(psi * (r - rF))
}

# ------------------------------------------------------------
# compute_at_rF(row, rFj)
#   row: list or 1-row tibble/data.frame (candidate)
#   rFj: scalar r_F for Hopf metric and stability at
# ------------------------------------------------------------
compute_at_rF <- function(row, rFj) {
  
  if (is.data.frame(row)) row <- as.list(row[1, , drop = FALSE])
  if (!is.list(row))      row <- as.list(row)
  
  p <- par_base
  if (!is.null(row$kappa_max) && is.finite(as_num1(row$kappa_max))) p$kappa_max <- as_num1(row$kappa_max)
  if (!is.null(row$kappa_min) && is.finite(as_num1(row$kappa_min))) p$kappa_min <- as_num1(row$kappa_min)
  
  sigma <- as_num1(row$sigma)
  g_n   <- as_num1(row$g_n)
  i     <- as_num1(row$i)
  delta <- as_num1(row$delta)
  psi   <- as_num1(row$psi)
  phi2  <- as_num1(row$phi2)
  
  if (!is.finite(g_n) || g_n <= 0) return(list(ok=FALSE, reason="bad_gn"))
  if (!is.finite(sigma) || sigma <= 0) return(list(ok=FALSE, reason="bad_sigma"))
  
  omega <- as_num1(row$omega_star)
  d     <- as_num1(row$d_star)
  if (!is.finite(omega) || !is.finite(d)) return(list(ok=FALSE, reason="missing_omega_or_d"))
  
  # implied r from identity
  r <- (1 - omega - i*d) / sigma
  if (!is.finite(r)) return(list(ok=FALSE, reason="nonfinite_r"))
  
  kap <- kappa_fun(r, p)
  kp  <- kappa_prime(r, p)
  if (!is.finite(kap) || !is.finite(kp)) return(list(ok=FALSE, reason="nonfinite_kappa_or_kp"))
  
  g  <- kap / sigma - delta
  gr <- kp  / sigma
  
  lam   <- lambda_fun(r, rFj, psi)
  denom <- 1 - lam
  if (!is.finite(lam) || !is.finite(denom) || abs(denom) < 1e-12) {
    return(list(ok=FALSE, reason="lambda_denominator_near_zero_or_nonfinite"))
  }
  
  q     <- lam / denom
  iotaF <- r * q
  f     <- iotaF / g_n
  if (!is.finite(f)) return(list(ok=FALSE, reason="nonfinite_f_star"))
  
  Z <- Z_fun(d, f, p)
  if (!is.finite(Z)) return(list(ok=FALSE, reason="nonfinite_Z_star"))
  
  # ----------------------------
  # PHI1 LOCK (World 1)
  # ----------------------------
  phi1 <- as_num1(row$phi1_fixed)
  if (!is.finite(phi1) || phi1 <= 0) phi1 <- as_num1(row$phi1)
  if (!is.finite(phi1) || phi1 <= 0) return(list(ok=FALSE, reason="missing_or_bad_phi1_locked"))
  
  # implied e*
  e_star <- (p$alpha - p$phi0 + phi2 * Z) / phi1
  if (!is.finite(e_star)) return(list(ok=FALSE, reason="nonfinite_e_star"))
  
  # ----------------------------
  # Derivatives and Jacobian
  # ----------------------------
  eps <- 1e-6
  lam_p <- lambda_fun(r + eps, rFj, psi)
  lam_m <- lambda_fun(r - eps, rFj, psi)
  if (!is.finite(lam_p) || !is.finite(lam_m)) return(list(ok=FALSE, reason="nonfinite_lambda_fd"))
  lam_r <- (lam_p - lam_m) / (2*eps)
  
  q_r     <- lam_r / (denom^2)
  iotaF_r <- q + r*q_r
  f_r     <- iotaF_r / g_n
  if (!is.finite(f_r)) return(list(ok=FALSE, reason="nonfinite_f_r"))
  
  Z_s <- Z * (1 - Z)
  
  dr_domega <- -1 / sigma
  dr_dd     <- -i / sigma
  
  ds_domega <- p$phi3 * (p$phi4 * f_r * dr_domega)
  ds_dd     <- p$phi3 * (1 + p$phi4 * f_r * dr_dd)
  
  Z_omega <- Z_s * ds_domega
  Z_d     <- Z_s * ds_dd
  
  J11 <- (g - g_n)
  J12 <- e_star * gr * dr_domega
  J13 <- e_star * gr * dr_dd
  
  J21 <- omega * phi1
  J22 <- omega * (-phi2 * Z_omega)
  J23 <- omega * (-phi2 * Z_d)
  
  J31 <- 0
  J32 <- 1 + kp * dr_domega - d * gr * dr_domega
  J33 <- i - g + kp * dr_dd - d * gr * dr_dd
  
  J <- matrix(c(J11,J12,J13,
                J21,J22,J23,
                J31,J32,J33), nrow=3, byrow=TRUE)
  
  if (any(!is.finite(J))) return(list(ok=FALSE, reason="nonfinite_jacobian"))
  
  ev <- tryCatch(eigen(J)$values, error=function(e) NA_complex_)
  if (all(is.na(ev))) return(list(ok=FALSE, reason="eigen_failed"))
  
  maxRe <- max(Re(ev))
  maxIm <- max(abs(Im(ev)))
  
  # For det(λI - J) = λ^3 + a1 λ^2 + a2 λ + a3
  a1 <- -sum(diag(J))
  a2 <- (J11*J22 - J12*J21) + (J11*J33 - J13*J31) + (J22*J33 - J23*J32)
  a3 <- -det(J)
  hopf_val <- a1*a2 - a3
  
  RH_ok <- is.finite(a1) && is.finite(a2) && is.finite(a3) && is.finite(hopf_val) &&
    (a1 > 0 && a2 > 0 && a3 > 0)
  
  list(
    ok=TRUE,
    e_star=e_star, omega_star=omega, d_star=d,
    r=r, lambda_star=lam, f_star=f, Z_star=Z,
    phi1_fixed=phi1,
    a1=a1, a2=a2, a3=a3,
    hopf_val=hopf_val, RH_ok=RH_ok,
    maxRe=maxRe, maxIm=maxIm,
    stable = (is.finite(maxRe) && maxRe < 0),
    has_complex = (is.finite(maxIm) && maxIm > eig_im_tol)
  )
}

# ============================================================
# Simulation hook for Stage 5 (reduced 3D) – PHI1 LOCKED
# ============================================================
simulate_model <- function(row, rF_sim, t_end, dt, x0 = NULL) {
  if (!requireNamespace("deSolve", quietly = TRUE)) stop("simulate_model(): requires 'deSolve'.")
  if (!requireNamespace("tibble", quietly = TRUE)) stop("simulate_model(): requires 'tibble'.")
  if (!requireNamespace("dplyr", quietly = TRUE))  stop("simulate_model(): requires 'dplyr'.")
  
  if (is.data.frame(row)) row <- as.list(row[1, , drop = FALSE])
  if (!is.list(row)) stop("simulate_model(): 'row' must be list-like.")
  
  p <- par_base
  if (!is.null(row$kappa_max) && is.finite(as_num1(row$kappa_max))) p$kappa_max <- as_num1(row$kappa_max)
  if (!is.null(row$kappa_min) && is.finite(as_num1(row$kappa_min))) p$kappa_min <- as_num1(row$kappa_min)
  
  sigma <- as_num1(row$sigma)
  g_n   <- as_num1(row$g_n)
  i     <- as_num1(row$i)
  delta <- as_num1(row$delta)
  psi   <- as_num1(row$psi)
  phi2  <- as_num1(row$phi2)
  
  if (!is.finite(sigma) || !is.finite(g_n) || !is.finite(i) || !is.finite(delta) ||
      !is.finite(psi)   || !is.finite(phi2) || sigma <= 0 || g_n <= 0) {
    out <- tibble::tibble(time=NA_real_, e=NA_real_, omega=NA_real_, d=NA_real_)
    attr(out, "stop_reason") <- "bad_params"
    return(out)
  }
  
  # PHI1 LOCK: use candidate phi1_fixed OR phi1; otherwise treat as invalid
  phi1 <- as_num1(row$phi1_fixed)
  if (!is.finite(phi1) || phi1 <= 0) phi1 <- as_num1(row$phi1)
  if (!is.finite(phi1) || phi1 <= 0) {
    out <- tibble::tibble(time=NA_real_, e=NA_real_, omega=NA_real_, d=NA_real_)
    attr(out, "stop_reason") <- "missing_phi1_locked"
    return(out)
  }
  
  e_target <- cfg$targets$e_target %||% 0.94
  lam_eps  <- cfg$clamps$lambda_eps %||% 1e-10
  
  # initial conditions
  if (is.null(x0)) {
    e0 <- if (is.finite(as_num1(row$e_star))) as_num1(row$e_star) else
      if (is.finite(as_num1(row$e_implied))) as_num1(row$e_implied) else e_target
    w0 <- if (is.finite(as_num1(row$omega_star))) as_num1(row$omega_star) else 0.65
    d0 <- if (is.finite(as_num1(row$d_star))) as_num1(row$d_star) else 0.50
    x0 <- c(e=e0, omega=w0, d=d0)
  } else {
    x0 <- as.numeric(x0)
    if (is.null(names(x0)) || !all(c("e","omega","d") %in% names(x0))) {
      out <- tibble::tibble(time=NA_real_, e=NA_real_, omega=NA_real_, d=NA_real_)
      attr(out, "stop_reason") <- "bad_x0"
      return(out)
    }
  }
  
  rhs <- function(t, state, parms) {
    e     <- state[["e"]]
    omega <- state[["omega"]]
    d     <- state[["d"]]
    
    r   <- (1 - omega - i*d) / sigma
    kap <- kappa_fun(r, p)
    g   <- kap / sigma - delta
    
    lam <- lambda_fun(r, rF = rF_sim, psi = psi)
    lam <- min(max(lam, lam_eps), 1 - lam_eps)
    
    iotaF <- r * (lam / (1 - lam))
    f     <- iotaF / g_n
    
    Z <- Z_fun(d, f, p)
    
    de     <- (g - g_n) * e
    domega <- omega * (p$phi0 + phi1 * e - p$alpha - phi2 * Z)
    dd     <- kap - (1 - omega) + i*d - d*g
    
    list(c(de, domega, dd))
  }
  
  times <- seq(0, t_end, by = dt)
  
  sol <- tryCatch(
    deSolve::ode(y=x0, times=times, func=rhs, parms=NULL, method="lsoda"),
    error = function(e) {
      out <- tibble::tibble(time=times, e=NA_real_, omega=NA_real_, d=NA_real_)
      attr(out, "stop_reason") <- paste0("ode_error:", conditionMessage(e))
      out
    }
  )
  
  if (inherits(sol, "tbl_df")) return(sol)
  
  sol <- as.data.frame(sol)
  out <- tibble::as_tibble(sol) %>% dplyr::select("time","e","omega","d")
  attr(out, "stop_reason") <- "ok"
  out
}
