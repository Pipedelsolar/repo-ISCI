library(splines2)
library(survival)
library(dplyr)
library(ggplot2)
library(readr)
library(scales)

# ============================================================
# 0) DATOS
# ============================================================

setwd("C:/Users/ntrig/Desktop/ISCI/Proyectos/Efectividad_Nirse/Code")

df_raw <- read_csv(
  "../Data/df_tv_matched_11_03_2026.csv",
  col_types = cols(
    DIAG9  = col_character(),
    DIAG10 = col_character(),
    DIAG11 = col_character()
  ),
  show_col_types = FALSE
)

df_raw$RUN        <- as.character(df_raw$RUN)
df_raw$Group      <- as.factor(df_raw$Group)
df_raw$start      <- as.numeric(df_raw$start)
df_raw$stop       <- as.numeric(df_raw$stop)
df_raw$event_vrs  <- as.integer(df_raw$event_vrs)
df_raw$inmunizado <- as.integer(df_raw$inmunizado)
df_raw            <- df_raw %>% filter(stop > start)

stopifnot(all(df_raw$event_vrs  %in% c(0, 1)))
stopifnot(all(df_raw$inmunizado %in% c(0, 1)))

# Tiempo de inmunización
if (!"t_inm" %in% names(df_raw)) {
  ref_date        <- as.Date("2024-04-01")
  df_raw$fechaInm <- as.Date(df_raw$fechaInm)
  df_raw$t_inm    <- as.numeric(df_raw$fechaInm - ref_date)
}
df_raw$t_inm <- as.numeric(df_raw$t_inm)

# tau = meses desde vacunación
df_raw <- df_raw %>%
  mutate(
    tau_start_days = ifelse(inmunizado == 1, pmax(0, start - t_inm), 0),
    tau_stop_days  = ifelse(inmunizado == 1, pmax(0, stop  - t_inm), 0),
    tau_eval_days  = ifelse(inmunizado == 1,
                            0.5 * (tau_start_days + tau_stop_days), 0),
    tau_eval_m     = pmax(0, tau_eval_days / 30.4375)
  )

cat("── Resumen ──\n")
cat("Filas:", nrow(df_raw), "\n")
cat("Personas:", n_distinct(df_raw$RUN), "\n")
cat("Eventos:", sum(df_raw$event_vrs), "\n")
cat("Eventos inmunizados:", sum(df_raw$event_vrs[df_raw$inmunizado == 1]), "\n")

# ============================================================
# 1) VERIFICAR PROPIEDAD CLAVE DE LA NCS
# ============================================================
# En la parametrización splines2::naturalSpline con intercept=TRUE:
# - Columna 1 = intercepto (no-cero en todo el rango, constante más allá de tau*)
# - Última columna = "forward extrapolation variable" (no-cero solo más allá de tau*)
# - Columnas intermedias = 0 más allá del boundary knot superior
#
# Si fijamos gamma_intercept = 0 y gamma_fev = 0,
# la spline es exactamente 0 para tau > tau_star

cat("\n── Verificando propiedad NCS ──\n")
tau_test <- c(5, 10, 14.9, 15.001, 16, 18)
B_test   <- as.matrix(splines2::naturalSpline(
  tau_test,
  knots          = c(5, 10),
  Boundary.knots = c(0, 15),
  intercept      = TRUE
))
cat("Base NCS evaluada en tau = ", tau_test, ":\n")
print(round(B_test, 5))
cat("\nPara tau > 15 (filas 4-6), solo col 1 y última deben ser != 0:\n")
print(round(B_test[4:6, ], 5))

# ============================================================
# 2) FUNCIONES DE CONSTRUCCIÓN
# ============================================================

make_ncs_jennings <- function(tau, knots_internal, tau_star_m) {
  as.matrix(splines2::naturalSpline(
    tau,
    knots          = knots_internal,
    Boundary.knots = c(0, tau_star_m),
    intercept      = TRUE
  ))
}

build_design_matrix <- function(df, knots_internal, tau_star_m) {
  # tau: usar tau_eval_m pero sin clipear — la spline se encarga de tau > tau*
  tau <- df$tau_eval_m
  inm <- as.numeric(df$inmunizado == 1)

  # Base NCS multiplicada por indicador de inmunizado
  B       <- make_ncs_jennings(tau, knots_internal, tau_star_m)
  n_basis <- ncol(B)

  # Diseño: [inmunizado, inmunizado * B1, ..., inmunizado * Bk]
  # - inmunizado: efecto PH de referencia (beta_inm)
  # - inmunizado * B: efecto time-varying via NCS
  X <- cbind(
    inmunizado = inm,
    sweep(B, 1, inm, `*`)   # inm * B columna a columna
  )

  # Nombres de columnas
  col_names <- c(
    "beta_inm",
    "gamma_intercept",                          # col 1 NCS — se fija a 0
    paste0("gamma_", seq_len(n_basis - 2L)),    # columnas internas — libres
    "gamma_fev"                                 # última col NCS — se fija a 0
  )
  colnames(X) <- col_names

  list(X = X, n_basis = n_basis, col_names = col_names)
}

# ============================================================
# 3) LIKELIHOOD + GRADIENTE + HESSIAN (vectorizado)
# ============================================================

cox_loglik_grad_hess <- function(theta, X, start, stop, event, strata) {

  p       <- length(theta)
  eta     <- as.numeric(X %*% theta)
  exp_eta <- exp(eta)

  ll   <- 0
  grad <- numeric(p)
  H    <- matrix(0, p, p)

  for (s in unique(strata)) {
    idx_s <- which(strata == s)
    s_    <- start[idx_s]
    t_    <- stop[idx_s]
    ev_   <- event[idx_s]
    ee_   <- exp_eta[idx_s]
    X_    <- X[idx_s, , drop = FALSE]

    event_times <- sort(unique(t_[ev_ == 1]))

    for (t in event_times) {
      ev_idx   <- which(t_ == t & ev_ == 1)
      risk_idx <- which(s_ < t & t_ >= t)
      if (length(risk_idx) == 0) next

      ee_r <- ee_[risk_idx]
      X_r  <- X_[risk_idx, , drop = FALSE]

      S0 <- sum(ee_r)
      S1 <- colSums(ee_r * X_r)
      S2 <- crossprod(X_r, ee_r * X_r)

      n_ev <- length(ev_idx)
      ll   <- ll   + sum(as.numeric(X_[ev_idx, , drop = FALSE] %*% theta)) -
                     n_ev * log(S0)
      grad <- grad + colSums(X_[ev_idx, , drop = FALSE]) - n_ev * S1 / S0
      H    <- H    - n_ev * (S2 / S0 - tcrossprod(S1 / S0))
    }
  }

  list(ll = ll, grad = grad, hess = H)
}

# ============================================================
# 4) NEWTON-RAPHSON CON CONSTRAINTS DIRECTOS
# ============================================================

fit_cox_nr <- function(df_prepared,
                        knots_internal,
                        tau_star_m  = 15,
                        tau_grid_by = 0.05,
                        max_iter    = 50,
                        tol         = 1e-8,
                        verbose     = TRUE) {

  df <- df_prepared

  # Construir diseño
  dm      <- build_design_matrix(df, knots_internal, tau_star_m)
  X       <- dm$X
  n_total <- ncol(X)

  # Indices constrained (fijos a 0) y libres
  constrained_idx <- c(2L, as.integer(n_total))
  free_idx        <- setdiff(seq_len(n_total), constrained_idx)
  n_free          <- length(free_idx)

  if (verbose) {
    cat(sprintf("Parámetros totales : %d\n", n_total))
    cat(sprintf("Constrained a 0   : %s\n",
                paste(colnames(X)[constrained_idx], collapse = ", ")))
    cat(sprintf("Parámetros libres  : %d\n", n_free))
  }

  # Warm start: modelo PH constante
  fit_init <- coxph(
    Surv(start, stop, event_vrs) ~ inmunizado + strata(Group) + cluster(RUN),
    data = df, ties = "efron", robust = TRUE
  )
  alpha    <- rep(0, n_free)
  idx_inm  <- match("beta_inm", colnames(X)[free_idx])
  if (!is.na(idx_inm)) alpha[idx_inm] <- as.numeric(coef(fit_init)["inmunizado"])

  ll_prev <- -Inf

  # ── Iteraciones Newton-Raphson ──
  for (iter in seq_len(max_iter)) {

    theta        <- numeric(n_total)
    theta[free_idx] <- alpha

    lgh <- cox_loglik_grad_hess(
      theta  = theta,
      X      = X,
      start  = df$start,
      stop   = df$stop,
      event  = df$event_vrs,
      strata = df$Group
    )

    grad_free <- lgh$grad[free_idx]
    hess_free <- lgh$hess[free_idx, free_idx, drop = FALSE]

    # Paso Newton: delta = -H^{-1} g
    step <- tryCatch(
      solve(hess_free, grad_free),
      error = function(e) {
        if (verbose) message("  Hessian singular — gradiente con paso pequeño")
        grad_free * 1e-3
      }
    )

    # Backtracking line search
    lr        <- 1.0
    alpha_new <- alpha - lr * step

    for (ls in seq_len(25)) {
      theta_try        <- numeric(n_total)
      theta_try[free_idx] <- alpha_new
      ll_try <- tryCatch(
        cox_loglik_grad_hess(
          theta  = theta_try,
          X      = X,
          start  = df$start,
          stop   = df$stop,
          event  = df$event_vrs,
          strata = df$Group
        )$ll,
        error = function(e) -Inf
      )
      if (is.finite(ll_try) && ll_try > lgh$ll) break
      lr        <- lr * 0.5
      alpha_new <- alpha - lr * step
    }

    delta_ll <- lgh$ll - ll_prev

    if (verbose) {
      cat(sprintf("  Iter %2d | ll = %10.3f | Δll = %+.6f | lr = %.4f\n",
                  iter, lgh$ll, delta_ll, lr))
    }

    if (iter > 1 && abs(delta_ll) < tol) {
      if (verbose) cat(sprintf("  ✓ Convergió en iteración %d\n", iter))
      break
    }

    alpha   <- alpha_new
    ll_prev <- lgh$ll
  }

  # ── Estimados finales ──
  theta_final        <- numeric(n_total)
  theta_final[free_idx] <- alpha
  names(theta_final) <- colnames(X)

  final_lgh <- cox_loglik_grad_hess(
    theta  = theta_final,
    X      = X,
    start  = df$start,
    stop   = df$stop,
    event  = df$event_vrs,
    strata = df$Group
  )

  vcov_free <- tryCatch(
    solve(-final_lgh$hess[free_idx, free_idx]),
    error = function(e) {
      if (verbose) message("  Advertencia: hessian no invertible — SE no disponibles")
      NULL
    }
  )

  # ── Predicción en grilla ──
  tau_grid <- seq(0, tau_star_m, by = tau_grid_by)

  df_grid <- data.frame(
    inmunizado = 1L,
    tau_eval_m = tau_grid
  )
  X_grid <- build_design_matrix(df_grid, knots_internal, tau_star_m)$X
  X_grid[, constrained_idx] <- 0   # forzar constraints

  eta <- as.numeric(X_grid %*% theta_final)
  eta[tau_grid >= tau_star_m] <- 0   # HR = 1 después de tau*

  if (!is.null(vcov_free)) {
    X_gf   <- X_grid[, free_idx, drop = FALSE]
    se_eta <- sqrt(pmax(0, rowSums((X_gf %*% vcov_free) * X_gf)))
    se_eta[tau_grid >= tau_star_m] <- 0
  } else {
    se_eta <- rep(NA_real_, length(tau_grid))
  }

  HR <- exp(eta)
  VE <- 1 - HR

  ve_curve <- data.frame(
    tau_m = tau_grid,
    logHR = eta,
    se    = se_eta,
    HR    = HR,
    HR_lo = exp(eta - 1.96 * se_eta),
    HR_hi = exp(eta + 1.96 * se_eta),
    VE    = VE,
    VE_lo = 1 - exp(eta + 1.96 * se_eta),
    VE_hi = 1 - exp(eta - 1.96 * se_eta)
  )

  # ── Model summary ──
  n_events <- sum(df$event_vrs)
  fit_const <- coxph(
    Surv(start, stop, event_vrs) ~ inmunizado + strata(Group) + cluster(RUN),
    data = df, ties = "efron", robust = TRUE
  )

  model_summary <- data.frame(
    n_knots_int  = length(knots_internal),
    knot_pos     = ifelse(length(knots_internal) == 0, "—",
                          paste(round(knots_internal, 1), collapse = ", ")),
    n_free       = n_free,
    n_events     = n_events,
    EPV          = round(n_events / n_free, 1),
    loglik       = round(final_lgh$ll, 3),
    loglik_const = round(fit_const$loglik[2], 3),
    LRT_stat     = round(2 * (final_lgh$ll - fit_const$loglik[2]), 3),
    LRT_p        = pchisq(
      2 * (final_lgh$ll - fit_const$loglik[2]),
      df = n_free - 1, lower.tail = FALSE
    ),
    AIC          = round(-2 * final_lgh$ll + 2 * n_free, 2)
  )

  list(
    theta_hat       = theta_final,
    vcov_free       = vcov_free,
    free_idx        = free_idx,
    constrained_idx = constrained_idx,
    loglik          = final_lgh$ll,
    ve_curve        = ve_curve,
    model_summary   = model_summary,
    knots_internal  = knots_internal,
    tau_star_m      = tau_star_m
  )
}

# ============================================================
# 5) PREPARAR DATOS Y CORRER MODELOS
# ============================================================

# Percentiles de tau para definir knots
tau_imm_all <- df_raw %>%
  filter(inmunizado == 1, tau_eval_m < 15, is.finite(tau_eval_m)) %>%
  pull(tau_eval_m)

q <- function(p) as.numeric(quantile(tau_imm_all, p))

cat("\nPercentiles de tau (brazo inmunizado):\n")
print(round(quantile(tau_imm_all, probs = seq(0.1, 0.9, by = 0.1)), 2))



# ============================================================
# MODELO COMPLETO:
# log-HR(tau) = beta0 + beta1*tau + sum_j gamma_j * B_j(tau)
#
# Con restricciones:
#   R1: log-HR(tau_star) = 0   (valor)
#   R2: d/dtau log-HR(tau_star) = 0  (pendiente)
#   R3: log-HR(tau) <= 0 para todo tau  (VE >= 0)
# ============================================================

build_design_matrix_full <- function(df, knots_internal, tau_star_m) {

  tau <- pmin(df$tau_eval_m, tau_star_m)
  inm <- as.numeric(df$inmunizado == 1)

  # Base NCS residualizada (igual que antes)
  # Pero aquí usamos naturalSpline sin intercept para las gammas
  # porque beta0 y beta1 ya capturan la parte lineal
  B_raw <- as.matrix(splines2::naturalSpline(
    tau,
    knots          = knots_internal,
    Boundary.knots = c(0, tau_star_m),
    intercept      = FALSE   # sin intercepto — la parte lineal va separada
  ))

  # Residualizar B contra {1, tau} para que gamma capture solo no-linealidad
  tau_imm <- tau[inm == 1 & df$tau_eval_m < tau_star_m]
  B_imm   <- B_raw[inm == 1 & df$tau_eval_m < tau_star_m, , drop = FALSE]

  M_imm <- cbind(1, tau_imm)
  XtX   <- crossprod(M_imm)
  proj  <- matrix(0, 2, ncol(B_raw))
  for (j in seq_len(ncol(B_raw))) {
    proj[, j] <- solve(XtX, crossprod(M_imm, B_imm[, j]))
  }

  M_all  <- cbind(1, tau)
  B_nl   <- B_raw
  for (j in seq_len(ncol(B_raw))) {
    B_nl[, j] <- B_raw[, j] - M_all %*% proj[, j]
  }

  # Diseño: beta0*inm + beta1*inm*tau + gamma_j*inm*B_nl_j
  X <- cbind(
    beta0  = inm,
    beta1  = inm * tau,
    sweep(B_nl, 1, inm, `*`)
  )
  colnames(X) <- c("beta0", "beta1",
                   paste0("gamma_", seq_len(ncol(B_nl))))

  list(X = X, B_nl = B_nl, proj = proj,
       knots_internal = knots_internal,
       tau_star_m = tau_star_m)
}

eval_nl_basis_at <- function(tau_val, knots_internal, tau_star_m, proj) {
  tau_clip <- pmin(pmax(tau_val, 0), tau_star_m)
  B_raw <- as.matrix(splines2::naturalSpline(
    tau_clip,
    knots          = knots_internal,
    Boundary.knots = c(0, tau_star_m),
    intercept      = FALSE
  ))
  M <- cbind(1, tau_clip)
  B_nl <- B_raw
  for (j in seq_len(ncol(B_raw))) {
    B_nl[, j] <- B_raw[, j] - M %*% proj[, j]
  }
  B_nl
}

fit_cox_nr_full <- function(df_prepared,
                             knots_internal,
                             tau_star_m   = 15,
                             tau_grid_by  = 0.05,
                             max_iter     = 50,
                             tol          = 1e-8,
                             lambda_ineq  = 0,      # 0 = sin constraint VE>=0
                             lambda_ridge=0,
                             alpha_init = NULL,
                             n_check      = 100,    # puntos grilla para constraint
                             verbose      = TRUE) {

  df <- df_prepared

  dm      <- build_design_matrix_full(df, knots_internal, tau_star_m)
  X       <- dm$X
  proj    <- dm$proj
  n_total <- ncol(X)

  # ── Restricciones de igualdad R1 y R2 via null-space ──
  eps         <- 1e-6
  tau_at_star <- tau_star_m - eps

  B_star   <- eval_nl_basis_at(tau_at_star,       knots_internal, tau_star_m, proj)
  B_star_p <- eval_nl_basis_at(tau_at_star + eps,  knots_internal, tau_star_m, proj)
  B_star_m <- eval_nl_basis_at(tau_at_star - eps,  knots_internal, tau_star_m, proj)
  dB_star  <- (B_star_p - B_star_m) / (2 * eps)

  A <- rbind(
    c(1, tau_at_star, as.numeric(B_star)),
    c(0, 1,           as.numeric(dB_star))
  )
  colnames(A) <- colnames(X)

  svd_A  <- svd(A, nu = 0, nv = n_total)
  r      <- sum(svd_A$d > 1e-10 * max(svd_A$d))
  Nmat   <- svd_A$v[, (r + 1):n_total, drop = FALSE]
  n_free <- ncol(Nmat)
  Z      <- X %*% Nmat
  colnames(Z) <- paste0("z", seq_len(n_free))

    alpha <- if (!is.null(alpha_init) && length(alpha_init) == n_free) {
    alpha_init
  } else {
    rep(0, n_free)
  }

  # ── Grilla para evaluar constraint VE >= 0 (log-HR <= 0) ──
  # Evaluamos en puntos densos de [0, tau_star)
  tau_chk     <- seq(0.01, tau_star_m - 0.01, length.out = n_check)
  B_chk_nl    <- eval_nl_basis_at(tau_chk, knots_internal, tau_star_m, proj)
  active_chk  <- rep(1, n_check)
  Xchk_orig   <- cbind(
    beta0 = active_chk,
    beta1 = active_chk * tau_chk,
    sweep(B_chk_nl, 1, active_chk, `*`)
  )
  colnames(Xchk_orig) <- colnames(X)
  # Mapear al espacio null-space
  Zchk <- Xchk_orig %*% Nmat   # [n_check x n_free]

  if (verbose) {
    cat(sprintf("Parámetros originales  : %d\n", n_total))
    cat(sprintf("Parámetros libres      : %d\n", n_free))
    cat(sprintf("Penalización VE>=0     : lambda = %g\n", lambda_ineq))
    cat(sprintf("Puntos grilla constraint: %d\n", n_check))
  }

  # ── Warm start ──
  fit_init <- coxph(
    Surv(start, stop, event_vrs) ~ inmunizado + strata(Group) + cluster(RUN),
    data = df, ties = "efron", robust = TRUE
  )
  alpha <- rep(0, n_free)

  # ── Función que evalúa penalización y su gradiente/hessian ──
  eval_penalty <- function(alpha) {

    # ── Penalización R3 (VE >= 0) ──
    pen_r3 <- 0; grad_r3 <- numeric(n_free); hess_r3 <- matrix(0, n_free, n_free)
    if (lambda_ineq > 0) {
      eta_chk     <- as.numeric(Zchk %*% alpha)
      viol        <- pmax(0, eta_chk)
      active_viol <- viol > 0
      pen_r3      <- lambda_ineq * sum(viol^2)
      grad_r3     <- 2 * lambda_ineq * as.numeric(t(Zchk) %*% viol)
      if (any(active_viol)) {
        hess_r3 <- 2 * lambda_ineq * crossprod(Zchk[active_viol, , drop=FALSE])
      }
    }

    # ── Penalización Ridge sobre gamma (estabiliza con muchos knots) ──
    pen_ridge <- 0; grad_ridge <- numeric(n_free); hess_ridge <- matrix(0, n_free, n_free)
    if (lambda_ridge > 0) {
      # Solo penalizar los gammas, no beta0 y beta1
      # En el espacio alpha es más difícil separar, así que penalizamos todo alpha
      # excepto el primer componente (que corresponde principalmente a beta0/beta1)
      pen_ridge   <- lambda_ridge * sum(alpha^2)
      grad_ridge  <- 2 * lambda_ridge * alpha
      hess_ridge  <- 2 * lambda_ridge * diag(n_free)
    }

    list(
      pen      = pen_r3 + pen_ridge,
      grad_pen = grad_r3 + grad_ridge,
      hess_pen = hess_r3 + hess_ridge
    )
  }

  # ── Newton-Raphson con penalización ──
  ll_prev    <- -Inf
  pen_prev   <- Inf

  for (iter in seq_len(max_iter)) {

    # Likelihood + gradiente + hessian
    lgh <- cox_loglik_grad_hess(
      theta  = alpha,
      X      = Z,
      start  = df$start,
      stop   = df$stop,
      event  = df$event_vrs,
      strata = df$Group
    )

    # Penalización
    pen_obj <- eval_penalty(alpha)

    # Objetivo total: maximizar ll - pen
    obj_total  <- lgh$ll       - pen_obj$pen
    grad_total <- lgh$grad     - pen_obj$grad_pen
    hess_total <- lgh$hess     - pen_obj$hess_pen   # hess es negativo semidefinido

    # Paso Newton
    step <- tryCatch(
      solve(hess_total, grad_total),
      error = function(e) {
        if (verbose) message("  Hessian singular — gradiente con paso pequeño")
        -grad_total * 1e-3
      }
    )

    # Backtracking line search sobre objetivo total
    lr        <- 1.0
    alpha_new <- alpha - lr * step

    for (ls in seq_len(30)) {
      lgh_try <- tryCatch(
        cox_loglik_grad_hess(
          theta = alpha_new, X = Z,
          start = df$start, stop = df$stop,
          event = df$event_vrs, strata = df$Group
        )$ll,
        error = function(e) -Inf
      )
      pen_try <- eval_penalty(alpha_new)$pen
      obj_try <- lgh_try - pen_try

      if (is.finite(obj_try) && obj_try > obj_total) break
      lr        <- lr * 0.5
      alpha_new <- alpha - lr * step
    }

    delta_obj <- obj_total - (ll_prev - pen_prev)

    if (verbose) {
      n_viol <- if (lambda_ineq > 0) {
        sum(as.numeric(Zchk %*% alpha) > 1e-6)
      } else 0
      cat(sprintf(
        "  Iter %2d | ll = %9.3f | pen = %7.3f | obj = %9.3f | viol = %d | lr = %.3f\n",
        iter, lgh$ll, pen_obj$pen, obj_total, n_viol, lr
      ))
    }

    if (iter > 1 && abs(delta_obj) < tol) {
      if (verbose) cat(sprintf("  ✓ Convergió en iteración %d\n", iter))
      break
    }

    alpha    <- alpha_new
    ll_prev  <- lgh$ll
    pen_prev <- pen_obj$pen
  }

  # ── Verificar violaciones finales ──
  eta_final_chk <- as.numeric(Zchk %*% alpha)
  n_viol_final  <- sum(eta_final_chk > 1e-6)
  if (verbose) {
    cat(sprintf("\nViolaciones finales VE < 0: %d / %d puntos\n",
                n_viol_final, n_check))
    if (n_viol_final > 0) {
      cat("  Rango de violaciones (log-HR):",
          round(range(eta_final_chk[eta_final_chk > 1e-6]), 4), "\n")
    }
  }

  # ── Recuperar theta original y varianza ──
  theta_hat        <- as.vector(Nmat %*% alpha)
  names(theta_hat) <- colnames(X)

  final_lgh <- cox_loglik_grad_hess(
    theta = alpha, X = Z,
    start = df$start, stop = df$stop,
    event = df$event_vrs, strata = df$Group
  )

  # Varianza: solo del término likelihood (no penalizado)
  # porque la penalización es un dispositivo numérico, no parte del modelo
  vcov_alpha <- tryCatch(
    solve(-final_lgh$hess),
    error = function(e) NULL
  )
  V_theta <- if (!is.null(vcov_alpha)) Nmat %*% vcov_alpha %*% t(Nmat) else NULL

  # ── Predicción ──
  tau_grid    <- seq(0, tau_star_m, by = tau_grid_by)
  B_grid_nl   <- eval_nl_basis_at(tau_grid, knots_internal, tau_star_m, proj)
  active_grid <- as.numeric(tau_grid < tau_star_m)

  Xg <- cbind(
    beta0 = active_grid,
    beta1 = active_grid * tau_grid,
    sweep(B_grid_nl, 1, active_grid, `*`)
  )
  colnames(Xg) <- colnames(X)

  eta <- as.numeric(Xg %*% theta_hat)
  eta[tau_grid >= tau_star_m] <- 0

  if (!is.null(V_theta)) {
    se_eta <- sqrt(pmax(0, rowSums((Xg %*% V_theta) * Xg)))
    near   <- tau_grid >= (tau_star_m - 1)
    w      <- (tau_grid[near] - (tau_star_m - 1)) / 1
    se_eta[near] <- se_eta[near] * (1 - w)^2
    se_eta[tau_grid >= tau_star_m] <- 0
  } else {
    se_eta <- rep(NA_real_, length(tau_grid))
  }

  HR <- exp(eta); VE <- 1 - HR

  ve_curve <- data.frame(
    tau_m  = tau_grid, logHR = eta, se = se_eta,
    HR     = HR,
    HR_lo  = exp(eta - 1.96 * se_eta),
    HR_hi  = exp(eta + 1.96 * se_eta),
    VE     = VE,
    VE_lo  = 1 - exp(eta + 1.96 * se_eta),
    VE_hi  = 1 - exp(eta - 1.96 * se_eta)
  )

  fit_const <- coxph(
    Surv(start, stop, event_vrs) ~ inmunizado + strata(Group) + cluster(RUN),
    data = df, ties = "efron", robust = TRUE
  )

  list(
    theta_hat      = theta_hat,
    alpha          = alpha,
    V_theta        = V_theta,
    Nmat           = Nmat,
    A_constraints  = A,
    loglik         = final_lgh$ll,
    AIC            = -2 * final_lgh$ll + 2 * n_free,
    LRT_stat       = 2 * (final_lgh$ll - fit_const$loglik[2]),
    LRT_p          = pchisq(2 * (final_lgh$ll - fit_const$loglik[2]),
                             df = n_free - 1, lower.tail = FALSE),
    n_viol_final   = n_viol_final,
    ve_curve       = ve_curve,
    knots_internal = knots_internal,
    tau_star_m     = tau_star_m,
    n_free         = n_free,
    lambda_ineq    = lambda_ineq
  )
}


configs_knots <- list(
  "1k" = list(knots = q(0.50),
              label = "1 internal knot — p50"),
  "2k" = list(knots = q(c(0.33, 0.67)),
              label = "2 internal knots — p33, p67"),
  "3k" = list(knots = q(c(0.10, 0.50, 0.90)),
              label = "3 internal knots — p10, p50, p90"),
  "4k" = list(knots = q(c(0.10, 0.50,0.75, 0.95)),
              label = "4 internal knots — p10, p50, p75, p90")
)

res_diag <- list()
alpha_prev <- NULL

for (nm in names(configs_knots)) {
  cfg <- configs_knots[[nm]]

  res <- fit_cox_nr_full(
    df_raw, cfg$knots, tau_star_m = 15,
    lambda_ineq  = 200,
    alpha_init   = alpha_prev,   # pasar alpha del modelo anterior
    verbose      = TRUE
  )
  results_knots[[nm]] <- res
  alpha_prev <- c(res$alpha, rep(0, 1))  # el siguiente tiene 1 param más
}


# Figura: facetas con IC
bind_rows(lapply(names(configs_knots), function(nm) {
  x <- results_knots[[nm]]
  if (is.null(x)) return(NULL)
  x$ve_curve %>% mutate(modelo = configs_knots[[nm]]$label)
})) %>%
  mutate(modelo = factor(modelo, levels = sapply(configs_knots, `[[`, "label"))) %>%
  ggplot(aes(tau_m, VE, color = modelo, fill = modelo)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
  geom_ribbon(aes(ymin = VE_lo, ymax = VE_hi),
              alpha = 0.12, linewidth = 0) +
  geom_line(linewidth = 1.0) +
  facet_wrap(~ modelo, ncol = 1) +
  scale_color_manual(values = c("#1B4F72", "#C0392B", "#1E8449", "#5b1e84")) +
  scale_fill_manual(values  = c("#1B4F72", "#C0392B", "#1E8449", "#5b1e84")) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(breaks = seq(0, 15, by = 4)) +
  coord_cartesian(xlim = c(0, 15), ylim = c(-0.05, 1.05)) +
  labs(
    x        = "Months since immunization",
    y        = "VE(τ) = 1 − HR(τ)",
    color    = NULL, fill = NULL,
    title    = "Sensitivity to number of internal knots",
    subtitle = "λ = 200  ·  R1: VE(τ₀)=0  ·  R2: VE'(τ₀)=0  ·  R3: VE(τ)≥0"
  ) +
  theme_bw(base_size = 20, base_family = "serif") +
  theme(
    plot.title       = element_text(face = "bold", size = 22),
    plot.subtitle    = element_text(color = "grey35", size = 15),
    legend.position  = "none",
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    strip.text       = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank(),
    aspect.ratio     = 0.45
  )

  