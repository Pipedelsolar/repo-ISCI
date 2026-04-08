library(splines2)
library(survival)
library(dplyr)
library(ggplot2)
library(readr)
library(scales)
library(tidyr)

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

if (!"t_inm" %in% names(df_raw)) {
  ref_date        <- as.Date("2024-04-01")
  df_raw$fechaInm <- as.Date(df_raw$fechaInm)
  df_raw$t_inm    <- as.numeric(df_raw$fechaInm - ref_date)
}
df_raw$t_inm <- as.numeric(df_raw$t_inm)

df_raw <- df_raw %>%
  mutate(
    tau_start_days = ifelse(inmunizado == 1, pmax(0, start - t_inm), 0),
    tau_stop_days  = ifelse(inmunizado == 1, pmax(0, stop  - t_inm), 0),
    tau_eval_days  = ifelse(inmunizado == 1,
                            0.5 * (tau_start_days + tau_stop_days), 0),
    tau_eval_m     = pmax(0, tau_eval_days / 30.4375)
  )

cat("── Resumen datos ──\n")
cat("Personas:          ", n_distinct(df_raw$RUN), "\n")
cat("Eventos totales:   ", sum(df_raw$event_vrs), "\n")
cat("Eventos inmunizados:", sum(df_raw$event_vrs[df_raw$inmunizado == 1]), "\n")

# Percentiles de tau para usar en knots
tau_imm_all <- df_raw %>%
  filter(inmunizado == 1, is.finite(tau_eval_m), tau_eval_m > 0) %>%
  pull(tau_eval_m)

q <- function(p, tau_max = Inf) {
  as.numeric(quantile(tau_imm_all[tau_imm_all < tau_max], p, na.rm = TRUE))
}

# ============================================================
# 1) FUNCIONES CORE: LIKELIHOOD + GRADIENTE + HESSIAN
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

    for (t in sort(unique(t_[ev_ == 1]))) {
      ev_idx   <- which(t_ == t & ev_ == 1)
      risk_idx <- which(s_ < t & t_ >= t)
      if (length(risk_idx) == 0) next

      ee_r <- ee_[risk_idx]
      X_r  <- X_[risk_idx, , drop = FALSE]
      S0   <- sum(ee_r)
      S1   <- colSums(ee_r * X_r)
      S2   <- crossprod(X_r, ee_r * X_r)
      n_ev <- length(ev_idx)

      ll   <- ll   + sum(as.numeric(X_[ev_idx, , drop=FALSE] %*% theta)) -
                     n_ev * log(S0)
      grad <- grad + colSums(X_[ev_idx, , drop=FALSE]) - n_ev * S1 / S0
      H    <- H    - n_ev * (S2/S0 - tcrossprod(S1/S0))
    }
  }
  list(ll = ll, grad = grad, hess = H)
}

# ============================================================
# 2) FAMILIAS FUNCIONALES f(tau) para tau in [0, tau_0]
# ============================================================
# Cada familia devuelve:
#   X      — matriz de diseño [n x p_total] (con beta0, beta1 si aplica, gammas)
#   proj   — proyección para eval posterior (NULL si no aplica)
#   p_names — nombres de columnas

# ── A) Lineal + Spline residualizado ──
make_family_spline <- function(df, tau_0, knots_pct = c(0.50)) {

  tau <- pmin(df$tau_eval_m, tau_0)
  inm <- as.numeric(df$inmunizado == 1)

  # Knots proporcionales a tau_0
  knots_internal <- q(knots_pct, tau_max = tau_0)
  knots_internal <- knots_internal[knots_internal > 0 &
                                   knots_internal < tau_0]

  # Base NCS sin intercepto
  B_raw <- as.matrix(splines2::naturalSpline(
    tau,
    knots          = knots_internal,
    Boundary.knots = c(0, tau_0),
    intercept      = FALSE
  ))
  p_spl <- ncol(B_raw)

  # Residualizar contra {1, tau} en soporte inmunizado activo
  idx_imm <- inm == 1 & df$tau_eval_m < tau_0
  tau_imm <- tau[idx_imm]
  B_imm   <- B_raw[idx_imm, , drop = FALSE]

  proj <- matrix(0, 2, p_spl)
  if (length(tau_imm) > 2) {
    M_imm <- cbind(1, tau_imm)
    XtX   <- crossprod(M_imm)
    for (j in seq_len(p_spl)) {
      proj[, j] <- tryCatch(
        solve(XtX, crossprod(M_imm, B_imm[, j])),
        error = function(e) c(0, 0)
      )
    }
  }

  M_all <- cbind(1, tau)
  B_nl  <- B_raw
  for (j in seq_len(p_spl)) {
    B_nl[, j] <- B_raw[, j] - M_all %*% proj[, j]
  }

  # Para tau >= tau_0: contribución = 0 (activo solo en [0, tau_0))
  active <- as.numeric(inm == 1 & df$tau_eval_m < tau_0)

  X <- cbind(
    beta0 = inm,                        # efecto en tau=0
    beta1 = inm * tau,                  # tendencia lineal
    sweep(B_nl, 1, active, `*`)         # no-linealidad solo en [0, tau_0)
  )
  colnames(X) <- c("beta0", "beta1",
                   paste0("gamma_", seq_len(p_spl)))

  list(X = X, proj = proj, knots_internal = knots_internal,
       family = "spline", tau_0 = tau_0,
       p_names = colnames(X))
}

# ── B) Solo lineal: f(tau) = beta0 + beta1 * tau ──
make_family_linear <- function(df, tau_0, ...) {

  tau <- pmin(df$tau_eval_m, tau_0)
  inm <- as.numeric(df$inmunizado == 1)

  X <- cbind(
    beta0 = inm,
    beta1 = inm * tau
  )
  colnames(X) <- c("beta0", "beta1")

  list(X = X, proj = NULL, knots_internal = NULL,
       family = "linear", tau_0 = tau_0,
       p_names = colnames(X))
}

# ── C) Logística decreciente: f(tau) = beta0 / (1 + exp(beta1*(tau - beta2))) ──
# No lineal en parámetros → se maneja distinto (ver nota abajo)
# Aquí la aproximamos como: beta0 + beta1*tau + beta2*tau^2  (cuadrática)
make_family_quadratic <- function(df, tau_0, ...) {

  tau <- pmin(df$tau_eval_m, tau_0)
  inm <- as.numeric(df$inmunizado == 1)

  X <- cbind(
    beta0 = inm,
    beta1 = inm * tau,
    beta2 = inm * tau^2
  )
  colnames(X) <- c("beta0", "beta1", "beta2")

  list(X = X, proj = NULL, knots_internal = NULL,
       family = "quadratic", tau_0 = tau_0,
       p_names = colnames(X))
}


make_family_cubic <- function(df, tau_0, ...) {

  tau <- pmin(df$tau_eval_m, tau_0)
  inm <- as.numeric(df$inmunizado == 1)

  X <- cbind(
    beta0 = inm,
    # beta1 = inm * tau,
    beta_cub = inm * tau^3
  )
  colnames(X) <- c("beta0", "beta_cub")

  list(X = X, proj = NULL, knots_internal = NULL,
       family = "cubic", tau_0 = tau_0,
       p_names = colnames(X))
}

# ── D) Raíz cuadrada: f(tau) = beta0 + beta1*sqrt(tau) ──
make_family_sqrt <- function(df, tau_0, ...) {

  tau <- pmin(df$tau_eval_m, tau_0)
  inm <- as.numeric(df$inmunizado == 1)

  X <- cbind(
    beta0 = inm,
    beta1 = inm * sqrt(tau)
  )
  colnames(X) <- c("beta0", "beta1_sqrt")

  list(X = X, proj = NULL, knots_internal = NULL,
       family = "sqrt", tau_0 = tau_0,
       p_names = colnames(X))
}

# ── E) Logarítmica: f(tau) = beta0 + beta1*log(tau+1) ──
make_family_log <- function(df, tau_0, ...) {

  tau <- pmin(df$tau_eval_m, tau_0)
  inm <- as.numeric(df$inmunizado == 1)

  X <- cbind(
    beta0 = inm,
    beta1 = inm * log(tau + 1)
  )
  colnames(X) <- c("beta0", "beta1_log")

  list(X = X, proj = NULL, knots_internal = NULL,
       family = "log", tau_0 = tau_0,
       p_names = colnames(X))
}

# ============================================================
# 3) CONSTRUCCIÓN DE RESTRICCIONES R1 (solo valor = 0 en tau_0)
#    NOTA: sin R2 (derivada), por eso f() puede tener quiebre en tau_0
# ============================================================

build_constraint_R1 <- function(X_dm, tau_0, family_obj) {
  # R1: f(tau_0) = 0
  # Para las familias lineales/cuadráticas/sqrt/log:
  # La última fila de X evaluada en tau = tau_0 da el vector de restricción

  p_names <- family_obj$p_names
  p       <- length(p_names)

  if (family_obj$family == "spline") {
    # Evaluar en tau_0 - eps
    proj <- family_obj$proj
    knots_internal <- family_obj$knots_internal
    eps <- 1e-6
    tau_at <- tau_0 - eps

    B_at <- as.matrix(splines2::naturalSpline(
      tau_at,
      knots          = knots_internal,
      Boundary.knots = c(0, tau_0),
      intercept      = FALSE
    ))
    M_at <- cbind(1, tau_at)
    B_nl_at <- B_at
    for (j in seq_len(ncol(B_at))) {
      B_nl_at[, j] <- B_at[, j] - M_at %*% proj[, j]
    }
    # R1: beta0 + beta1*tau_0 + sum gamma_j * B_nl_j(tau_0) = 0
    a_R1 <- c(1, tau_at, as.numeric(B_nl_at))

  } else if (family_obj$family == "linear") {
    a_R1 <- c(1, tau_0)

  } else if (family_obj$family == "quadratic") {
    a_R1 <- c(1, tau_0, tau_0^2)

  } else if (family_obj$family == "sqrt") {
    a_R1 <- c(1, sqrt(tau_0))

  } else if (family_obj$family == "log") {
    a_R1 <- c(1, log(tau_0 + 1))
  }
    else if (family_obj$family == "cubic") {
    a_R1 <- c(1, tau_0^3)
  }

  A <- matrix(a_R1, nrow = 1)
  colnames(A) <- p_names
  rownames(A) <- "R1_value_at_tau0"
  A
}

# ============================================================
# 4) NEWTON-RAPHSON CON NULL-SPACE (R1) + PENALIZACIÓN (R3: VE>=0)
# ============================================================

fit_one_model <- function(df,
                           family_fn,
                           tau_0,
                           family_args  = list(),
                           lambda_ineq  = 200,
                           n_check      = 150,
                           max_iter     = 60,
                           tol          = 1e-9,
                           verbose      = FALSE) {

  # Construir diseño
  dm <- do.call(family_fn, c(list(df = df, tau_0 = tau_0), family_args))
  X  <- dm$X
  p  <- ncol(X)

  # Restricción R1
  A <- build_constraint_R1(X, tau_0, dm)

  # Null-space de A (rank 1 → n_free = p - 1)
  svd_A  <- svd(A, nu = 0, nv = p)
  r      <- sum(svd_A$d > 1e-10 * max(svd_A$d))
  Nmat   <- svd_A$v[, (r+1):p, drop = FALSE]
  n_free <- ncol(Nmat)

  # Covariables reparametrizadas
  Z <- X %*% Nmat
  colnames(Z) <- paste0("z", seq_len(n_free))

  # Grilla para R3 (VE >= 0, log-HR <= 0)
  tau_chk    <- seq(0.01, tau_0 - 0.01, length.out = n_check)
  Xchk_list  <- do.call(family_fn,
                         c(list(df = data.frame(
                           inmunizado  = rep(1L, n_check),
                           tau_eval_m  = tau_chk,
                           start       = rep(0, n_check),
                           stop        = tau_chk,
                           event_vrs   = rep(0L, n_check)
                         ), tau_0 = tau_0), family_args))
  Xchk       <- Xchk_list$X
  Zchk       <- Xchk %*% Nmat

  # Penalización R3
  eval_pen <- function(alpha) {
    if (lambda_ineq == 0)
      return(list(pen=0, grad=numeric(n_free), hess=matrix(0,n_free,n_free)))
    eta_c <- as.numeric(Zchk %*% alpha)
    viol  <- pmax(0, eta_c)
    pen   <- lambda_ineq * sum(viol^2)
    gp    <- 2 * lambda_ineq * as.numeric(t(Zchk) %*% viol)
    av    <- viol > 0
    hp    <- if (any(av)) 2*lambda_ineq*crossprod(Zchk[av,,drop=FALSE]) else
             matrix(0, n_free, n_free)
    list(pen = pen, grad = gp, hess = hp)
  }

  # Warm start
  alpha   <- rep(0, n_free)
  ll_prev <- -Inf
  pp_prev <- Inf

  for (iter in seq_len(max_iter)) {
    lgh <- cox_loglik_grad_hess(
      theta  = alpha, X = Z,
      start  = df$start, stop = df$stop,
      event  = df$event_vrs, strata = df$Group
    )
    pen_obj <- eval_pen(alpha)

    obj   <- lgh$ll - pen_obj$pen
    grad  <- lgh$grad - pen_obj$grad
    hess  <- lgh$hess - pen_obj$hess

    step <- tryCatch(
      solve(hess, grad),
      error = function(e) -grad * 1e-3
    )

    lr        <- 1.0
    alpha_new <- alpha - lr * step

    for (ls in seq_len(30)) {
      ll_t <- tryCatch(
        cox_loglik_grad_hess(
          alpha_new, Z, df$start, df$stop, df$event_vrs, df$Group
        )$ll,
        error = function(e) -Inf
      )
      pt   <- eval_pen(alpha_new)$pen
      if (is.finite(ll_t - pt) && (ll_t - pt) > obj) break
      lr        <- lr * 0.5
      alpha_new <- alpha - lr * step
    }

    delta <- obj - (ll_prev - pp_prev)
    if (verbose) {
      cat(sprintf("  [tau0=%.1f] Iter %2d | ll=%9.3f | pen=%7.3f | lr=%.3f\n",
                  tau_0, iter, lgh$ll, pen_obj$pen, lr))
    }
    if (iter > 1 && abs(delta) < tol) break

    alpha   <- alpha_new
    ll_prev <- lgh$ll
    pp_prev <- pen_obj$pen
  }

  # Estimados finales
  theta_hat        <- as.vector(Nmat %*% alpha)
  names(theta_hat) <- colnames(X)

  final_lgh <- cox_loglik_grad_hess(
    alpha, Z, df$start, df$stop, df$event_vrs, df$Group
  )
  loglik <- final_lgh$ll

  vcov_alpha <- tryCatch(solve(-final_lgh$hess), error = function(e) NULL)
  V_theta    <- if (!is.null(vcov_alpha)) Nmat %*% vcov_alpha %*% t(Nmat) else NULL

  # Predicción
  tau_grid    <- seq(0, tau_0, by = 0.05)
  Xg_list     <- do.call(family_fn,
                          c(list(df = data.frame(
                            inmunizado = rep(1L, length(tau_grid)),
                            tau_eval_m = tau_grid,
                            start      = rep(0, length(tau_grid)),
                            stop       = tau_grid,
                            event_vrs  = rep(0L, length(tau_grid))
                          ), tau_0 = tau_0), family_args))
  Xg  <- Xg_list$X
  eta <- as.numeric(Xg %*% theta_hat)
  eta[tau_grid >= tau_0] <- 0

  if (!is.null(V_theta)) {
    se_eta <- sqrt(pmax(0, rowSums((Xg %*% V_theta) * Xg)))
    # Suavizar SE cerca de tau_0
    near <- tau_grid >= (tau_0 - 0.5)
    w    <- (tau_grid[near] - (tau_0 - 0.5)) / 0.5
    se_eta[near] <- se_eta[near] * (1 - w)^2
    se_eta[tau_grid >= tau_0] <- 0
  } else {
    se_eta <- rep(NA_real_, length(tau_grid))
  }

  HR <- exp(eta); VE <- 1 - HR
  n_viol <- sum(as.numeric(Zchk %*% alpha) > 1e-6)

  list(
    loglik     = loglik,
    AIC        = -2 * loglik + 2 * n_free,
    n_free     = n_free,
    n_viol     = n_viol,
    theta_hat  = theta_hat,
    V_theta    = V_theta,
    tau_0      = tau_0,
    family     = dm$family,
    ve_curve   = data.frame(
      tau_m = tau_grid, logHR = eta, se = se_eta,
      HR = HR, HR_lo = exp(eta - 1.96*se_eta), HR_hi = exp(eta + 1.96*se_eta),
      VE = VE, VE_lo = 1 - exp(eta + 1.96*se_eta), VE_hi = 1 - exp(eta - 1.96*se_eta)
    )
  )
}

# ============================================================
# 5) PROFILE LIKELIHOOD SOBRE tau_0
# ============================================================

profile_tau0 <- function(df,
                          family_fn,
                          family_label,
                          family_args  = list(),
                          tau0_grid    = seq(6, 17, by = 0.5),
                          lambda_ineq  = 200,
                          verbose      = TRUE) {

  results <- vector("list", length(tau0_grid))

  for (i in seq_along(tau0_grid)) {
    tau_0 <- tau0_grid[i]
    if (verbose) cat(sprintf("  [%s] tau_0 = %.1f ...\n", family_label, tau_0))

    results[[i]] <- tryCatch(
      fit_one_model(
        df           = df,
        family_fn    = family_fn,
        tau_0        = tau_0,
        family_args  = family_args,
        lambda_ineq  = lambda_ineq,
        verbose      = FALSE
      ),
      error = function(e) {
        if (verbose) message("    ERROR: ", e$message)
        NULL
      }
    )
  }

  # Profile likelihood
  logliks <- sapply(results, function(x) if (!is.null(x)) x$loglik else NA_real_)
  aics    <- sapply(results, function(x) if (!is.null(x)) x$AIC    else NA_real_)

  best_idx <- which.max(logliks)
  best_tau <- tau0_grid[best_idx]

  # IC al 95% via profile: log-lik > max - 1.92 (chi2_1/2)
  ci_idx <- which(logliks >= max(logliks, na.rm=TRUE) - 1.92)
  ci_lo  <- min(tau0_grid[ci_idx], na.rm=TRUE)
  ci_hi  <- max(tau0_grid[ci_idx], na.rm=TRUE)

  if (verbose) {
    cat(sprintf("\n  [%s] tau_0* = %.1f mo  95%% CI: [%.1f, %.1f]\n",
                family_label, best_tau, ci_lo, ci_hi))
  }

  list(
    family       = family_label,
    tau0_grid    = tau0_grid,
    logliks      = logliks,
    aics         = aics,
    best_tau     = best_tau,
    best_loglik  = logliks[best_idx],
    best_aic     = aics[best_idx],
    best_result  = results[[best_idx]],
    ci_lo        = ci_lo,
    ci_hi        = ci_hi,
    all_results  = results
  )
}

# ============================================================
# 6) CORRER TODAS LAS FAMILIAS
# ============================================================

tau0_grid <- seq(6, 17, by = 0.5)

families <- list(
  list(fn = make_family_linear,    label = "Linear",
       args = list()),
  list(fn = make_family_quadratic, label = "Quadratic",
       args = list()),
  list(fn = make_family_sqrt,      label = "Sqrt(τ)",
       args = list()),
  list(fn = make_family_log,       label = "Log(τ+1)",
       args = list()),
  list(fn = make_family_spline,    label = "Spline (1 knot)",
       args = list(knots_pct = 0.50)),
  list(fn = make_family_spline,    label = "Spline (2 knots)",
       args = list(knots_pct = c(0.33, 0.67)))
)

profiles <- list()

for (fam in families) {
  cat(sprintf("\n══ Familia: %s ══\n", fam$label))
  profiles[[fam$label]] <- profile_tau0(
    df           = df_raw,
    family_fn    = fam$fn,
    family_label = fam$label,
    family_args  = fam$args,
    tau0_grid    = tau0_grid,
    lambda_ineq  = 200,
    verbose      = TRUE
  )
}

leo_cubic =  list(fn = make_family_cubic, label = "Cubic", args = list())

cat(sprintf("\n══ Familia: %s ══\n", leo_cubic$label))
profiles[[leo_cubic$label]] <- profile_tau0(
    df           = df_raw,
    family_fn    = leo_cubic$fn,
    family_label = leo_cubic$label,
    family_args  = leo_cubic$args,
    tau0_grid    = tau0_grid,
    lambda_ineq  = 200,
    verbose      = TRUE
    )

# ============================================================
# 7) TABLA RESUMEN: mejor tau_0 por familia
# ============================================================

tab_summary <- bind_rows(lapply(profiles, function(p) {
  data.frame(
    family      = p$family,
    best_tau_0  = round(p$best_tau, 1),
    ci_lo       = round(p$ci_lo, 1),
    ci_hi       = round(p$ci_hi, 1),
    loglik      = round(p$best_loglik, 3),
    AIC         = round(p$best_aic, 2),
    n_free      = p$best_result$n_free,
    n_viol      = p$best_result$n_viol
  )
})) %>%
  mutate(delta_AIC = round(AIC - min(AIC), 2)) %>%
  arrange(AIC)

cat("\n══ Resumen: mejor modelo por familia ══\n")
print(tab_summary, row.names = FALSE)

# ============================================================
# 8) FIGURA A: Profile likelihood por familia
# ============================================================

df_profile <- bind_rows(lapply(profiles, function(p) {
  data.frame(
    family  = p$family,
    tau_0   = p$tau0_grid,
    loglik  = p$logliks,
    delta_ll = p$logliks - p$best_loglik   # relativo al máximo de cada familia
  )
}))



p_profile <- df_profile %>%
  ggplot(aes(tau_0, delta_ll, color = family, group = family)) +
  geom_hline(yintercept = -1.92, linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
  geom_line(linewidth = 0.9) +
  geom_point(data = tab_summary,
             aes(x = best_tau_0, y = 0),
             size = 3, shape = 21, fill = "white", stroke = 1.5) +
  annotate("text", x = min(tau0_grid), y = -1.75,
           label = "95% CI threshold (Δll = -1.92)",
           hjust = 0, size = 5, color = "grey40", family = "serif") +
  scale_x_continuous(breaks = seq(6, 18, by = 1)) +
  scale_y_continuous(breaks = seq(-6, 0, by = 1)) +
  coord_cartesian(ylim = c(-6, 0.3)) +
  labs(
    x        = "τ₀ (months)",
    y        = "Δ log-likelihood (relative to maximum)",
    color    = "Family",
    title    = "Profile likelihood over τ₀",
    subtitle = "Each curve: max over shape parameters given τ₀  ·  λ = 200"
  ) +
  theme_classic(base_size = 20, base_family = "serif") +
  theme(
    plot.title      = element_text(face = "bold", size = 22),
    plot.subtitle   = element_text(color = "grey35", size = 15),
    axis.title      = element_text(size = 18),
    axis.text       = element_text(size = 16),
    legend.position = "bottom",
    legend.text     = element_text(size = 13),
    aspect.ratio    = 0.55
  )

print(p_profile)
ggsave("profile_tau0.pdf", p_profile, width = 13, height = 7, device = cairo_pdf)

# ============================================================
# 9) FIGURA B: Curvas VE en el tau_0 óptimo por familia
# ============================================================

df_best_curves <- bind_rows(lapply(profiles, function(p) {
  if (is.null(p$best_result)) return(NULL)
  p$best_result$ve_curve %>%
    mutate(
      family    = p$family,
      best_tau0 = p$best_tau
    ) %>%
    filter(family != "Sqrt(τ)")
}))

# Extender curva a tau_max con VE = 0 después de tau_0
tau_max <- max(tau0_grid) + 1
df_best_curves_ext <- bind_rows(
  df_best_curves,
  df_best_curves %>%
    group_by(family, best_tau0) %>%
    summarise(.groups = "drop") %>%
    rowwise() %>%
    do(data.frame(
      tau_m     = seq(.$best_tau0, tau_max, by = 0.05),
      logHR     = 0, se = 0,
      HR        = 1, HR_lo = 1, HR_hi = 1,
      VE        = 0, VE_lo = 0, VE_hi = 0,
      family    = .$family,
      best_tau0 = .$best_tau0
    ))
) %>%
  arrange(family, tau_m)

p_best_curves <- df_best_curves_ext %>%
  mutate(
    label = paste0(family, "  (τ₀* = ", best_tau0, " mo)")
  ) %>%
  ggplot(aes(tau_m, VE, color = label, fill = label)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
#   geom_ribbon(aes(ymin = VE_lo, ymax = VE_hi),
#               alpha = 0.10, linewidth = 0) +
  geom_line(linewidth = 1.0) +
  # Marcar tau_0* de cada familia
  geom_vline(
    data = tab_summary,
    aes(xintercept = best_tau_0, color = paste0(family, "  (τ₀* = ", best_tau_0, " mo)")),
    linetype = "dotted", linewidth = 0.6, alpha = 0.7
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     breaks = seq(0, 1, by = 0.25)) +
  scale_x_continuous(breaks = seq(0, tau_max, by = 3)) +
  coord_cartesian(xlim = c(0, tau_max), ylim = c(-0.05, 1.05)) +
  labs(
    x        = "Months since immunization",
    y        = "VE(τ) = 1 − HR(τ)",
    color    = NULL, fill = NULL,
    title    = "Vaccine effectiveness at optimal τ₀",
    subtitle = "Each curve uses its MLE τ₀  ·  R1: VE(τ₀)=0  ·  R3: VE≥0"
  ) +
  theme_classic(base_size = 20, base_family = "serif") +
  theme(
    plot.title      = element_text(face = "bold", size = 22),
    plot.subtitle   = element_text(color = "grey35", size = 15),
    axis.title      = element_text(size = 18),
    axis.text       = element_text(size = 16),
    legend.position = "bottom",
    legend.text     = element_text(size = 11),
    aspect.ratio    = 0.5
  )

print(p_best_curves)
ggsave("ve_optimal_tau0.pdf", p_best_curves,
       width = 13, height = 7, device = cairo_pdf)

# ============================================================
# 10) FIGURA C: Facetas — curva VE en función de tau_0
#     (para la mejor familia según AIC)
# ============================================================

best_family_label <- tab_summary$family[1]
best_profile      <- profiles[[best_family_label]]
best_fam_obj      <- families[[which(sapply(families, `[[`, "label") ==
                                       best_family_label)]]

# Subset de tau_0 para no saturar la figura
tau0_show <- tau0_grid[seq(1, length(tau0_grid), by = 2)]

df_all_tau0 <- bind_rows(lapply(seq_along(best_profile$tau0_grid), function(i) {
  tau_0 <- best_profile$tau0_grid[i]
  if (!tau_0 %in% tau0_show) return(NULL)
  res   <- best_profile$all_results[[i]]
  if (is.null(res)) return(NULL)
  res$ve_curve %>%
    mutate(
      tau_0   = tau_0,
      delta_ll = best_profile$logliks[i] - best_profile$best_loglik
    )
}))

p_tau0_facet <- df_all_tau0 %>%
  mutate(
    tau_label = paste0("τ₀ = ", tau_0, " mo\n(Δll = ",
                       round(delta_ll, 2), ")")
  ) %>%
  ggplot(aes(tau_m, VE)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  geom_ribbon(aes(ymin = VE_lo, ymax = VE_hi),
              fill = "#1B4F72", alpha = 0.12) +
  geom_line(color = "#1B4F72", linewidth = 0.9) +
  facet_wrap(~ tau_label, ncol = 4) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     breaks = c(0, 0.5, 1)) +
  scale_x_continuous(breaks = seq(0, 18, by = 5)) +
  coord_cartesian(ylim = c(-0.05, 1.05)) +
  labs(
    x        = "Months since immunization",
    y        = "VE(τ)",
    title    = paste0("VE curves across τ₀ grid — ", best_family_label),
    subtitle = "Shaded: 95% CI  ·  Δll relative to best τ₀"
  ) +
  theme_bw(base_size = 16, base_family = "serif") +
  theme(
    plot.title       = element_text(face = "bold", size = 18),
    plot.subtitle    = element_text(color = "grey35", size = 13),
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    strip.text       = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank()
  )

print(p_tau0_facet)
ggsave("ve_tau0_facet.pdf", p_tau0_facet,
       width = 14, height = 10, device = cairo_pdf)



# La forma más directa
best_spline2 <- profiles[["Spline (2 knots)"]]$best_result
tau_0_opt    <- best_spline2$tau_0

knots_usados <- q(c(0.33, 0.67), tau_max = tau_0_opt)

cat(sprintf("tau_0 óptimo : %.1f meses\n", tau_0_opt))
cat(sprintf("Knot p33     : %.2f meses\n", knots_usados[1]))
cat(sprintf("Knot p67     : %.2f meses\n", knots_usados[2]))