# ============================================================
# ANÁLISIS DE SENSIBILIDAD: NÚMERO Y UBICACIÓN DE KNOTS
# ============================================================


library(survival)
library(splines)
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

# ============================================================
# 0) DATA
# ============================================================
setwd("C:/Users/ntrig/Desktop/ISCI/Proyectos/Efectividad_Nirse/Code")

df_tv_matched <- read_csv(
  "../Data/df_tv_matched_11_03_2026.csv",
  col_types = cols(
    DIAG9  = col_character(),
    DIAG10 = col_character(),
    DIAG11 = col_character()
  ),
  show_col_types = FALSE
)

df_raw <- df_tv_matched

need <- c("RUN", "start", "stop", "event_vrs", "Group", "inmunizado")
stopifnot(all(need %in% names(df_raw)))

df_raw$RUN        <- as.character(df_raw$RUN)
df_raw$Group      <- as.factor(df_raw$Group)
df_raw$start      <- as.numeric(df_raw$start)
df_raw$stop       <- as.numeric(df_raw$stop)
df_raw$event_vrs  <- as.integer(df_raw$event_vrs)
df_raw$inmunizado <- as.integer(df_raw$inmunizado)

if (any(df_raw$stop < df_raw$start)) {
  message("Advertencia: hay filas con stop <= start. Se eliminarán.")
  df_raw <- df_raw %>% filter(stop >= start)
}

stopifnot(all(df_raw$event_vrs %in% c(0, 1)))
stopifnot(all(df_raw$inmunizado %in% c(0, 1)))

# ============================================================
# 1) TIEMPO DE INMUNIZACIÓN
# ============================================================
if (!"t_inm" %in% names(df_raw)) {
  ref_date <- as.Date("2024-04-01")
  if (!"fechaInm" %in% names(df_raw)) {
    stop("No encuentro ni 't_inm' ni 'fechaInm'.")
  }
  df_raw$fechaInm <- as.Date(df_raw$fechaInm)
  df_raw$t_inm    <- as.numeric(df_raw$fechaInm - ref_date)
}

df_raw$t_inm <- as.numeric(df_raw$t_inm)

if (any(df_raw$inmunizado == 1 & is.na(df_raw$t_inm))) {
  stop("Hay filas con inmunizado=1 pero sin t_inm/fechaInm.")
}

# ============================================================
# DIAGNÓSTICO: estructura de temporadas
# ============================================================
ref_date <- as.Date("2024-04-01")

diag_temporadas <- df_raw %>%
  filter(inmunizado == 1) %>%
  group_by(RUN) %>%
  summarise(t_inm = first(t_inm), .groups = "drop") %>%
  mutate(
    fecha_vac = ref_date + t_inm,
    temporada = case_when(
      fecha_vac < as.Date("2025-03-01") ~ "Temporada 1 (2024)",
      TRUE                               ~ "Temporada 2 (2025)"
    )
  )

cat("── Distribución de vacunaciones por temporada ──\n")
print(table(diag_temporadas$temporada))

# Añadir temporada al df_raw para uso posterior
df_raw <- df_raw %>%
  left_join(
    diag_temporadas %>% select(RUN, temporada),
    by = "RUN"
  ) %>%
  mutate(temporada = ifelse(is.na(temporada), "Control", temporada))

cat("\n── Eventos por temporada ──\n")
df_raw %>%
  group_by(temporada) %>%
  summarise(
    n_personas = n_distinct(RUN),
    n_eventos  = sum(event_vrs),
    .groups    = "drop"
  ) %>%
  print()

# ============================================================
# 2) FUNCIONES AUXILIARES
# ============================================================

collapse_contiguous <- function(data) {
  data %>%
    arrange(RUN, start, stop) %>%
    group_by(RUN) %>%
    mutate(
      t_inm_aux = ifelse(is.na(t_inm), -999999, t_inm),
      new_block =
        row_number() == 1 |
        start != lag(stop) |
        inmunizado != lag(inmunizado) |
        Group != lag(Group) |
        t_inm_aux != lag(t_inm_aux),
      block = cumsum(new_block)
    ) %>%
    group_by(RUN, block) %>%
    summarise(
      Group      = first(Group),
      inmunizado = first(inmunizado),
      t_inm      = first(t_inm),
      start      = first(start),
      stop       = last(stop),
      event_vrs  = as.integer(last(event_vrs)),
      .groups    = "drop"
    ) %>%
    arrange(RUN, start, stop)
}

split_one_interval <- function(start, stop, event, width) {
  if (stop <= start) {
    return(data.frame(
      start     = numeric(0),
      stop      = numeric(0),
      event_vrs = integer(0)
    ))
  }
  cuts <- seq(start, stop, by = width)
  if (length(cuts) == 0 || tail(cuts, 1) < stop) cuts <- c(cuts, stop)
  if (cuts[1] != start) cuts <- c(start, cuts)
  cuts <- unique(cuts)
  out  <- data.frame(start = head(cuts, -1), stop = tail(cuts, -1))
  out  <- out %>% filter(stop > start)
  out$event_vrs <- 0L
  if (nrow(out) > 0) out$event_vrs[nrow(out)] <- as.integer(event)
  out
}

rebuild_intervals <- function(data, split_width_days) {
  data_collapsed <- collapse_contiguous(data)
  split_list     <- vector("list", nrow(data_collapsed))

  for (i in seq_len(nrow(data_collapsed))) {
    row_i <- data_collapsed[i, ]
    if (row_i$inmunizado == 1L) {
      tmp <- split_one_interval(
        start = row_i$start,
        stop  = row_i$stop,
        event = row_i$event_vrs,
        width = split_width_days
      )
      tmp$RUN        <- row_i$RUN
      tmp$Group      <- row_i$Group
      tmp$inmunizado <- row_i$inmunizado
      tmp$t_inm      <- row_i$t_inm
      split_list[[i]] <- tmp[, c("RUN","Group","inmunizado","t_inm","start","stop","event_vrs")]
    } else {
      split_list[[i]] <- row_i[, c("RUN","Group","inmunizado","t_inm","start","stop","event_vrs")]
    }
  }

  bind_rows(split_list) %>%
    arrange(RUN, start, stop) %>%
    filter(stop > start)
}

add_tau_variables <- function(df) {
  df$tau_start_days <- ifelse(df$inmunizado == 1, pmax(0, df$start - df$t_inm), 0)
  df$tau_stop_days  <- ifelse(df$inmunizado == 1, pmax(0, df$stop  - df$t_inm), 0)
  df$tau_eval_days  <- ifelse(
    df$inmunizado == 1,
    0.5 * (df$tau_start_days + df$tau_stop_days), 0
  )
  df$tau_eval_m              <- df$tau_eval_days / 30.4375
  df$tau_eval_m[df$tau_eval_m < 0] <- 0
  df
}

# ============================================================
# HELPERS SPLINE + NULL-SPACE
# ============================================================

null_space <- function(A, tol = 1e-10) {
  s    <- svd(A, nu = 0, nv = ncol(A))
  if (length(s$d) == 0) return(diag(ncol(A)))
  dmax <- max(s$d)
  if (!is.finite(dmax) || dmax == 0) return(diag(ncol(A)))
  r <- sum(s$d > tol * dmax)
  if (r >= ncol(A)) stop("Las restricciones eliminan todo el espacio paramétrico.")
  s$v[, seq.int(r + 1, ncol(A)), drop = FALSE]
}

make_bs_basis <- function(x, spline_df, tau_star_m, knots = NULL) {
  if (is.null(knots)) {
    B <- bs(x, df = spline_df, degree = 3,
            Boundary.knots = c(0, tau_star_m), intercept = FALSE)
    list(B = B, knots = attr(B, "knots"))
  } else {
    B <- bs(x, knots = knots, degree = 3,
            Boundary.knots = c(0, tau_star_m), intercept = FALSE)
    list(B = B, knots = knots)
  }
}

residualize_basis_against_linear <- function(B_support, tau_support, B_all, tau_all) {
  M_support <- cbind(1, tau_support)
  M_all     <- cbind(1, tau_all)
  B_nl_all  <- B_all
  proj_coef <- matrix(NA_real_, nrow = ncol(M_support), ncol = ncol(B_support))
  XtX       <- crossprod(M_support)
  for (j in seq_len(ncol(B_support))) {
    coef_j         <- solve(XtX, crossprod(M_support, B_support[, j]))
    proj_coef[, j] <- coef_j
    B_nl_all[, j]  <- B_all[, j] - M_all %*% coef_j
  }
  list(B_nl_all = B_nl_all, proj_coef = proj_coef)
}

eval_nonlinear_basis <- function(x, knots, spline_df, tau_star_m, proj_coef) {
  x_clip <- pmin(pmax(x, 0), tau_star_m)
  B_raw  <- make_bs_basis(x_clip, spline_df, tau_star_m, knots)$B
  if (ncol(B_raw) != ncol(proj_coef)) {
    stop(paste0("bs() devolvió ", ncol(B_raw), " columnas pero proj_coef tiene ",
                ncol(proj_coef), "."))
  }
  M    <- cbind(1, x_clip)
  B_nl <- B_raw
  for (j in seq_len(ncol(B_raw))) {
    B_nl[, j] <- B_raw[, j] - M %*% proj_coef[, j]
  }
  B_nl
}

# ============================================================
# FIT_WANING_MODEL — sin cambios estructurales
# CORRECCIÓN: tau_grid_m acotado a tau_star_m
# ============================================================

fit_waning_model <- function(df_prepared,
                              split_width_days,
                              spline_df,
                              tau_grid_by = 0.05,
                              tau_star_m  = 15) {

  df            <- df_prepared
  tau_eval_clip <- pmin(df$tau_eval_m, tau_star_m)
  active_idx    <- df$inmunizado == 1 & df$tau_eval_m < tau_star_m
  tau_imm       <- tau_eval_clip[active_idx]
  tau_imm       <- tau_imm[is.finite(tau_imm)]

  if (length(unique(tau_imm)) <= spline_df) {
    stop(paste0("No hay suficiente soporte único en tau para spline_df=", spline_df,
                " con malla ", split_width_days, " días."))
  }

  # Base cúbica en soporte activo
  bs_support_obj <- make_bs_basis(tau_imm, spline_df, tau_star_m, knots = NULL)
  knots          <- bs_support_obj$knots
  B_support      <- bs_support_obj$B

  # Base para todas las filas
  B_all <- make_bs_basis(tau_eval_clip, spline_df, tau_star_m, knots)$B

  # Residualizar
  rez       <- residualize_basis_against_linear(B_support, tau_imm, B_all, tau_eval_clip)
  B_nl_all  <- rez$B_nl_all
  proj_coef <- rez$proj_coef

  # Diseño original
  active <- as.numeric(df$inmunizado == 1 & df$tau_eval_m < tau_star_m)
  X_orig <- cbind(
    beta0 = active,
    beta1 = active * df$tau_eval_m,
    active * B_nl_all
  )
  colnames(X_orig) <- c("beta0", "beta1", paste0("gamma", seq_len(ncol(B_nl_all))))

  # Restricciones en tau_star
  tau_at_star <- tau_star_m - 1e-6
  eps         <- 1e-6

  B_at_star   <- eval_nonlinear_basis(tau_at_star,       knots, spline_df, tau_star_m, proj_coef)
  B_at_star_p <- eval_nonlinear_basis(tau_at_star + eps, knots, spline_df, tau_star_m, proj_coef)
  B_at_star_m <- eval_nonlinear_basis(tau_at_star - eps, knots, spline_df, tau_star_m, proj_coef)
  dB_at_star  <- (B_at_star_p - B_at_star_m) / (2 * eps)

  A <- rbind(
    c(1, tau_at_star, as.numeric(B_at_star)),
    c(0, 1,           as.numeric(dB_at_star))
  )
  colnames(A) <- colnames(X_orig)
  rownames(A) <- c("value_at_tau0", "slope_at_tau0")

  # Null-space
  Nmat <- null_space(A)
  Z    <- X_orig %*% Nmat
  colnames(Z) <- paste0("z", seq_len(ncol(Z)))
  df_fit <- cbind(df, as.data.frame(Z))

  rhs_ns <- paste(c(colnames(Z), "strata(Group)", "cluster(RUN)"), collapse = " + ")
  f_ns   <- as.formula(paste0("Surv(start, stop, event_vrs) ~ ", rhs_ns))

  fit_ns <- coxph(f_ns, data = df_fit, ties = "efron",
                  robust = TRUE, x = TRUE, model = TRUE)

  fit_const <- coxph(
    Surv(start, stop, event_vrs) ~ inmunizado + strata(Group) + cluster(RUN),
    data = df, ties = "efron", robust = TRUE, x = TRUE, model = TRUE
  )

  # Recuperar theta
  alpha_full <- coef(fit_ns)
  keep_alpha <- !is.na(alpha_full)
  alpha      <- alpha_full[keep_alpha]
  V_alpha    <- vcov(fit_ns)[keep_alpha, keep_alpha, drop = FALSE]
  N_keep     <- Nmat[, keep_alpha, drop = FALSE]

  theta_hat <- as.vector(N_keep %*% alpha)
  names(theta_hat) <- colnames(X_orig)
  V_theta <- N_keep %*% V_alpha %*% t(N_keep)
  rownames(V_theta) <- colnames(V_theta) <- colnames(X_orig)

  # ── CORRECCIÓN CLAVE: grilla acotada a tau_star_m ──
  # Antes llegaba hasta max(tau_eval_m) ~ 18 meses, lo que causaba
  # extrapolación fuera del boundary knot y curvas inestables
  tau_grid_m <- seq(0, tau_star_m, by = tau_grid_by)

  B_grid_nl <- eval_nonlinear_basis(tau_grid_m, knots, spline_df, tau_star_m, proj_coef)
  active_grid <- as.numeric(tau_grid_m < tau_star_m)

  Xg_orig <- cbind(
    beta0 = active_grid,
    beta1 = active_grid * tau_grid_m,
    active_grid * B_grid_nl
  )
  colnames(Xg_orig) <- colnames(X_orig)

  eta    <- as.vector(Xg_orig %*% theta_hat)
  se_eta <- sqrt(rowSums((Xg_orig %*% V_theta) * Xg_orig))

  HR    <- exp(eta)
  HR_lo <- exp(eta - 1.96 * se_eta)
  HR_hi <- exp(eta + 1.96 * se_eta)
  VE    <- 1 - HR
  VE_lo <- 1 - HR_hi
  VE_hi <- 1 - HR_lo

  ve_curve <- data.frame(
    split_width_days = split_width_days,
    spline_df        = spline_df,
    tau_0            = tau_star_m,
    tau_m            = tau_grid_m,
    logHR            = eta,
    se               = se_eta,
    HR = HR, HR_lo = HR_lo, HR_hi = HR_hi,
    VE = VE, VE_lo = VE_lo, VE_hi = VE_hi
  )

  # Soporte
  support           <- subset(df, inmunizado == 1)
  support$month_bin <- floor(support$tau_eval_m)
  support_summary   <- aggregate(
    cbind(person_time_days = stop - start, events = event_vrs) ~ month_bin,
    data = support, FUN = sum
  )
  n_rows_tab <- table(support$month_bin)
  support_summary$n_rows          <- as.numeric(n_rows_tab[as.character(support_summary$month_bin)])
  support_summary$split_width_days <- split_width_days
  support_summary$spline_df        <- spline_df
  support_summary$tau_0            <- tau_star_m

  # Summary
  model_summary <- data.frame(
    split_width_days = split_width_days,
    spline_df        = spline_df,
    tau_0            = tau_star_m,
    n_rows           = nrow(df),
    n_subjects       = dplyr::n_distinct(df$RUN),
    n_events         = sum(df$event_vrs),
    loglik_ns        = fit_ns$loglik[2],
    loglik_const     = fit_const$loglik[2],
    LRT_stat         = 2 * (fit_ns$loglik[2] - fit_const$loglik[2]),
    df_diff          = sum(!is.na(alpha_full)) - sum(!is.na(coef(fit_const))),
    p_value          = pchisq(
      2 * (fit_ns$loglik[2] - fit_const$loglik[2]),
      df         = sum(!is.na(alpha_full)) - sum(!is.na(coef(fit_const))),
      lower.tail = FALSE
    )
  )

  list(
    df_model        = df_fit,
    fit_ns          = fit_ns,
    fit_const       = fit_const,
    ve_curve        = ve_curve,
    support_summary = support_summary,
    model_summary   = model_summary,
    theta_hat       = theta_hat,
    V_theta         = V_theta,
    A_constraints   = A,
    N_null          = Nmat,
    knots           = knots,
    proj_coef       = proj_coef
  )
}

# ============================================================
# 3) PREPARAR MALLAS
# ============================================================

split_grid  <- c(14)
spline_grid <- c(3, 4, 5, 6)
tau_star_gird <- c(13, 14, 15)

prepared_by_width <- list()
for (w in split_grid) {
  message("Preparando malla de ", w, " días...")
  tmp <- rebuild_intervals(df_raw, split_width_days = w)
  tmp <- add_tau_variables(tmp)
  prepared_by_width[[paste0("w", w)]] <- tmp
}

# ============================================================
# 4) FIT DE LOS MODELOS
# ============================================================

results <- list()

for (w in split_grid) {
  for (sdf in spline_grid) {
    for (tau_0 in tau_star_gird) {
      key <- paste0("w", w, "_df", sdf, "_tau_", tau_0)
      message("Ajustando modelo: malla=", w, " días ; spline_df=", sdf, " tau=", tau_0)

      results[[key]] <- tryCatch(
        fit_waning_model(
          df_prepared      = prepared_by_width[[paste0("w", w)]],
          split_width_days = w,
          spline_df        = sdf,
          tau_grid_by      = 0.05,
          tau_star_m       = tau_0
        ),
        error = function(e) {
          message("  ERROR: ", conditionMessage(e))
          NULL
        }
      )
    }
  }
}

# ============================================================
# 5) TABLAS CONSOLIDADAS
# ============================================================

all_curves <- imap_dfr(results, function(res, key) {
  if (is.null(res)) return(NULL)
  m <- str_match(key, "^w(\\d+)_df(\\d+)_tau_(\\d+)$")
  res$ve_curve %>%
    mutate(
      split_width_days = as.integer(m[2]),
      spline_df        = as.integer(m[3]),
      tau_0            = as.integer(m[4]),
      model_key        = key
    )
})

all_support       <- bind_rows(lapply(results, function(x) if (!is.null(x)) x$support_summary))
all_model_summary <- bind_rows(lapply(results, function(x) if (!is.null(x)) x$model_summary))

# ── Extraer tau_imm del dataset preparado ──
df_w14      <- prepared_by_width[["w14"]]
tau_imm_w14 <- with(df_w14,
  pmin(tau_eval_m, 15)[inmunizado == 1 & tau_eval_m < 15]
)
tau_imm_w14 <- tau_imm_w14[is.finite(tau_imm_w14) & tau_imm_w14 > 0]

cat("── Percentiles de tau_imm ──\n")
print(round(quantile(tau_imm_w14, probs = seq(0, 1, by = 0.05)), 2))

# ============================================================
# DEFINIR CONFIGURACIONES DE KNOTS
# Convención papers: "k knots" = k knots INTERNOS
# En R: spline_df = k_internos + degree = k_internos + 3
# ============================================================

q <- function(p) as.numeric(quantile(tau_imm_w14, probs = p))

# Redefinir knot_configs eliminando 4 y 5 knots
knot_configs <- list(

  # 0 knots — un solo polinomio cúbico
  "0k_global" = list(knots = NULL, spline_df = 3,
                     label = "0 knots"),

  # 1 knot interno
  "1k_p33"  = list(knots = q(0.33), spline_df = 4, label = "1 knot — p33"),
  "1k_p50"  = list(knots = q(0.50), spline_df = 4, label = "1 knot — p50"),
  "1k_p67"  = list(knots = q(0.67), spline_df = 4, label = "1 knot — p67"),
  "1k_p90"  = list(knots = q(0.90), spline_df = 4, label = "1 knot — p90"),
  "1k_p95"  = list(knots = q(0.95), spline_df = 4, label = "1 knot — p95"),

  # 2 knots internos
  "2k_p25_p75" = list(knots = q(c(0.25, 0.75)), spline_df = 5,
                       label = "2 knots — p25, p75"),
  "2k_p33_p67" = list(knots = q(c(0.33, 0.67)), spline_df = 5,
                       label = "2 knots — p33, p67"),
  "2k_p25_p95" = list(knots = q(c(0.25, 0.95)), spline_df = 5,
                       label = "2 knots — p25, p95"),
  "2k_p50_p95" = list(knots = q(c(0.50, 0.95)), spline_df = 5,
                       label = "2 knots — p50, p95"),

  # 3 knots internos — límite razonable con 160 eventos
  "3k_p10_p50_p90" = list(knots = q(c(0.10, 0.50, 0.90)), spline_df = 6,
                            label = "3 knots — p10, p50, p90"),
  "3k_p05_p50_p95" = list(knots = q(c(0.05, 0.50, 0.95)), spline_df = 6,
                            label = "3 knots — p05, p50, p95"),
  "3k_p25_p50_p75" = list(knots = q(c(0.25, 0.50, 0.75)), spline_df = 6,
                            label = "3 knots — p25, p50, p75")
)

cat(sprintf("\nTotal de configuraciones a ajustar: %d\n", length(knot_configs)))

# ============================================================
# FUNCIÓN DE FIT CON KNOTS MANUALES
# ============================================================

fit_waning_knots <- function(df_prepared,
                              split_width_days,
                              spline_df,
                              knots_manual = NULL,
                              label        = "",
                              tau_grid_by  = 0.05,
                              tau_star_m   = 15) {

  df            <- df_prepared
  tau_eval_clip <- pmin(df$tau_eval_m, tau_star_m)
  active_idx    <- df$inmunizado == 1 & df$tau_eval_m < tau_star_m
  tau_imm       <- tau_eval_clip[active_idx]
  tau_imm       <- tau_imm[is.finite(tau_imm)]

  # ── Base spline con knots especificados ──
  bs_support_obj <- make_bs_basis(tau_imm, spline_df, tau_star_m,
                                   knots = knots_manual)
  knots     <- bs_support_obj$knots
  B_support <- bs_support_obj$B
  B_all     <- make_bs_basis(tau_eval_clip, spline_df, tau_star_m, knots)$B

  rez       <- residualize_basis_against_linear(B_support, tau_imm, B_all, tau_eval_clip)
  B_nl_all  <- rez$B_nl_all
  proj_coef <- rez$proj_coef

  active <- as.numeric(df$inmunizado == 1 & df$tau_eval_m < tau_star_m)
  X_orig <- cbind(
    beta0 = active,
    beta1 = active * df$tau_eval_m,
    active * B_nl_all
  )
  colnames(X_orig) <- c("beta0", "beta1",
                         paste0("gamma", seq_len(ncol(B_nl_all))))

  # Restricciones
  tau_at_star <- tau_star_m - 1e-6
  eps         <- 1e-6
  B_at_star   <- eval_nonlinear_basis(tau_at_star,       knots, spline_df, tau_star_m, proj_coef)
  B_at_star_p <- eval_nonlinear_basis(tau_at_star + eps, knots, spline_df, tau_star_m, proj_coef)
  B_at_star_m <- eval_nonlinear_basis(tau_at_star - eps, knots, spline_df, tau_star_m, proj_coef)
  dB_at_star  <- (B_at_star_p - B_at_star_m) / (2 * eps)

  A <- rbind(
    c(1, tau_at_star, as.numeric(B_at_star)),
    c(0, 1,           as.numeric(dB_at_star))
  )
  colnames(A) <- colnames(X_orig)
  rownames(A) <- c("value_at_tau0", "slope_at_tau0")

  Nmat <- null_space(A)
  Z    <- X_orig %*% Nmat
  colnames(Z) <- paste0("z", seq_len(ncol(Z)))
  df_fit <- cbind(df, as.data.frame(Z))

  rhs_ns <- paste(c(colnames(Z), "strata(Group)", "cluster(RUN)"), collapse = " + ")
  f_ns   <- as.formula(paste0("Surv(start, stop, event_vrs) ~ ", rhs_ns))

  fit_ns <- coxph(f_ns, data = df_fit, ties = "efron",
                  robust = TRUE, x = TRUE, model = TRUE)

  fit_const <- coxph(
    Surv(start, stop, event_vrs) ~ inmunizado + strata(Group) + cluster(RUN),
    data = df, ties = "efron", robust = TRUE
  )

  alpha_full <- coef(fit_ns)
  keep_alpha <- !is.na(alpha_full)
  alpha      <- alpha_full[keep_alpha]
  V_alpha    <- vcov(fit_ns)[keep_alpha, keep_alpha, drop = FALSE]
  N_keep     <- Nmat[, keep_alpha, drop = FALSE]

  theta_hat        <- as.vector(N_keep %*% alpha)
  names(theta_hat) <- colnames(X_orig)
  V_theta          <- N_keep %*% V_alpha %*% t(N_keep)
  rownames(V_theta) <- colnames(V_theta) <- colnames(X_orig)

  tau_grid_m  <- seq(0, tau_star_m, by = tau_grid_by)
  B_grid_nl   <- eval_nonlinear_basis(tau_grid_m, knots, spline_df, tau_star_m, proj_coef)
  active_grid <- as.numeric(tau_grid_m < tau_star_m)

  Xg_orig <- cbind(
    beta0 = active_grid,
    beta1 = active_grid * tau_grid_m,
    active_grid * B_grid_nl
  )
  colnames(Xg_orig) <- colnames(X_orig)

  eta    <- as.vector(Xg_orig %*% theta_hat)
  se_eta <- sqrt(rowSums((Xg_orig %*% V_theta) * Xg_orig))
  HR     <- exp(eta)
  VE     <- 1 - HR

  n_free <- sum(keep_alpha)

  ve_curve <- data.frame(
    config    = label,
    n_knots   = length(knots),
    spline_df = spline_df,
    tau_m     = tau_grid_m,
    logHR     = eta,   se = se_eta,
    HR  = HR,
    HR_lo = exp(eta - 1.96 * se_eta),
    HR_hi = exp(eta + 1.96 * se_eta),
    VE    = VE,
    VE_lo = 1 - exp(eta + 1.96 * se_eta),
    VE_hi = 1 - exp(eta - 1.96 * se_eta)
  )

  list(
    fit_ns    = fit_ns,
    ve_curve  = ve_curve,
    knots     = knots,
    n_free    = n_free,
    loglik    = fit_ns$loglik[2],
    loglik_const = fit_const$loglik[2],
    AIC       = -2 * fit_ns$loglik[2] + 2 * n_free,
    LRT_stat  = 2 * (fit_ns$loglik[2] - fit_const$loglik[2]),
    LRT_p     = pchisq(2 * (fit_ns$loglik[2] - fit_const$loglik[2]),
                       df = n_free - 1, lower.tail = FALSE)
  )
}

# ============================================================
# CORRER TODAS LAS CONFIGURACIONES
# ============================================================

results_knot_grid <- list()

for (cfg_name in names(knot_configs)) {
  cfg <- knot_configs[[cfg_name]]
  message(sprintf("Ajustando: %s", cfg$label))

  results_knot_grid[[cfg_name]] <- tryCatch(
    fit_waning_knots(
      df_prepared      = prepared_by_width[["w14"]],
      split_width_days = 14,
      spline_df        = cfg$spline_df,
      knots_manual     = cfg$knots,
      label            = cfg$label,
      tau_star_m       = 15
    ),
    error = function(e) {
      message("  ERROR: ", e$message)
      NULL
    }
  )
}

# ============================================================
# TABLA COMPARATIVA
# ============================================================

tab_knot_grid <- bind_rows(lapply(names(results_knot_grid), function(k) {
  x <- results_knot_grid[[k]]
  if (is.null(x)) return(NULL)
  cfg <- knot_configs[[k]]
  data.frame(
    config      = cfg$label,
    n_knots_int = length(x$knots),
    spline_df   = cfg$spline_df,
    knot_pos    = ifelse(length(x$knots) == 0, "—",
                         paste(round(x$knots, 1), collapse = ", ")),
    n_free      = x$n_free,
    loglik      = round(x$loglik, 2),
    AIC         = round(x$AIC, 2),
    LRT_p       = ifelse(x$LRT_p < 0.001, "<0.001",
                         sprintf("%.3f", x$LRT_p))
  )
})) %>%
  mutate(delta_AIC = round(AIC - min(AIC), 2)) %>%
  arrange(AIC)

cat("\n══ Comparación completa de configuraciones de knots ══\n")
print(tab_knot_grid, row.names = FALSE)

# Guardar tabla HTML
library(kableExtra)
tab_knot_grid %>%
  kable("html",
        caption = "Knot placement sensitivity — spline df=3 to 8, τ₀=15 mo",
        col.names = c("Configuration", "Int. knots", "df", "Knot positions (mo)",
                      "Free params", "Log-lik", "AIC", "LRT p", "ΔAIC")) %>%
  kable_styling(full_width = FALSE,
                bootstrap_options = c("striped","hover","condensed"),
                font_size = 13) %>%
  row_spec(which.min(tab_knot_grid$AIC),
           bold = TRUE, background = "#D6EAF8") %>%
  row_spec(0, bold = TRUE, background = "#2C3E50", color = "white") %>%
  as.character() %>%
  {writeLines(paste0(
    '<!DOCTYPE html><html><head><meta charset="utf-8">',
    '<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css">',
    '<style>body{padding:30px;font-family:Georgia,serif;}</style>',
    '</head><body>', ., '</body></html>'
  ), "knot_sensitivity_table.html")}
cat("✓ Tabla guardada en knot_sensitivity_table.html\n")

# ============================================================
# FIGURAS: una por número de knots internos
# ============================================================

all_knot_curves <- bind_rows(
  lapply(results_knot_grid, function(x) if (!is.null(x)) x$ve_curve)
)

# Paleta por número de knots
knot_colors <- c(
  "0" = "grey50",
  "1" = "#1B4F72",
  "2" = "#C0392B",
  "3" = "#1E8449",
  "4" = "#7D3C98",
  "5" = "#E67E22"
)

# ── Figura A: overlay por n_knots (sin IC) ──
p_knot_overlay <- all_knot_curves %>%
  mutate(
    n_knots_f = factor(n_knots,
                       levels = 0:5,
                       labels = paste0(0:5, " internal knots"))
  ) %>%
  ggplot(aes(x = tau_m, y = VE, color = n_knots_f, group = config)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.5) +
  geom_line(linewidth = 0.8, alpha = 0.75) +
  scale_color_manual(values = setNames(knot_colors, paste0(0:5, " internal knots"))) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(-0.25, 1, by = 0.25)) +
  scale_x_continuous(breaks = seq(0, 15, by = 3)) +
  coord_cartesian(xlim = c(0, 15), ylim = c(-0.15, 1.05)) +
  labs(
    x        = "Months since immunization",
    y        = "VE(τ) = 1 − HR(τ)",
    color    = "Knot count",
    title    = "Sensitivity to knot number and placement",
    subtitle = "Each line = one knot configuration  ·  τ₀ = 15 mo  ·  split = 14 d"
  ) +
  theme_classic(base_size = 20, base_family = "serif") +
  theme(
    plot.title      = element_text(face = "bold", size = 22),
    plot.subtitle   = element_text(color = "grey35", size = 16),
    axis.title      = element_text(size = 18),
    axis.text       = element_text(size = 16),
    legend.position = "bottom",
    legend.text     = element_text(size = 13),
    aspect.ratio    = 0.55
  )

print(p_knot_overlay)
ggsave("knot_sensitivity_overlay.pdf", p_knot_overlay,
       width = 12, height = 7, device = cairo_pdf)

p_knot_facet <- all_knot_curves %>%
  mutate(
    n_knots_f = factor(
      paste0(n_knots, " internal knots"),
      levels = c("0 internal knots", "1 internal knots",
                 "2 internal knots", "3 internal knots")
    )
  ) %>%
  ggplot(aes(x = tau_m, y = VE,
             color = config, fill = config, group = config)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  geom_ribbon(aes(ymin = VE_lo, ymax = VE_hi),
              alpha = 0.10, linewidth = 0) +
  geom_line(linewidth = 0.85) +
  facet_wrap(~ n_knots_f, ncol = 2) +   # 2 columnas — más legible
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = c(0, 0.5, 1)) +
  scale_x_continuous(breaks = seq(0, 15, by = 5)) +
  coord_cartesian(xlim = c(0, 15), ylim = c(-0.15, 1.05)) +
  labs(
    x        = "Months since immunization",
    y        = "VE(τ) = 1 − HR(τ)",
    title    = "Knot placement sensitivity by knot count",
    subtitle = sprintf("160 events  ·  Max 3 internal knots (EPV ≥ 40)  ·  τ₀ = 15 mo")
  ) +
  theme_bw(base_size = 18, base_family = "serif") +
  theme(
    plot.title       = element_text(face = "bold", size = 20),
    plot.subtitle    = element_text(color = "grey35", size = 14),
    axis.title       = element_text(size = 16),
    axis.text        = element_text(size = 14),
    legend.position  = "none",
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    strip.text       = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank()
  )

print(p_knot_facet)
ggsave("knot_sensitivity_facet.pdf", p_knot_facet,
       width = 14, height = 9, device = cairo_pdf)

# ── Figura C: AIC por configuración ──
p_aic <- tab_knot_grid %>%
  mutate(
    config    = factor(config, levels = rev(config)),
    n_knots_f = factor(paste0(n_knots_int, " knots"),
                       levels = paste0(0:5, " knots"))
  ) %>%
  ggplot(aes(x = delta_AIC, y = config, fill = n_knots_f)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 2, linetype = "dashed",
             color = "#C0392B", linewidth = 0.7) +
  annotate("text", x = 2.3, y = 1,
           label = "ΔAIC = 2", hjust = 0, size = 6,
           color = "#C0392B", family = "serif") +
  scale_fill_manual(values = knot_colors,
                    labels = paste0(0:5, " internal knots")) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(
    x     = "ΔAIC (vs best model)",
    y     = NULL,
    fill  = "Knot count",
    title = "Model comparison — knot configurations",
    subtitle = "Models with ΔAIC < 2 are equivalent"
  ) +
  theme_classic(base_size = 18, base_family = "serif") +
  theme(
    plot.title      = element_text(face = "bold", size = 20),
    plot.subtitle   = element_text(color = "grey35", size = 14),
    axis.title      = element_text(size = 16),
    axis.text       = element_text(size = 13),
    legend.position = "bottom",
    legend.text     = element_text(size = 12)
  )

print(p_aic)
ggsave("knot_aic_comparison.pdf", p_aic,
       width = 12, height = 9, device = cairo_pdf)