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

# Columnas mínimas
need <- c("RUN", "start", "stop", "event_vrs", "Group", "inmunizado")
stopifnot(all(need %in% names(df_raw)))

# Tipos
df_raw$RUN        <- as.character(df_raw$RUN)
df_raw$Group      <- as.factor(df_raw$Group)
df_raw$start      <- as.numeric(df_raw$start)
df_raw$stop       <- as.numeric(df_raw$stop)
df_raw$event_vrs  <- as.integer(df_raw$event_vrs)
df_raw$inmunizado <- as.integer(df_raw$inmunizado)

# OJO: en Cox counting-process idealmente stop > start.
# Si tienes igualdad exacta, mejor revisar o filtrar.
if (any(df_raw$stop < df_raw$start)) {
  message("Advertencia: hay filas con stop <= start. Se eliminarán.")
  df_raw <- df_raw %>% filter(stop >= start)
}

stopifnot(all(df_raw$event_vrs %in% c(0, 1)))
stopifnot(all(df_raw$inmunizado %in% c(0, 1)))

# ============================================================
# 1) TIEMPO DE INMUNIZACIÓN EN LA MISMA ESCALA QUE start/stop
# ============================================================
if (!"t_inm" %in% names(df_raw)) {
  ref_date <- as.Date("2024-04-01")

  if (!"fechaInm" %in% names(df_raw)) {
    stop("No encuentro ni 't_inm' ni 'fechaInm'. Necesito una de las dos.")
  }

  df_raw$fechaInm <- as.Date(df_raw$fechaInm)
  df_raw$t_inm    <- as.numeric(df_raw$fechaInm - ref_date)
}

df_raw$t_inm <- as.numeric(df_raw$t_inm)

if (any(df_raw$inmunizado == 1 & is.na(df_raw$t_inm))) {
  stop("Hay filas con inmunizado=1 pero sin t_inm/fechaInm.")
}

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
      start = numeric(0),
      stop = numeric(0),
      event_vrs = integer(0)
    ))
  }

  cuts <- seq(start, stop, by = width)

  if (length(cuts) == 0 || tail(cuts, 1) < stop) {
    cuts <- c(cuts, stop)
  }

  if (cuts[1] != start) {
    cuts <- c(start, cuts)
  }

  cuts <- unique(cuts)

  out <- data.frame(
    start = head(cuts, -1),
    stop  = tail(cuts, -1)
  )

  out <- out %>% filter(stop > start)

  out$event_vrs <- 0L
  if (nrow(out) > 0) {
    out$event_vrs[nrow(out)] <- as.integer(event)
  }

  out
}

rebuild_intervals <- function(data, split_width_days) {
  data_collapsed <- collapse_contiguous(data)

  split_list <- vector("list", nrow(data_collapsed))

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

      split_list[[i]] <- tmp[, c("RUN", "Group", "inmunizado", "t_inm", "start", "stop", "event_vrs")]
    } else {
      split_list[[i]] <- row_i[, c("RUN", "Group", "inmunizado", "t_inm", "start", "stop", "event_vrs")]
    }
  }

  bind_rows(split_list) %>%
    arrange(RUN, start, stop) %>%
    filter(stop > start)
}

add_tau_variables <- function(df) {
  df$tau_start_days <- ifelse(df$inmunizado == 1, pmax(0, df$start - df$t_inm), 0)
  df$tau_stop_days  <- ifelse(df$inmunizado == 1, pmax(0, df$stop  - df$t_inm), 0)

  df$tau_eval_days <- ifelse(
    df$inmunizado == 1,
    0.5 * (df$tau_start_days + df$tau_stop_days),
    0
  )

  df$tau_eval_m <- df$tau_eval_days / 30.4375
  df$tau_eval_m[df$tau_eval_m < 0] <- 0

  df
}

# ============================================================
# HELPERS PARA MODELO LINEAL + CÚBICO + RESTRICCIONES
# ============================================================

null_space <- function(A, tol = 1e-10) {
  s <- svd(A, nu = 0, nv = ncol(A))

  if (length(s$d) == 0) {
    return(diag(ncol(A)))
  }

  dmax <- max(s$d)
  if (!is.finite(dmax) || dmax == 0) {
    return(diag(ncol(A)))
  }

  r <- sum(s$d > tol * dmax)

  if (r >= ncol(A)) {
    stop("Las restricciones eliminan todo el espacio paramétrico.")
  }

  s$v[, seq.int(r + 1, ncol(A)), drop = FALSE]
}

make_bs_basis <- function(x, spline_df, tau_star_m, knots = NULL) {
  if (is.null(knots)) {
    B <- bs(
      x,
      df = spline_df,
      degree = 3,
      Boundary.knots = c(0, tau_star_m),
      intercept = FALSE
    )
    list(B = B, knots = attr(B, "knots"))
  } else {
    B <- bs(
      x,
      knots = knots,
      degree = 3,
      Boundary.knots = c(0, tau_star_m),
      intercept = FALSE
    )
    list(B = B, knots = knots)
  }
}

# residualiza base cúbica respecto de 1 y tau para separar parte lineal explícita
residualize_basis_against_linear <- function(B_support, tau_support, B_all, tau_all) {
  M_support <- cbind(1, tau_support)
  M_all     <- cbind(1, tau_all)

  B_nl_all <- B_all
  proj_coef <- matrix(NA_real_, nrow = ncol(M_support), ncol = ncol(B_support))

  XtX <- crossprod(M_support)

  for (j in seq_len(ncol(B_support))) {
    coef_j <- solve(XtX, crossprod(M_support, B_support[, j]))
    proj_coef[, j] <- coef_j
    B_nl_all[, j] <- B_all[, j] - M_all %*% coef_j
  }

  list(B_nl_all = B_nl_all, proj_coef = proj_coef)
}

# ── FIX 1: eval_nonlinear_basis — forzar ncol consistente ──
eval_nonlinear_basis <- function(x, knots, spline_df, tau_star_m, proj_coef) {
  
  # Clipear x para no extrapolar más allá del boundary
  x_clip <- pmin(pmax(x, 0), tau_star_m)
  
  B_raw <- make_bs_basis(
    x          = x_clip, 
    spline_df  = spline_df, 
    tau_star_m = tau_star_m, 
    knots      = knots
  )$B
  
  # Forzar dimensión correcta por si bs() colapsa alguna columna
  expected_cols <- ncol(proj_coef)
  if (ncol(B_raw) != expected_cols) {
    stop(paste0(
      "bs() devolvió ", ncol(B_raw), " columnas pero proj_coef tiene ", 
      expected_cols, ". Revisar knots o boundary."
    ))
  }
  
  M   <- cbind(1, x_clip)
  B_nl <- B_raw
  for (j in seq_len(ncol(B_raw))) {
    B_nl[, j] <- B_raw[, j] - M %*% proj_coef[, j]
  }
  B_nl
}

fit_waning_model <- function(df_prepared,
                             split_width_days,
                             spline_df,
                             tau_grid_by = 0.05,
                             tau_star_m = 15) {

  df <- df_prepared
  
  tau_eval_clip <- pmin(df$tau_eval_m, tau_star_m)

  # ── índice de filas "activas" (inmunizado y dentro del rango) ──
  active_idx <- df$inmunizado == 1 & df$tau_eval_m < tau_star_m

  tau_imm <- tau_eval_clip[active_idx]
  tau_imm <- tau_imm[is.finite(tau_imm)]

  if (length(unique(tau_imm)) <= spline_df) {
    stop(paste0(
      "No hay suficiente soporte único en tau para spline_df=", spline_df,
      " con malla ", split_width_days, " días."
    ))
  }

  # ----------------------------------------------------------
  # 1) Base cúbica cruda SOLO en soporte observado (activas)
  # ----------------------------------------------------------
  bs_support_obj <- make_bs_basis(
    x          = tau_imm,        # <-- solo filas activas
    spline_df  = spline_df,
    tau_star_m = tau_star_m,
    knots      = NULL
  )

  knots     <- bs_support_obj$knots
  B_support <- bs_support_obj$B   # nrow = length(tau_imm)

  # base cúbica cruda para TODAS las filas del modelo
  B_all <- make_bs_basis(
    x          = tau_eval_clip,
    spline_df  = spline_df,
    tau_star_m = tau_star_m,
    knots      = knots
  )$B                             # nrow = nrow(df)

  # ----------------------------------------------------------
  # 2) Residualizar: pasar soporte y tau consistentes
  # ----------------------------------------------------------
  rez <- residualize_basis_against_linear(
    B_support   = B_support,      # nrow = n_activas
    tau_support = tau_imm,        # length = n_activas  ✓
    B_all       = B_all,          # nrow = nrow(df)
    tau_all     = tau_eval_clip   # length = nrow(df)   ✓
  )

  B_nl_all  <- rez$B_nl_all
  proj_coef <- rez$proj_coef

  # ----------------------------------------------------------
  # 3) Construir diseño ORIGINAL:
  #    eta(tau)=beta0 + beta1*tau + sum gamma_j * Btilde_j(tau)
  #
  #    y para tau >= tau_star dejamos efecto = 0 (cola plana)
  # ----------------------------------------------------------
  active <- as.numeric(df$inmunizado == 1 & df$tau_eval_m < tau_star_m)

  X_orig <- cbind(
    beta0 = active,
    beta1 = active * df$tau_eval_m,
    active * B_nl_all
  )

  colnames(X_orig) <- c(
    "beta0",
    "beta1",
    paste0("gamma", seq_len(ncol(B_nl_all)))
  )

  # ----------------------------------------------------------
  # 4) Restricciones lineales:
  #    eta(tau_star)=0
  #    eta'(tau_star)=0
  # ----------------------------------------------------------

  tau_at_star <- tau_star_m - 1e-6   # justo antes del boundary, evita extrapolación

  B_at_star   <- eval_nonlinear_basis(
    x          = tau_at_star,
    knots      = knots,
    spline_df  = spline_df,
    tau_star_m = tau_star_m,
    proj_coef  = proj_coef
  )

  eps <- 1e-6
  B_at_star_p <- eval_nonlinear_basis(
    x          = tau_at_star + eps,
    knots      = knots,
    spline_df  = spline_df,
    tau_star_m = tau_star_m,
    proj_coef  = proj_coef
  )

  B_at_star_m <- eval_nonlinear_basis(
    x          = tau_at_star - eps,
    knots      = knots,
    spline_df  = spline_df,
    tau_star_m = tau_star_m,
    proj_coef  = proj_coef
  )

  dB_at_star <- (B_at_star_p - B_at_star_m) / (2 * eps)

  A <- rbind(
    c(1, tau_at_star, as.numeric(B_at_star)),
    c(0, 1,           as.numeric(dB_at_star))
  )

  colnames(A) <- colnames(X_orig)
  rownames(A) <- c("value_at_tau0", "slope_at_tau0")

  # ----------------------------------------------------------
  # 5) Null-space reparametrization: theta = N alpha
  # ----------------------------------------------------------
  Nmat <- null_space(A)

  Z <- X_orig %*% Nmat
  colnames(Z) <- paste0("z", seq_len(ncol(Z)))

  df_fit <- cbind(df, as.data.frame(Z))

  rhs_terms <- c(
    colnames(Z),
    "strata(Group)",
    "cluster(RUN)"
  )

  rhs_ns <- paste(rhs_terms, collapse = " + ")
  f_ns <- as.formula(paste0("Surv(start, stop, event_vrs) ~ ", rhs_ns))

  fit_ns <- coxph(
    formula = f_ns,
    data    = df_fit,
    ties    = "efron",
    robust  = TRUE,
    x       = TRUE,
    model   = TRUE
  )

  fit_const <- coxph(
    Surv(start, stop, event_vrs) ~ inmunizado + strata(Group) + cluster(RUN),
    data   = df,
    ties   = "efron",
    robust = TRUE,
    x      = TRUE,
    model  = TRUE
  )

  # ----------------------------------------------------------
  # 6) Recuperar coeficientes originales theta
  # ----------------------------------------------------------
  alpha_full <- coef(fit_ns)
  keep_alpha <- !is.na(alpha_full)

  alpha <- alpha_full[keep_alpha]
  V_alpha <- vcov(fit_ns)[keep_alpha, keep_alpha, drop = FALSE]

  N_keep <- Nmat[, keep_alpha, drop = FALSE]

  theta_hat <- as.vector(N_keep %*% alpha)
  names(theta_hat) <- colnames(X_orig)

  V_theta <- N_keep %*% V_alpha %*% t(N_keep)
  rownames(V_theta) <- colnames(V_theta) <- colnames(X_orig)

  # ----------------------------------------------------------
  # 7) Predicción en grilla
  # ----------------------------------------------------------
  tau_grid_m <- seq(
    0,
    ceiling(max(df$tau_eval_m[df$inmunizado == 1], na.rm = TRUE)),
    by = tau_grid_by
  )

  B_grid_raw <- make_bs_basis(
    x = tau_grid_m,
    spline_df = spline_df,
    tau_star_m = tau_star_m,
    knots = knots
  )$B

  B_grid_nl <- eval_nonlinear_basis(
    x = tau_grid_m,
    knots = knots,
    spline_df = spline_df,
    tau_star_m = tau_star_m,
    proj_coef = proj_coef
  )

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
    HR               = HR,
    HR_lo            = HR_lo,
    HR_hi            = HR_hi,
    VE               = VE,
    VE_lo            = VE_lo,
    VE_hi            = VE_hi
  )

  # ----------------------------------------------------------
  # 8) Soporte
  # ----------------------------------------------------------
  support <- subset(df, inmunizado == 1)
  support$month_bin <- floor(support$tau_eval_m)

  support_summary <- aggregate(
    cbind(person_time_days = stop - start, events = event_vrs) ~ month_bin,
    data = support,
    FUN = sum
  )

  n_rows_tab <- table(support$month_bin)
  support_summary$n_rows <- as.numeric(n_rows_tab[as.character(support_summary$month_bin)])
  support_summary$split_width_days <- split_width_days
  support_summary$spline_df <- spline_df
  support_summary$tau_0 <- tau_star_m

  # ----------------------------------------------------------
  # 9) Summary
  # ----------------------------------------------------------
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
      df = sum(!is.na(alpha_full)) - sum(!is.na(coef(fit_const))),
      lower.tail = FALSE
    )
  )

  list(
    df_model         = df_fit,
    fit_ns           = fit_ns,
    fit_const        = fit_const,
    ve_curve         = ve_curve,
    support_summary  = support_summary,
    model_summary    = model_summary,
    theta_hat        = theta_hat,
    V_theta          = V_theta,
    A_constraints    = A,
    N_null           = Nmat,
    knots            = knots,
    proj_coef        = proj_coef
  )
}
# ============================================================
# 3) PREPARAR LAS 3 MALLAS UNA SOLA VEZ
# ============================================================

split_grid  <- c(14) #c(7, 14, 30)
spline_grid <- c(3,4,5,6) #c(3, 4, 5)

tau_star_gird <- c(13,14,15)

prepared_by_width <- list()

for (w in split_grid) {
  message("Preparando malla de ", w, " días...")
  tmp <- rebuild_intervals(df_raw, split_width_days = w)
  tmp <- add_tau_variables(tmp)
  prepared_by_width[[paste0("w", w)]] <- tmp
}

# ============================================================
# 4) FIT DE LOS 9 MODELOS
# ============================================================
results <- list()

for (w in split_grid) {
  for (sdf in spline_grid) {
    for (tau_0 in tau_star_gird) {
    key <- paste0("w", w, "_df", sdf, "_tau_",tau_0)
    message("Ajustando modelo: malla=", w, " días ; spline_df=", sdf, ' tau=', tau_0)

    results[[key]] <- fit_waning_model(
      df_prepared      = prepared_by_width[[paste0("w", w)]],
      split_width_days = w,
      spline_df        = sdf,
      tau_grid_by      = 0.05,
      tau_star_m       = tau_0
      # use_tail_constraint = TRUE
    )
  }
  }
}

# ============================================================
# 5) TABLAS CONSOLIDADAS
# ============================================================
# all_curves <- bind_rows(lapply(results, function(x) x$ve_curve))


all_curves <- imap_dfr(results, function(res, key) {
  
  m <- str_match(key, "^w(\\d+)_df(\\d+)_tau_(\\d+)$")
  
  curve_df <- res$ve_curve %>%   # cambia esto si tu objeto usa otro nombre
    mutate(
      split_width_days = as.integer(m[2]),
      spline_df        = as.integer(m[3]),
      tau_0            = as.integer(m[4]),
      model_key        = key
    )
  
  curve_df
})

all_support <- bind_rows(lapply(results, function(x) x$support_summary))
all_model_summary <- bind_rows(lapply(results, function(x) x$model_summary))

print(all_model_summary)

breaks_01 <- function(x) {
  lo <- floor((min(x, na.rm = TRUE) + 1e-8) * 10) / 10
  hi <- ceiling((max(x, na.rm = TRUE) - 1e-8) * 10) / 10
  seq(lo, hi, by = 0.1)
}

labels_01 <- function(x) {
  x[abs(x) < 1e-8] <- 0
  sprintf("%.1f", x)
}


p_tau_grid <- all_curves %>%
  mutate(
    split_width_days = factor(split_width_days, levels = c(14)),
    spline_df = factor(spline_df, levels = c(3, 4, 5,6),
                       labels = c("ns(df = 3)", "ns(df = 4)", "ns(df = 5)", "ns(df = 6)")),
    tau_0 = factor(tau_0, levels = c(13, 14, 15),
                   labels = c(expression(tau[0] == 13),
                              expression(tau[0] == 14),
                              expression(tau[0] == 15)))
  ) %>%
  ggplot(aes(x = tau_m, y = VE, color = split_width_days, fill = split_width_days)) +
  geom_ribbon(aes(ymin = VE_lo, ymax = VE_hi), alpha = 0.08, linewidth = 0) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_grid(rows = vars(spline_df), cols = vars(tau_0), scales = "free_y", labeller = label_parsed) +
  scale_y_continuous(
    breaks = c(0, 0.5, 1),
    labels = labels_01
  ) +
  coord_cartesian(xlim = c(0, 15)) +
  labs(
    x = "Months since immunization",
    y = "VE(t) = 1 - HR(t)",
    color = "Split width (days)",
    fill  = "Split width (days)",
    title = "Sensitivity to spline flexibility and tail constraint"
  ) +
  theme_bw(base_size = 11) +
  theme(
    aspect.ratio = 0.55,
    legend.position = "bottom"
  )

print(p_tau_grid)




# ============================================================
# FIGURA PRINCIPAL: VE(t) — paper-ready
# ============================================================
library(ggplot2)
library(dplyr)

# Extraer resultado del modelo seleccionado
key_sel  <- "w14_df3_tau_15"
res_sel  <- results[[key_sel]]
curve    <- res_sel$ve_curve
support  <- res_sel$support_summary

# Histograma de soporte (persona-tiempo por mes) para fondo
support_hist <- support %>%
  mutate(rel_pt = person_time_days / max(person_time_days, na.rm = TRUE) * 0.18)

# Paleta
col_line  <- "#1B4F72"
col_band  <- "#2E86C1"
col_zero  <- "grey40"
col_hist  <- "grey82"

p_main <- ggplot() +
  # Fondo: barras de persona-tiempo (escaladas al eje VE)
  geom_rect(
    data = support_hist,
    aes(xmin = month_bin, xmax = month_bin + 0.9,
        ymin = -0.04, ymax = -0.04 + rel_pt),
    fill  = col_hist,
    color = NA,
    alpha = 0.9
  ) +
  # Línea de referencia VE = 0
  geom_hline(yintercept = 0, linetype = "dashed",
             color = col_zero, linewidth = 0.45) +
  # Banda de confianza
  geom_ribbon(
    data = curve,
    aes(x = tau_m, ymin = VE_lo, ymax = VE_hi),
    fill  = col_band,
    alpha = 0.15
  ) +
  # Línea central VE
  geom_line(
    data = curve,
    aes(x = tau_m, y = VE),
    color     = col_line,
    linewidth = 1.1
  ) +
  # Marca en tau_star
  geom_vline(xintercept = 15, linetype = "dotted",
             color = "grey50", linewidth = 0.5) +
  annotate("text", x = 15.15, y = 0.92, label = expression(tau[0] == 15),
           hjust = 0, size = 3.2, color = "grey40", family = "serif") +
  # Anotación de soporte
  annotate("text", x = 0.1, y = -0.01,
           label = "Person-time\n(relative)",
           hjust = 0, vjust = 1, size = 2.6,
           color = "grey50", family = "serif") +
  scale_x_continuous(
    breaks = 0:15,
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(
    breaks = c(-0.25, 0, 0.25, 0.5, 0.75, 1.0),
    labels = scales::percent_format(accuracy = 1),
    limits = c(-0.27, 1.02)
  ) +
  coord_cartesian(xlim = c(0, 15.6)) +
  labs(
    x     = "Months since immunization",
    y     = "Vaccine effectiveness  VE(t) = 1 − HR(t)",
    title = "Waning of vaccine effectiveness over time",
    subtitle = "Cox model with constrained cubic spline  ·  split = 14 d, df = 3, τ₀ = 15 mo",
    caption  = "Shaded area: pointwise 95% CI.  Bars: relative person-time at risk (immunized arm)."
  ) +
  theme_classic(base_size = 12, base_family = "serif") +
  theme(
    plot.title    = element_text(face = "bold", size = 20),
    plot.subtitle = element_text(color = "grey35", size = 9.5),
    plot.caption  = element_text(color = "grey45", size = 8, hjust = 0),
    axis.title    = element_text(size = 10.5),
    axis.text     = element_text(size = 9.5),
    panel.grid.major.y = element_line(color = "grey93", linewidth = 0.35),
    panel.grid.major.x = element_blank(),
    plot.margin   = margin(12, 16, 8, 8)
  )

print(p_main)

# Guardar
ggsave("ve_curve_paper.pdf", p_main, width = 6.5, height = 4.2, device = cairo_pdf)
ggsave("ve_curve_paper.png", p_main, width = 6.5, height = 4.2, dpi = 320)


# ============================================================
# FUNCIÓN GENÉRICA — solo cambia make_f_basis()
# ============================================================

# Cada "familia" es una función que recibe tau y devuelve
# una matriz de columnas [f1(tau), f2(tau), ...]
# más su derivada para la restricción R2.

# ── Familia spline (lo que tienes ahora) ──
make_f_spline <- function(tau, tau_star_m, spline_df, knots = NULL) {
  obj <- make_bs_basis(tau, spline_df = spline_df,
                       tau_star_m = tau_star_m, knots = knots)
  list(B = obj$B, knots = obj$knots)
}

# ── Familia logarítmica: f(tau) = log(1 + tau) ──
make_f_log <- function(tau, tau_star_m, ...) {
  B <- matrix(log(1 + tau), ncol = 1)
  list(B = B, knots = NULL)
}

# ── Familia Weibull: f(tau) = tau^gamma  (gamma fijo) ──
make_f_weibull <- function(tau, tau_star_m, gamma = 2, ...) {
  B <- matrix(tau^gamma, ncol = 1)
  list(B = B, knots = NULL)
}

# ── Familia exponencial: f(tau) = exp(-lambda * tau) ──
make_f_exponential <- function(tau, tau_star_m, lambda = 0.3, ...) {
  B <- matrix(exp(-lambda * tau), ncol = 1)
  list(B = B, knots = NULL)
}

# ── Familia potencia: f(tau) = tau^p  (p fijo) ──
make_f_power <- function(tau, tau_star_m, p = 0.5, ...) {
  B <- matrix(tau^p, ncol = 1)
  list(B = B, knots = NULL)
}

# ============================================================
# FIT GENÉRICO — recibe cualquier make_f_*
# ============================================================

fit_waning_generic <- function(df_prepared,
                                split_width_days,
                                tau_star_m    = 15,
                                tau_grid_by   = 0.05,
                                f_basis_fn    = make_f_spline,
                                f_basis_args  = list(spline_df = 3),
                                family_label  = "spline") {
  df <- df_prepared
  tau_eval_clip <- pmin(df$tau_eval_m, tau_star_m)
  active <- as.numeric(df$inmunizado == 1 & df$tau_eval_m < tau_star_m)

  # ── 1) Construir f(tau) en soporte activo para estimar knots ──
  tau_imm <- tau_eval_clip[df$inmunizado == 1 & df$tau_eval_m < tau_star_m]
  tau_imm <- tau_imm[is.finite(tau_imm)]

  f_support <- do.call(f_basis_fn, c(list(tau = tau_imm,
                                           tau_star_m = tau_star_m),
                                      f_basis_args))
  knots_out <- f_support$knots   # NULL para familias paramétricas

  # ── 2) f(tau) para todas las filas ──
  f_all <- do.call(f_basis_fn, c(list(tau = tau_eval_clip,
                                       tau_star_m = tau_star_m),
                                  c(f_basis_args, list(knots = knots_out))))$B

  n_f <- ncol(f_all)
  f_names <- paste0("f", seq_len(n_f))

  # ── 3) Matriz de diseño completa ──
  # log-HR(tau) = beta0*active + beta1*active*tau + beta2*active*f1(tau) + ...
  X_orig <- cbind(
    beta0 = active,
    beta1 = active * tau_eval_clip,
    sweep(f_all, 1, active, `*`)   # active * f_j(tau) para cada columna j
  )
  colnames(X_orig) <- c("beta0", "beta1", f_names)

  # ── 4) Restricciones A*theta = 0 en tau_star ──
  eps         <- 1e-6
  tau_at_star <- tau_star_m - eps

  eval_f_at <- function(t) {
    do.call(f_basis_fn, c(list(tau = t, tau_star_m = tau_star_m),
                           c(f_basis_args, list(knots = knots_out))))$B
  }

  f_star   <- eval_f_at(tau_at_star)           # valor en tau*
  f_star_p <- eval_f_at(tau_at_star + eps)     # para derivada
  f_star_m <- eval_f_at(tau_at_star - eps)
  df_star  <- (f_star_p - f_star_m) / (2 * eps)  # derivada numérica

  # R1: beta0 + beta1*tau* + sum(beta_j * f_j(tau*)) = 0
  # R2: beta1 + sum(beta_j * f_j'(tau*)) = 0
  A <- rbind(
    c(1, tau_at_star, as.numeric(f_star)),
    c(0, 1,           as.numeric(df_star))
  )
  colnames(A) <- colnames(X_orig)
  rownames(A) <- c("value_at_tau0", "slope_at_tau0")

  # ── 5) Null-space + ajuste ──
  Nmat <- null_space(A)
  Z    <- X_orig %*% Nmat
  colnames(Z) <- paste0("z", seq_len(ncol(Z)))

  df_fit <- cbind(df, as.data.frame(Z))

  f_ns <- as.formula(paste0(
    "Surv(start, stop, event_vrs) ~ ",
    paste(c(colnames(Z), "strata(Group)", "cluster(RUN)"), collapse = " + ")
  ))

  fit_ns <- coxph(f_ns, data = df_fit, ties = "efron",
                  robust = TRUE, x = TRUE, model = TRUE)

  fit_const <- coxph(
    Surv(start, stop, event_vrs) ~ inmunizado + strata(Group) + cluster(RUN),
    data = df, ties = "efron", robust = TRUE
  )

  # ── 6) Recuperar theta ──
  alpha_full <- coef(fit_ns)
  keep_alpha <- !is.na(alpha_full)
  alpha      <- alpha_full[keep_alpha]
  V_alpha    <- vcov(fit_ns)[keep_alpha, keep_alpha, drop = FALSE]
  N_keep     <- Nmat[, keep_alpha, drop = FALSE]

  theta_hat <- as.vector(N_keep %*% alpha)
  names(theta_hat) <- colnames(X_orig)
  V_theta <- N_keep %*% V_alpha %*% t(N_keep)
  rownames(V_theta) <- colnames(V_theta) <- colnames(X_orig)

  # ── 7) Predicción en grilla ──
  tau_grid_m <- seq(0,
    ceiling(max(df$tau_eval_m[df$inmunizado == 1], na.rm = TRUE)),
    by = tau_grid_by)

  f_grid <- do.call(f_basis_fn, c(list(tau = tau_grid_m,
                                        tau_star_m = tau_star_m),
                                   c(f_basis_args, list(knots = knots_out))))$B

  active_g <- as.numeric(tau_grid_m < tau_star_m)

  Xg <- cbind(
    beta0 = active_g,
    beta1 = active_g * tau_grid_m,
    sweep(f_grid, 1, active_g, `*`)
  )
  colnames(Xg) <- colnames(X_orig)

  eta    <- as.vector(Xg %*% theta_hat)
  se_eta <- sqrt(rowSums((Xg %*% V_theta) * Xg))

  ve_curve <- data.frame(
    family   = family_label,
    tau_m    = tau_grid_m,
    logHR    = eta,
    se       = se_eta,
    HR       = exp(eta),
    HR_lo    = exp(eta - 1.96 * se_eta),
    HR_hi    = exp(eta + 1.96 * se_eta),
    VE       = 1 - exp(eta),
    VE_lo    = 1 - exp(eta + 1.96 * se_eta),
    VE_hi    = 1 - exp(eta - 1.96 * se_eta)
  )

  # ── 8) Model summary ──
  model_summary <- data.frame(
    family           = family_label,
    tau_0            = tau_star_m,
    n_params_free    = sum(keep_alpha),
    loglik_ns        = fit_ns$loglik[2],
    loglik_const     = fit_const$loglik[2],
    LRT_stat         = 2 * (fit_ns$loglik[2] - fit_const$loglik[2]),
    AIC              = -2 * fit_ns$loglik[2] + 2 * sum(keep_alpha),
    p_value          = pchisq(
      2 * (fit_ns$loglik[2] - fit_const$loglik[2]),
      df = sum(keep_alpha) - sum(!is.na(coef(fit_const))),
      lower.tail = FALSE
    )
  )

  list(
    fit_ns        = fit_ns,
    fit_const     = fit_const,
    ve_curve      = ve_curve,
    model_summary = model_summary,
    theta_hat     = theta_hat,
    V_theta       = V_theta,
    A_constraints = A,
    N_null        = Nmat,
    knots         = knots_out
  )
}

# ============================================================
# CORRER Y COMPARAR LAS 4 FAMILIAS
# ============================================================
data_w14 <- prepared_by_width[["w14"]]

res_families <- list(
  spline = fit_waning_generic(
    data_w14, split_width_days = 14, tau_star_m = 15,
    f_basis_fn   = make_f_spline,
    f_basis_args = list(spline_df = 3),
    family_label = "Spline (df=3)"
  ),
  log = fit_waning_generic(
    data_w14, split_width_days = 14, tau_star_m = 15,
    f_basis_fn   = make_f_log,
    f_basis_args = list(),
    family_label = "Logarítmica"
  ),
  weibull = fit_waning_generic(
    data_w14, split_width_days = 14, tau_star_m = 15,
    f_basis_fn   = make_f_weibull,
    f_basis_args = list(gamma = 2),
    family_label = "Weibull (γ=2)"
  ),
  exponential = fit_waning_generic(
    data_w14, split_width_days = 14, tau_star_m = 15,
    f_basis_fn   = make_f_exponential,
    f_basis_args = list(lambda = 0.3),
    family_label = "Exponencial (λ=0.3)"
  )
)

# Tabla comparativa
comp_table <- bind_rows(lapply(res_families, `[[`, "model_summary")) %>%
  mutate(delta_AIC = AIC - min(AIC)) %>%
  arrange(AIC)

print(comp_table)

df_curves <- bind_rows(lapply(res_families, `[[`, "ve_curve"))

# ── Paleta y orden ──
familia_colors <- c(
  "Spline (df=3)"      = "#1B4F72",
  "Logarítmica"        = "#C0392B",
  "Weibull (γ=2)"      = "#1E8449",
  "Exponencial (λ=0.3)"= "#7D3C98"
)

# ============================================================
# FIGURA 1: Un solo panel, sin IC, todas las curvas juntas
# ============================================================
p_overlay <- df_curves %>%
  ggplot(aes(tau_m, VE, color = family)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = familia_colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(breaks = seq(0, 15, by = 3)) +
  coord_cartesian(xlim = c(0, 15), ylim = c(-0.05, 1.05)) +
  labs(
    x        = "Months since immunization",
    y        = "VE(τ) = 1 − HR(τ)",
    color    = NULL,
    title    = "Waning of vaccine effectiveness — functional form comparison",
    subtitle = "Point estimates only  ·  τ₀ = 15 mo  ·  split = 14 d"
  ) +
  theme_classic(base_size = 15, base_family = "serif") +
  theme(
    plot.title       = element_text(face = "bold", size = 20),
    plot.subtitle    = element_text(color = "grey35", size = 15),  # más grande
    legend.position  = "bottom",
    legend.text      = element_text(size = 12),
    legend.key.width = unit(1.1, "cm"),                              # líneas más largas en leyenda
    panel.grid.major.y = element_line(color = "grey93", linewidth = 0.3),
    plot.margin      = margin(12, 16, 8, 10),
    aspect.ratio     = 0.55                                          # más ancho que alto
  )

print(p_overlay)


library(dplyr)
library(knitr)
library(kableExtra)

# ============================================================
# EXTRAER COEFICIENTES DE CADA MODELO
# ============================================================

extract_coef_block <- function(res, family_label) {
  th   <- res$theta_hat
  vth  <- res$V_theta
  se_t <- sqrt(diag(vth))
  z_t  <- th / se_t
  p_t  <- 2 * pnorm(abs(z_t), lower.tail = FALSE)
  
  # Estrellas
  stars <- case_when(
    p_t < 0.001 ~ "***",
    p_t < 0.01  ~ "**",
    p_t < 0.05  ~ "*",
    TRUE        ~ ""
  )
  
  # Formato: "estimate (SE)***"
  cell <- paste0(
    sprintf("%.3f", th), stars, "\n",
    "(", sprintf("%.3f", se_t), ")"
  )
  
  data.frame(
    param  = names(th),
    value  = cell,
    stringsAsFactors = FALSE
  ) %>% setNames(c("param", family_label))
}

# Extraer de cada familia
tabs <- lapply(names(res_families), function(nm) {
  extract_coef_block(res_families[[nm]],
                     res_families[[nm]]$model_summary$family)
})

# Unir por param (full join para manejar distintos números de params)
tab_coef <- Reduce(function(a, b) full_join(a, b, by = "param"), tabs)

# Ordenar filas: beta0, beta1, f1, f2, ...
param_order <- c("beta0", "beta1",
                 paste0("f", 1:10))
tab_coef <- tab_coef %>%
  mutate(param = factor(param, levels = param_order)) %>%
  arrange(param) %>%
  mutate(param = as.character(param))

# Etiquetas más bonitas para las filas
param_labels <- c(
  beta0 = "β₀  (intercept)",
  beta1 = "β₁  (linear slope)",
  f1    = "f₁(τ)",
  f2    = "f₂(τ)",
  f3    = "f₃(τ)"
)
tab_coef$param <- dplyr::recode(tab_coef$param, !!!param_labels)

# ============================================================
# FILA DE MÉTRICAS (pie de tabla)
# ============================================================

metrics <- bind_rows(lapply(res_families, `[[`, "model_summary")) %>%
  mutate(
    AIC     = round(-2 * loglik_ns + 2 * n_params_free, 2),
    loglik  = round(loglik_ns, 2),
    LRT_p   = ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))
  )

# Una fila por métrica, una columna por modelo
make_metric_row <- function(label, values) {
  row <- c(label, as.character(values))
  setNames(as.data.frame(t(row), stringsAsFactors = FALSE),
           c("param", colnames(tab_coef)[-1]))
}

rows_metrics <- bind_rows(
  make_metric_row("N (free params)",  metrics$n_params_free),
  make_metric_row("Log-likelihood",   metrics$loglik),
  make_metric_row("AIC",              metrics$AIC),
  make_metric_row("LRT p-value",      metrics$LRT_p)
)

# NA → "—"
tab_coef[is.na(tab_coef)] <- "—"

# Tabla final
tab_final <- bind_rows(tab_coef, rows_metrics)

# ============================================================
# RENDER
# ============================================================
n_coef <- nrow(tab_coef)
n_met  <- nrow(rows_metrics)
col_names <- c("", colnames(tab_final)[-1])

tab_final %>%
  kable(
    format    = "html",
    col.names = col_names,
    align     = c("l", rep("c", ncol(tab_final) - 1)),
    caption   = "Model coefficients by functional form  ·  β₀ + β₁τ + β₂f(τ)  ·  τ₀ = 15 mo",
    escape    = FALSE
  ) %>%
  kable_styling(
    full_width        = FALSE,
    bootstrap_options = c("striped", "hover", "condensed"),
    font_size         = 13
  ) %>%
  # Separar coeficientes de métricas
  row_spec(n_coef, extra_css = "border-bottom: 2px solid #555;") %>%
  # Sombrear filas de métricas
  row_spec(
    seq(n_coef + 1, n_coef + n_met),
    background = "#f5f5f5",
    italic     = TRUE
  ) %>%
  # Cabecera en negrita
  row_spec(0, bold = TRUE, background = "#2C3E50", color = "white") %>%
  # Columna de parámetros en negrita
  column_spec(1, bold = TRUE, width = "12em") %>%
  # Guardar
  {
    html_out <- as.character(.)
    writeLines(paste0(
      '<!DOCTYPE html><html><head><meta charset="utf-8">',
      '<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css">',
      '<style>
         body { padding: 30px; font-family: Georgia, serif; }
         td { white-space: pre-line; vertical-align: middle; }
         caption { font-weight: bold; font-size: 14px; margin-bottom: 8px; }
       </style>',
      '</head><body>', html_out, '</body></html>'
    ), "coef_comparison_table.html")
    cat("✓ Tabla guardada en coef_comparison_table.html\n")
    invisible(.)
  }

# También imprimir en consola (versión simple)
cat("\n")
tab_final %>%
  kable("simple", col.names = col_names) %>%
  print()

# Nota al pie para la consola
cat("\nNota: estimate (SE). Significancia: * p<0.05  ** p<0.01  *** p<0.001\n")
cat("'—' indica que el parámetro no existe en ese modelo.\n")

# Guardar con dimensiones horizontales
ggsave("ve_overlay.pdf", p_overlay, width = 7.5, height = 5.0, device = cairo_pdf)
ggsave("ve_overlay.png", p_overlay, width = 7.5, height = 5.0, dpi = 320)

# ============================================================
# FIGURA 2: Facetas con IC
# ============================================================
p_facet <- df_curves %>%
  ggplot(aes(tau_m, VE, color = family, fill = family)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey55", linewidth = 0.4) +
  geom_ribbon(aes(ymin = VE_lo, ymax = VE_hi),
              alpha = 0.15, linewidth = 0) +
  geom_line(linewidth = 0.95) +
  facet_wrap(~ family, ncol = 2) +
  scale_color_manual(values = familia_colors) +
  scale_fill_manual(values  = familia_colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(breaks = seq(0, 15, by = 5)) +
  coord_cartesian(xlim = c(0, 15), ylim = c(-0.35, 1.1)) +
  labs(
    x        = "Months since immunization",
    y        = "VE(τ) = 1 − HR(τ)",
    title    = "Waning of vaccine effectiveness — functional form comparison",
    subtitle = "Shaded area: pointwise 95% CI  ·  τ₀ = 15 mo  ·  split = 14 d"
  ) +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(color = "grey40", size = 9),
    legend.position = "none",          # color ya está en el facet label
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    strip.text    = element_text(face = "bold", size = 9.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey93", linewidth = 0.3),
    plot.margin   = margin(10, 14, 6, 8)
  )

print(p_facet)

# ── Guardar ambas ──
ggsave("ve_overlay.pdf",  p_overlay, width = 6.5, height = 4.0, device = cairo_pdf)
ggsave("ve_overlay.png",  p_overlay, width = 6.5, height = 4.0, dpi = 320)
ggsave("ve_facet.pdf",    p_facet,   width = 7.5, height = 6.0, device = cairo_pdf)
ggsave("ve_facet.png",    p_facet,   width = 7.5, height = 6.0, dpi = 320)

# ============================================================
# PROFILE LIKELIHOOD GENÉRICO
# ============================================================

profile_shape <- function(df_prepared,
                          tau_star_m   = 15,
                          f_basis_fn   = make_f_weibull,
                          shape_name   = "gamma",
                          shape_grid   = seq(0.3, 4, by = 0.1),
                          family_label = "Weibull") {

  logliks <- numeric(length(shape_grid))

  for (i in seq_along(shape_grid)) {
    val <- shape_grid[i]
    args <- setNames(list(val), shape_name)

    res_i <- tryCatch(
      fit_waning_generic(
        df_prepared  = df_prepared,
        split_width_days = 14,
        tau_star_m   = tau_star_m,
        f_basis_fn   = f_basis_fn,
        f_basis_args = args,
        family_label = family_label
      ),
      error = function(e) NULL
    )

    logliks[i] <- if (is.null(res_i)) NA_real_ else res_i$fit_ns$loglik[2]
  }

  df_profile <- data.frame(
    shape  = shape_grid,
    loglik = logliks,
    delta  = logliks - max(logliks, na.rm = TRUE)  # relativo al máximo
  )

  # Estimador puntual
  best_idx <- which.max(logliks)
  best_val <- shape_grid[best_idx]

  # IC aproximado al 95%: delta > -1.92  (chi²₁/2)
  ci <- df_profile %>%
    filter(delta > -1.92) %>%
    summarise(lo = min(shape), hi = max(shape))

  list(
    profile   = df_profile,
    best_val  = best_val,
    best_ll   = logliks[best_idx],
    ci_lo     = ci$lo,
    ci_hi     = ci$hi,
    shape_name = shape_name,
    family    = family_label
  )
}

# ── Correr para Weibull (γ) y Exponencial (λ) ──
prof_wb  <- profile_shape(
  prepared_by_width[["w14"]],
  tau_star_m   = 15,
  f_basis_fn   = make_f_weibull,
  shape_name   = "gamma",
  shape_grid   = seq(0.3, 5, by = 0.1),
  family_label = "Weibull"
)

prof_exp <- profile_shape(
  prepared_by_width[["w14"]],
  tau_star_m   = 15,
  f_basis_fn   = make_f_exponential,
  shape_name   = "lambda",
  shape_grid   = seq(0.05, 1.5, by = 0.05),
  family_label = "Exponencial"
)

# Resultados
cat("── Weibull ──\n")
cat(sprintf("  γ̂ = %.2f   95%% CI: [%.2f, %.2f]\n",
            prof_wb$best_val, prof_wb$ci_lo, prof_wb$ci_hi))

cat("── Exponencial ──\n")
cat(sprintf("  λ̂ = %.2f   95%% CI: [%.2f, %.2f]\n",
            prof_exp$best_val, prof_exp$ci_lo, prof_exp$ci_hi))

# ============================================================
# GRÁFICO DE PROFILE LIKELIHOOD
# ============================================================
plot_profile <- function(prof) {
  ggplot(prof$profile, aes(shape, delta)) +
    geom_line(color = "#1B4F72", linewidth = 1) +
    geom_hline(yintercept = -1.92, linetype = "dashed",
               color = "#C0392B", linewidth = 0.5) +
    geom_vline(xintercept = prof$best_val, linetype = "dotted",
               color = "grey40", linewidth = 0.5) +
    annotate("text",
             x = prof$best_val, y = -0.1,
             label = sprintf("%s = %.2f", prof$shape_name, prof$best_val),
             hjust = -0.1, size = 3.5, family = "serif") +
    annotate("text",
             x = min(prof$profile$shape, na.rm = TRUE),
             y = -1.75,
             label = "95% CI threshold",
             hjust = 0, size = 3, color = "#C0392B", family = "serif") +
    labs(
      x     = prof$shape_name,
      y     = "Δ log-likelihood",
      title = sprintf("Profile likelihood — %s", prof$family)
    ) +
    theme_classic(base_size = 11, base_family = "serif") +
    theme(plot.title = element_text(face = "bold"))
}

library(patchwork)
p_profiles <- plot_profile(prof_wb) | plot_profile(prof_exp)
print(p_profiles)

# ============================================================
# REAJUSTAR CON γ̂ Y λ̂ ÓPTIMOS
# ============================================================
res_wb_opt <- fit_waning_generic(
  prepared_by_width[["w14"]],
  split_width_days = 14,
  tau_star_m   = 15,
  f_basis_fn   = make_f_weibull,
  f_basis_args = list(gamma = prof_wb$best_val),
  family_label = sprintf("Weibull (γ̂=%.2f)", prof_wb$best_val)
)

res_exp_opt <- fit_waning_generic(
  prepared_by_width[["w14"]],
  split_width_days = 14,
  tau_star_m   = 15,
  f_basis_fn   = make_f_exponential,
  f_basis_args = list(lambda = prof_exp$best_val),
  family_label = sprintf("Exponencial (λ̂=%.2f)", prof_exp$best_val)
)

# Comparación final incluyendo modelos optimizados
res_all <- c(res_families, list(wb_opt = res_wb_opt, exp_opt = res_exp_opt))

bind_rows(lapply(res_all, `[[`, "model_summary")) %>%
  mutate(AIC = round(-2 * loglik_ns + 2 * n_params_free, 2),
         delta_AIC = round(AIC - min(AIC), 2)) %>%
  arrange(AIC) %>%
  kable("simple")

# ============================================================
# TABLA DE COEFICIENTES
# ============================================================
library(knitr)
# install.packages("kableExtra")  # si no lo tienes
library(kableExtra)

res_sel <- results[["w14_df3_tau_15"]]

# ── (a) Coeficientes reparametrizados (z1, z2, ...) ──
alpha_v  <- coef(res_sel$fit_ns)
se_v     <- sqrt(diag(vcov(res_sel$fit_ns)))
z_stat   <- alpha_v / se_v
p_v      <- 2 * pnorm(abs(z_stat), lower.tail = FALSE)
ci_lo    <- alpha_v - 1.96 * se_v
ci_hi    <- alpha_v + 1.96 * se_v

tab_alpha <- data.frame(
  Parameter   = names(alpha_v),
  Estimate    = round(alpha_v, 4),
  SE          = round(se_v, 4),
  `95% CI`    = paste0("[", round(ci_lo, 3), ", ", round(ci_hi, 3), "]"),
  `z value`   = round(z_stat, 3),
  `Pr(>|z|)`  = ifelse(p_v < 0.001, "<0.001", round(p_v, 3)),
  check.names = FALSE
)

# ── (b) Coeficientes originales theta (beta0, beta1, gamma_j) ──
th   <- res_sel$theta_hat
vth  <- res_sel$V_theta
se_t <- sqrt(diag(vth))
z_t  <- th / se_t
p_t  <- 2 * pnorm(abs(z_t), lower.tail = FALSE)

tab_theta <- data.frame(
  Parameter  = names(th),
  Estimate   = round(th, 4),
  SE         = round(se_t, 4),
  `95% CI`   = paste0("[", round(th - 1.96*se_t, 3), ", ",
                           round(th + 1.96*se_t, 3), "]"),
  `z value`  = round(z_t, 3),
  `Pr(>|z|)` = ifelse(p_t < 0.001, "<0.001", round(p_t, 3)),
  check.names = FALSE
)

# ── (c) Modelo constante (referencia) ──
alpha_c  <- coef(res_sel$fit_const)
se_c     <- sqrt(diag(vcov(res_sel$fit_const)))
z_c      <- alpha_c / se_c
p_c      <- 2 * pnorm(abs(z_c), lower.tail = FALSE)

tab_const <- data.frame(
  Parameter  = names(alpha_c),
  Estimate   = round(alpha_c, 4),
  SE         = round(se_c, 4),
  `95% CI`   = paste0("[", round(alpha_c - 1.96*se_c, 3), ", ",
                           round(alpha_c + 1.96*se_c, 3), "]"),
  `z value`  = round(z_c, 3),
  `Pr(>|z|)` = ifelse(p_c < 0.001, "<0.001", round(p_c, 3)),
  check.names = FALSE
)

# ── (d) LRT ──
ms  <- res_sel$model_summary
cat("\n── Likelihood-ratio test (spline vs constante) ──\n")
cat(sprintf("  log-lik spline   : %.4f\n", ms$loglik_ns))
cat(sprintf("  log-lik constante: %.4f\n", ms$loglik_const))
cat(sprintf("  LRT stat (χ²)    : %.4f  (df = %d)\n", ms$LRT_stat, ms$df_diff))
cat(sprintf("  p-value          : %.4g\n", ms$p_value))

# Imprimir tablas
cat("\n══ Coeficientes reparametrizados (modelo espacio nulo) ══\n")
print(kable(tab_alpha, format = "simple", align = "lrrllr"))

cat("\n══ Coeficientes originales θ = (β₀, β₁, γ₁, ...) ══\n")
print(kable(tab_theta, format = "simple", align = "lrrllr"))

cat("\n══ Modelo constante (referencia) ══\n")
print(kable(tab_const, format = "simple", align = "lrrllr"))

# En vez de save_kable(), usar writeLines directamente
html_out <- tab_theta %>%
  kable("html", caption = "Original parameters θ — spline VE model (w14, df3, τ₀=15)",
        escape = FALSE) %>%
  kable_styling(
    full_width        = FALSE,
    bootstrap_options = c("striped", "hover", "condensed")
  ) %>%
  as.character()

# Envolver en HTML completo con Bootstrap CDN
html_page <- paste0(
  '<!DOCTYPE html><html><head><meta charset="utf-8">',
  '<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css">',
  '<style>body{padding:30px; font-family: Georgia, serif;}</style>',
  '</head><body>',
  html_out,
  '</body></html>'
)

writeLines(html_page, "coef_table.html")
cat("Tabla guardada en coef_table.html\n")

# ============================================================
# SPLINES: crudas · ponderadas · suma
# ============================================================
res_sel   <- results[["w14_df3_tau_15"]]
knots_    <- res_sel$knots
pc_       <- res_sel$proj_coef
th_       <- res_sel$theta_hat
tau_star_ <- 15
sdf_      <- 3

# Grilla densa
tau_g <- seq(0, tau_star_, by = 0.01)

# Base cruda bs()
B_raw <- make_bs_basis(tau_g, spline_df = sdf_,
                       tau_star_m = tau_star_, knots = knots_)$B

# Base no lineal residualizada
B_nl  <- eval_nonlinear_basis(tau_g, knots = knots_,
                               spline_df = sdf_,
                               tau_star_m = tau_star_,
                               proj_coef = pc_)

# Extraer gamma_j de theta
gamma_names <- grep("^gamma", names(th_), value = TRUE)
gamma_vals  <- th_[gamma_names]

# Para cada spline j: cruda, nl, ponderada
n_sp <- length(gamma_vals)

df_sp <- do.call(rbind, lapply(seq_len(n_sp), function(j) {
  data.frame(
    tau     = tau_g,
    spline  = paste0("B[", j, "](τ)"),
    raw     = B_raw[, j],
    nl      = B_nl[, j],
    wtd     = gamma_vals[j] * B_nl[, j]
  )
}))

# Suma de bases ponderadas + parte lineal
beta0_ <- th_["beta0"]
beta1_ <- th_["beta1"]
eta_g  <- beta0_ + beta1_ * tau_g + B_nl %*% gamma_vals
# Poner a 0 para tau >= tau_star (por construcción)
eta_g[tau_g >= tau_star_] <- 0

df_sum <- data.frame(tau = tau_g,
                     logHR = as.numeric(eta_g),
                     VE    = 1 - exp(as.numeric(eta_g)))

# ── Panel A: splines crudas ──
p_raw <- df_sp %>%
  ggplot(aes(tau, raw, color = spline)) +
  geom_line(linewidth = 0.85) +
  geom_vline(xintercept = knots_, linetype = "dotted", color = "grey60", linewidth = 0.4) +
  labs(title = "A  ·  Raw B-spline basis functions",
       x = "τ (months)", y = "B(τ)", color = NULL) +
  theme_classic(base_size = 11, base_family = "serif") +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

# ── Panel B: bases no lineales residualizadas ──
p_nl <- df_sp %>%
  ggplot(aes(tau, nl, color = spline)) +
  geom_hline(yintercept = 0, color = "grey80") +
  geom_line(linewidth = 0.85) +
  geom_vline(xintercept = knots_, linetype = "dotted", color = "grey60", linewidth = 0.4) +
  labs(title = "B  ·  Non-linear residualized basis  B̃(τ)",
       x = "τ (months)", y = "B̃(τ)", color = NULL) +
  theme_classic(base_size = 11, base_family = "serif") +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

# ── Panel C: bases ponderadas ──
p_wtd <- df_sp %>%
  ggplot(aes(tau, wtd, color = spline)) +
  geom_hline(yintercept = 0, color = "grey80") +
  geom_line(linewidth = 0.85) +
  geom_vline(xintercept = knots_, linetype = "dotted", color = "grey60", linewidth = 0.4) +
  labs(title = "C  ·  Weighted non-linear components  γⱼ·B̃ⱼ(τ)",
       x = "τ (months)", y = "γ·B̃(τ)", color = NULL) +
  theme_classic(base_size = 11, base_family = "serif") +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

# ── Panel D: suma = log-HR(τ) ──
p_sum <- ggplot(df_sum, aes(tau, logHR)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(color = "#1B4F72", linewidth = 1.1) +
  geom_vline(xintercept = tau_star_, linetype = "dotted",
             color = "grey50", linewidth = 0.4) +
  labs(title = "D  ·  Fitted log-HR(τ)  =  β₀ + β₁τ + Σγⱼ·B̃ⱼ(τ)",
       x = "τ (months)", y = "log HR(τ)") +
  theme_classic(base_size = 11, base_family = "serif") +
  theme(plot.title = element_text(face = "bold"))

# Combinar los 4 paneles
library(patchwork)
p_splines <- (p_raw | p_nl) / (p_wtd | p_sum) +
  plot_annotation(
    title   = "Spline decomposition — w14 · df=3 · τ₀=15",
    theme   = theme(plot.title = element_text(face = "bold", size = 13,
                                               family = "serif"))
  )

print(p_splines)
ggsave("spline_decomposition.pdf", p_splines,
       width = 10, height = 7, device = cairo_pdf)







# ============================================================
# DIAGNÓSTICOS ADICIONALES
# ============================================================

res_sel  <- results[["w14_df3_tau_15"]]
df_model <- res_sel$df_model
fit_ns   <- res_sel$fit_ns

# ── 4.1  Schoenfeld residuals (test PH) ──
# Solo aplica a los z-parámetros; el RUN está como cluster, no como covariate
ph_test <- cox.zph(fit_ns, transform = "km")
print(ph_test)        # tabla con p-values por variable y global
plot(ph_test)         # gráfico por variable — importante revisar

# ── 4.4  Número de eventos y persona-tiempo por mes (tabla) ──
cat("\n── Soporte observado por mes (arm inmunizado) ──\n")
print(
  res_sel$support_summary %>%
    mutate(
      rate_per_1000d = round(events / person_time_days * 1000, 3)
    ) %>%
    select(month_bin, n_rows, person_time_days, events, rate_per_1000d) %>%
    kable("simple",
          col.names = c("Month bin", "Rows", "Person-days", "Events", "Rate/1000d"))
)

# ── 4.5  Comparación log-verosimilitud todos los modelos ──
cat("\n── Resumen de todos los modelos ajustados ──\n")
print(
  all_model_summary %>%
    arrange(spline_df, tau_0) %>%
    mutate(across(where(is.numeric), ~ round(., 3))) %>%
    kable("simple")
)