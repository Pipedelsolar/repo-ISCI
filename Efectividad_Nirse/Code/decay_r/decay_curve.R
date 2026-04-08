library(survival)
library(splines)
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)

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

make_cubic_basis <- function(x, x_ref, spline_df, tau_star_m) {
  # Base cúbica cruda
  B_all <- bs(
    x,
    df = spline_df,
    degree = 3,
    Boundary.knots = c(0, tau_star_m),
    intercept = FALSE
  )

  B_ref <- bs(
    x_ref,
    df = spline_df,
    degree = 3,
    Boundary.knots = c(0, tau_star_m),
    intercept = FALSE
  )

  list(B_all = B_all, B_ref = B_ref)
}

fit_waning_model_v_past <- function(df_prepared, split_width_days, spline_df, tau_grid_by = 0.05) {

  df <- df_prepared

  tau_imm <- df$tau_eval_m[df$inmunizado == 1]
  tau_imm <- tau_imm[is.finite(tau_imm)]

  if (length(unique(tau_imm)) <= spline_df) {
    stop(paste0(
      "No hay suficiente soporte único en tau para spline_df=", spline_df,
      " con malla ", split_width_days, " días."
    ))
  }

  ns_template <- ns(
    tau_imm,
    df = spline_df,
    Boundary.knots = range(tau_imm, na.rm = TRUE)
  )

  knots  <- attr(ns_template, "knots")
  bknots <- attr(ns_template, "Boundary.knots")

  B_all <- ns(
    df$tau_eval_m,
    knots = knots,
    Boundary.knots = bknots,
    intercept = FALSE
  )

  B_ref <- as.numeric(predict(ns_template, newx = 0))
  B_ctr <- sweep(B_all, 2, B_ref, FUN = "-")

  df$imm_tau <- df$inmunizado * df$tau_eval_m

  for (j in seq_len(ncol(B_ctr))) {
    df[[paste0("imm_ns", j)]] <- df$inmunizado * B_ctr[, j]
  }

  rhs_ns <- paste(
  c(
    "inmunizado",
    "imm_tau",
    paste0("imm_ns", seq_len(ncol(B_ctr))),
    "strata(Group)",
    "cluster(RUN)"
  ),
  collapse = " + "
  )

  f_ns <- as.formula(
    paste0("Surv(start, stop, event_vrs) ~ ", rhs_ns)
  )

  fit_ns <- coxph(
    formula = f_ns,
    data    = df,
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

  tau_grid_m <- seq(
    0,
    ceiling(max(df$tau_eval_m[df$inmunizado == 1], na.rm = TRUE)),
    by = tau_grid_by
  )

  B_grid <- ns(
    tau_grid_m,
    knots = knots,
    Boundary.knots = bknots,
    intercept = FALSE
  )

  B_grid_ctr <- sweep(B_grid, 2, B_ref, FUN = "-")

  Xg <- cbind(
    inmunizado = 1,
    imm_tau = tau_grid_m,
    B_grid_ctr
  )
  colnames(Xg) <- c("inmunizado", "imm_tau", paste0("imm_ns", seq_len(ncol(B_grid_ctr))))

  beta_full <- coef(fit_ns)
  keep_beta <- !is.na(beta_full)

  beta <- beta_full[keep_beta]
  V    <- vcov(fit_ns)[keep_beta, keep_beta, drop = FALSE]
  Xg   <- Xg[, names(beta), drop = FALSE]

  eta    <- as.vector(Xg %*% beta)
  se_eta <- sqrt(rowSums((Xg %*% V) * Xg))

  HR    <- exp(eta)
  HR_lo <- exp(eta - 1.96 * se_eta)
  HR_hi <- exp(eta + 1.96 * se_eta)

  VE    <- 1 - HR
  VE_lo <- 1 - HR_hi
  VE_hi <- 1 - HR_lo

  ve_curve <- data.frame(
    split_width_days = split_width_days,
    spline_df        = spline_df,
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

  model_summary <- data.frame(
    split_width_days = split_width_days,
    spline_df        = spline_df,
    n_rows           = nrow(df),
    n_subjects       = dplyr::n_distinct(df$RUN),
    n_events         = sum(df$event_vrs),
    loglik_ns        = fit_ns$loglik[2],
    loglik_const     = fit_const$loglik[2],
    LRT_stat         = 2 * (fit_ns$loglik[2] - fit_const$loglik[2]),
    df_diff          = df_diff <- sum(!is.na(coef(fit_ns))) - sum(!is.na(coef(fit_const))),
    p_value          = pchisq(
      2 * (fit_ns$loglik[2] - fit_const$loglik[2]),
      df = sum(!is.na(coef(fit_ns))) - sum(!is.na(coef(fit_const))),
      lower.tail = FALSE
    )
  )

  list(
    df_model         = df,
    fit_ns           = fit_ns,
    fit_const        = fit_const,
    ve_curve         = ve_curve,
    support_summary  = support_summary,
    model_summary    = model_summary
  )
}

fit_waning_model_v_past_2 <- function(df_prepared, split_width_days, spline_df, tau_grid_by = 0.05, tau_star_m = 15, use_tail_constraint = TRUE) {

  df <- df_prepared

  tau_imm <- df$tau_eval_m[df$inmunizado == 1]
  tau_imm <- tau_imm[is.finite(tau_imm)]

  if (length(unique(tau_imm)) <= spline_df) {
    stop(paste0(
      "No hay suficiente soporte único en tau para spline_df=", spline_df,
      " con malla ", split_width_days, " días."
    ))
  }

  # ==========================================================
  # BASE SPLINE CRUDA
  # Si imponemos 2 restricciones en la cola, conviene más df.
  # ==========================================================
  ns_template <- ns(
    tau_imm,
    df = spline_df,
    Boundary.knots = c(0, tau_star_m)
  )

  knots  <- attr(ns_template, "knots")
  bknots <- attr(ns_template, "Boundary.knots")

  B_all <- ns(
    df$tau_eval_m,
    knots = knots,
    Boundary.knots = bknots,
    intercept = FALSE
  )

  # ==========================================================
  # SI use_tail_constraint=TRUE:
  # Construimos base C_j(tau) tal que:
  # C_j(tau_star)=0 y C_j'(tau_star)=0
  # ==========================================================
  if (use_tail_constraint) {
    eps <- 1e-6

    B_star  <- ns(tau_star_m, knots = knots, Boundary.knots = bknots, intercept = FALSE)
    B_plus  <- ns(tau_star_m + eps, knots = knots, Boundary.knots = bknots, intercept = FALSE)
    B_minus <- ns(tau_star_m - eps, knots = knots, Boundary.knots = bknots, intercept = FALSE)

    dB_star <- (B_plus - B_minus) / (2 * eps)

    B_ctr <- B_all
    for (j in seq_len(ncol(B_all))) {
      B_ctr[, j] <- B_all[, j] - B_star[1, j] - dB_star[1, j] * (df$tau_eval_m - tau_star_m)
    }

    # por estabilidad numérica, dejar exactamente 0 en la cola
    B_ctr[df$tau_eval_m >= tau_star_m, ] <- 0
  } else {
    # versión anterior: solo centrar en tau=0
    B_ref <- as.numeric(predict(ns_template, newx = 0))
    B_ctr <- sweep(B_all, 2, B_ref, FUN = "-")
  }

  # ==========================================================
  # INTERACCIONES
  # Si hay constraint de cola, NO meter imm_tau aparte.
  # Si no hay constraint, puedes seguir con imm_tau.
  # ==========================================================
  if (!use_tail_constraint) {
    df$imm_tau <- df$inmunizado * df$tau_eval_m
  }

  for (j in seq_len(ncol(B_ctr))) {
    df[[paste0("imm_ns", j)]] <- df$inmunizado * B_ctr[, j]
  }

  rhs_terms <- c()

  if (!use_tail_constraint) {
    rhs_terms <- c(rhs_terms, "inmunizado", "imm_tau")
  }

  rhs_terms <- c(
    rhs_terms,
    paste0("imm_ns", seq_len(ncol(B_ctr))),
    "strata(Group)",
    "cluster(RUN)"
  )

  rhs_ns <- paste(rhs_terms, collapse = " + ")

  f_ns <- as.formula(
    paste0("Surv(start, stop, event_vrs) ~ ", rhs_ns)
  )

  fit_ns <- coxph(
    formula = f_ns,
    data    = df,
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

  tau_grid_m <- seq(
    0,
    ceiling(max(df$tau_eval_m[df$inmunizado == 1], na.rm = TRUE)),
    by = tau_grid_by
  )

  B_grid <- ns(
    tau_grid_m,
    knots = knots,
    Boundary.knots = bknots,
    intercept = FALSE
  )

  if (use_tail_constraint) {
    eps <- 1e-6
    B_star  <- ns(tau_star_m, knots = knots, Boundary.knots = bknots, intercept = FALSE)
    B_plus  <- ns(tau_star_m + eps, knots = knots, Boundary.knots = bknots, intercept = FALSE)
    B_minus <- ns(tau_star_m - eps, knots = knots, Boundary.knots = bknots, intercept = FALSE)
    dB_star <- (B_plus - B_minus) / (2 * eps)

    B_grid_ctr <- B_grid
    for (j in seq_len(ncol(B_grid))) {
      B_grid_ctr[, j] <- B_grid[, j] - B_star[1, j] - dB_star[1, j] * (tau_grid_m - tau_star_m)
    }

    B_grid_ctr[tau_grid_m >= tau_star_m, ] <- 0

    Xg <- B_grid_ctr
    colnames(Xg) <- paste0("imm_ns", seq_len(ncol(B_grid_ctr)))
  } else {
    B_ref <- as.numeric(predict(ns_template, newx = 0))
    B_grid_ctr <- sweep(B_grid, 2, B_ref, FUN = "-")

    Xg <- cbind(
      inmunizado = 1,
      imm_tau = tau_grid_m,
      B_grid_ctr
    )
    colnames(Xg) <- c("inmunizado", "imm_tau", paste0("imm_ns", seq_len(ncol(B_grid_ctr))))
  }

  beta_full <- coef(fit_ns)
  keep_beta <- !is.na(beta_full)

  beta <- beta_full[keep_beta]
  V    <- vcov(fit_ns)[keep_beta, keep_beta, drop = FALSE]
  Xg   <- Xg[, names(beta), drop = FALSE]

  eta    <- as.vector(Xg %*% beta)
  se_eta <- sqrt(rowSums((Xg %*% V) * Xg))

  HR    <- exp(eta)
  HR_lo <- exp(eta - 1.96 * se_eta)
  HR_hi <- exp(eta + 1.96 * se_eta)

  VE    <- 1 - HR
  VE_lo <- 1 - HR_hi
  VE_hi <- 1 - HR_lo

  ve_curve <- data.frame(
    split_width_days = split_width_days,
    spline_df        = spline_df,
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

  model_summary <- data.frame(
    split_width_days = split_width_days,
    spline_df        = spline_df,
    n_rows           = nrow(df),
    n_subjects       = dplyr::n_distinct(df$RUN),
    n_events         = sum(df$event_vrs),
    loglik_ns        = fit_ns$loglik[2],
    loglik_const     = fit_const$loglik[2],
    LRT_stat         = 2 * (fit_ns$loglik[2] - fit_const$loglik[2]),
    df_diff          = sum(!is.na(coef(fit_ns))) - sum(!is.na(coef(fit_const))),
    p_value          = pchisq(
      2 * (fit_ns$loglik[2] - fit_const$loglik[2]),
      df = sum(!is.na(coef(fit_ns))) - sum(!is.na(coef(fit_const))),
      lower.tail = FALSE
    )
  )

  list(
    df_model         = df,
    fit_ns           = fit_ns,
    fit_const        = fit_const,
    ve_curve         = ve_curve,
    support_summary  = support_summary,
    model_summary    = model_summary
  )
}


fit_waning_model <- function(df_prepared, split_width_days, spline_df,
                             tau_grid_by = 0.05,
                             tau_star_m = 15,
                             use_tail_constraint = TRUE) {

  df <- df_prepared

  tau_imm <- df$tau_eval_m[df$inmunizado == 1]
  tau_imm <- tau_imm[is.finite(tau_imm)]

  if (length(unique(tau_imm)) <= spline_df) {
    stop(paste0(
      "No hay suficiente soporte único en tau para spline_df=", spline_df,
      " con malla ", split_width_days, " días."
    ))
  }

  # ==========================================================
  # BASE CÚBICA CRUDA
  # ==========================================================
  B_all <- bs(
    df$tau_eval_m,
    df = spline_df,
    degree = 3,
    Boundary.knots = c(0, tau_star_m),
    intercept = FALSE
  )

  # ==========================================================
  # TAIL CONSTRAINT:
  # C_j(tau_star)=0 y C_j'(tau_star)=0
  # ==========================================================
  if (use_tail_constraint) {
    eps <- 1e-6

    B_star  <- bs(tau_star_m, knots = attr(B_all, "knots"),
                  degree = 3,
                  Boundary.knots = c(0, tau_star_m),
                  intercept = FALSE)

    B_plus  <- bs(tau_star_m + eps, knots = attr(B_all, "knots"),
                  degree = 3,
                  Boundary.knots = c(0, tau_star_m),
                  intercept = FALSE)

    B_minus <- bs(tau_star_m - eps, knots = attr(B_all, "knots"),
                  degree = 3,
                  Boundary.knots = c(0, tau_star_m),
                  intercept = FALSE)

    dB_star <- (B_plus - B_minus) / (2 * eps)

    B_ctr <- B_all
    for (j in seq_len(ncol(B_all))) {
      B_ctr[, j] <- B_all[, j] - B_star[1, j] - dB_star[1, j] * (df$tau_eval_m - tau_star_m)
    }

    B_ctr[df$tau_eval_m >= tau_star_m, ] <- 0
  } else {
    B_ref <- bs(
      0,
      knots = attr(B_all, "knots"),
      degree = 3,
      Boundary.knots = c(0, tau_star_m),
      intercept = FALSE
    )
    B_ctr <- sweep(B_all, 2, B_ref, FUN = "-")
  }

  # ==========================================================
  # INTERACCIONES
  # En modo constrained NO metemos imm_tau libre.
  # ==========================================================
  if (!use_tail_constraint) {
    df$imm_tau <- df$inmunizado * df$tau_eval_m
  }

  for (j in seq_len(ncol(B_ctr))) {
    df[[paste0("imm_bs", j)]] <- df$inmunizado * B_ctr[, j]
  }

  rhs_terms <- c()

  if (!use_tail_constraint) {
    rhs_terms <- c(rhs_terms, "inmunizado", "imm_tau")
  }

  rhs_terms <- c(
    rhs_terms,
    paste0("imm_bs", seq_len(ncol(B_ctr))),
    "strata(Group)",
    "cluster(RUN)"
  )

  rhs_ns <- paste(rhs_terms, collapse = " + ")

  f_ns <- as.formula(
    paste0("Surv(start, stop, event_vrs) ~ ", rhs_ns)
  )

  fit_ns <- coxph(
    formula = f_ns,
    data    = df,
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

  tau_grid_m <- seq(
    0,
    ceiling(max(df$tau_eval_m[df$inmunizado == 1], na.rm = TRUE)),
    by = tau_grid_by
  )

  B_grid <- bs(
    tau_grid_m,
    knots = attr(B_all, "knots"),
    degree = 3,
    Boundary.knots = c(0, tau_star_m),
    intercept = FALSE
  )

  if (use_tail_constraint) {
    eps <- 1e-6

    B_star  <- bs(tau_star_m, knots = attr(B_all, "knots"),
                  degree = 3,
                  Boundary.knots = c(0, tau_star_m),
                  intercept = FALSE)

    B_plus  <- bs(tau_star_m + eps, knots = attr(B_all, "knots"),
                  degree = 3,
                  Boundary.knots = c(0, tau_star_m),
                  intercept = FALSE)

    B_minus <- bs(tau_star_m - eps, knots = attr(B_all, "knots"),
                  degree = 3,
                  Boundary.knots = c(0, tau_star_m),
                  intercept = FALSE)

    dB_star <- (B_plus - B_minus) / (2 * eps)

    B_grid_ctr <- B_grid
    for (j in seq_len(ncol(B_grid))) {
      B_grid_ctr[, j] <- B_grid[, j] - B_star[1, j] - dB_star[1, j] * (tau_grid_m - tau_star_m)
    }

    B_grid_ctr[tau_grid_m >= tau_star_m, ] <- 0

    Xg <- B_grid_ctr
    colnames(Xg) <- paste0("imm_bs", seq_len(ncol(B_grid_ctr)))
  } else {
    B_ref <- bs(
      0,
      knots = attr(B_all, "knots"),
      degree = 3,
      Boundary.knots = c(0, tau_star_m),
      intercept = FALSE
    )

    B_grid_ctr <- sweep(B_grid, 2, B_ref, FUN = "-")

    Xg <- cbind(
      inmunizado = 1,
      imm_tau = tau_grid_m,
      B_grid_ctr
    )
    colnames(Xg) <- c("inmunizado", "imm_tau", paste0("imm_bs", seq_len(ncol(B_grid_ctr))))
  }

  beta_full <- coef(fit_ns)
  keep_beta <- !is.na(beta_full)

  beta <- beta_full[keep_beta]
  V    <- vcov(fit_ns)[keep_beta, keep_beta, drop = FALSE]
  Xg   <- Xg[, names(beta), drop = FALSE]

  eta    <- as.vector(Xg %*% beta)
  se_eta <- sqrt(rowSums((Xg %*% V) * Xg))

  HR    <- exp(eta)
  HR_lo <- exp(eta - 1.96 * se_eta)
  HR_hi <- exp(eta + 1.96 * se_eta)

  VE    <- 1 - HR
  VE_lo <- 1 - HR_hi
  VE_hi <- 1 - HR_lo

  ve_curve <- data.frame(
    split_width_days = split_width_days,
    spline_df        = spline_df,
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

  model_summary <- data.frame(
    split_width_days = split_width_days,
    spline_df        = spline_df,
    n_rows           = nrow(df),
    n_subjects       = dplyr::n_distinct(df$RUN),
    n_events         = sum(df$event_vrs),
    loglik_ns        = fit_ns$loglik[2],
    loglik_const     = fit_const$loglik[2],
    LRT_stat         = 2 * (fit_ns$loglik[2] - fit_const$loglik[2]),
    df_diff          = sum(!is.na(coef(fit_ns))) - sum(!is.na(coef(fit_const))),
    p_value          = pchisq(
      2 * (fit_ns$loglik[2] - fit_const$loglik[2]),
      df = sum(!is.na(coef(fit_ns))) - sum(!is.na(coef(fit_const))),
      lower.tail = FALSE
    )
  )

  list(
    df_model         = df,
    fit_ns           = fit_ns,
    fit_const        = fit_const,
    ve_curve         = ve_curve,
    support_summary  = support_summary,
    model_summary    = model_summary
  )
}

# ============================================================
# 3) PREPARAR LAS 3 MALLAS UNA SOLA VEZ
# ============================================================

split_grid  <- c(14) #c(7, 14, 30)
spline_grid <- c(3,4,5) #c(3, 4, 5)

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
    key <- paste0("w", w, "_df", sdf)
    message("Ajustando modelo: malla=", w, " días ; spline_df=", sdf)

    results[[key]] <- fit_waning_model(
      df_prepared      = prepared_by_width[[paste0("w", w)]],
      split_width_days = w,
      spline_df        = sdf,
      tau_grid_by      = 0.05,
      tau_star_m       = 15,
      use_tail_constraint = TRUE
    )
  }
}

# ============================================================
# 5) TABLAS CONSOLIDADAS
# ============================================================
all_curves <- bind_rows(lapply(results, function(x) x$ve_curve))
all_support <- bind_rows(lapply(results, function(x) x$support_summary))
all_model_summary <- bind_rows(lapply(results, function(x) x$model_summary))

print(all_model_summary)

# Si quieres exportarlas:
# write.csv(all_curves, "../Output/all_ve_curves.csv", row.names = FALSE)
# write.csv(all_support, "../Output/all_support_summary.csv", row.names = FALSE)
# write.csv(all_model_summary, "../Output/all_model_summary.csv", row.names = FALSE)

# ============================================================
# 6) UNA FIGURA POR spline_df, CON LAS 3 MALLAS + IC
# ============================================================
# plot_list <- list()

# for (sdf in spline_grid) {
#   p <- all_curves %>%
#     filter(spline_df == sdf) %>%
#     mutate(split_width_days = factor(split_width_days, levels = c(14))) %>%
#     ggplot(aes(x = tau_m, y = VE, color = split_width_days, fill = split_width_days)) +
#     geom_ribbon(aes(ymin = VE_lo, ymax = VE_hi), alpha = 0.10, linewidth = 0) +
#     geom_line(linewidth = 1) +
#     geom_hline(yintercept = 0, linetype = 2) +
#     scale_y_continuous(
#       limits = c(0, 1),
#       breaks = seq(0, 1, by = 0.1)
#     ) +
#     coord_cartesian(xlim = c(0, 15), ylim = c(0, 1)) +
#     labs(
#       x = "Months since immunization",
#       y = "VE(t) = 1 - HR(t)",
#       color = "Split width (days)",
#       fill  = "Split width (days)",
#       title = paste0("Spline-smoothed waning curve | ns(df = ", sdf, ")")
#     ) +
#     theme_bw() +
#     theme(
#       aspect.ratio = 0.6
#     )

#   plot_list[[paste0("df", sdf)]] <- p
#   print(p)
# }


# Si quieres guardar automáticamente:
# dir.create("../Output/waning_plots", showWarnings = FALSE)
# for (nm in names(plot_list)) {
#   ggsave(
#     filename = paste0("../Output/waning_plots/", nm, ".png"),
#     plot = plot_list[[nm]],
#     width = 9,
#     height = 5.5,
#     dpi = 300
#   )
# }

# ============================================================
# 7) OPCIONAL: FIGURA FACETEADA 3x1
# ============================================================
breaks_01 <- function(x) {
  lo <- floor((min(x, na.rm = TRUE) + 1e-8) * 10) / 10
  hi <- ceiling((max(x, na.rm = TRUE) - 1e-8) * 10) / 10
  seq(lo, hi, by = 0.1)
}

labels_01 <- function(x) {
  x[abs(x) < 1e-8] <- 0
  sprintf("%.1f", x)
}

p_facet <- all_curves %>%
  mutate(
    split_width_days = factor(split_width_days, levels = c(14)),
    spline_df = factor(spline_df, levels = c(3, 4, 5))
  ) %>%
  ggplot(aes(x = tau_m, y = VE, color = split_width_days, fill = split_width_days)) +
  geom_ribbon(aes(ymin = VE_lo, ymax = VE_hi), alpha = 0.08, linewidth = 0) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~ spline_df, ncol = 1, scales = "free_y") +
  scale_y_continuous(
    breaks = breaks_01,
    labels = labels_01
  ) +
  coord_cartesian(xlim = c(0, 15)) +
  labs(
    x = "Months since immunization",
    y = "VE(t) = 1 - HR(t)",
    color = "Split width (days)",
    fill  = "Split width (days)",
    title = "Sensitivity to split width and spline flexibility"
  ) +
  theme_bw(base_size = 12) +
  theme(
    aspect.ratio = 0.45,
    legend.position = "bottom"
  )

print(p_facet)


# coef(results[["w14_df3"]]$fit_ns)

# summary(results[["w14_df3"]]$fit_ns)
    
################### SPLINES GRAPHICS ################


# Elegir modelo
fit_obj <- results[["w14_df3"]]
tau_star_m <- 15
spline_df <- 3

# Reconstruir soporte observado
df_model <- fit_obj$df_model
tau_imm <- df_model$tau_eval_m[df_model$inmunizado == 1]
tau_imm <- tau_imm[is.finite(tau_imm)]

# Reconstruir la base cruda cúbica tal como en fit_waning_model
B_tmp <- bs(
  tau_imm,
  df = spline_df,
  degree = 3,
  Boundary.knots = c(0, tau_star_m),
  intercept = FALSE
)

knots_internal <- attr(B_tmp, "knots")

tau_grid <- seq(0, tau_star_m + 3, by = 0.05)

B_grid <- bs(
  tau_grid,
  knots = knots_internal,
  degree = 3,
  Boundary.knots = c(0, tau_star_m),
  intercept = FALSE
)

eps <- 1e-6
B_star  <- bs(tau_star_m, knots = knots_internal, degree = 3,
              Boundary.knots = c(0, tau_star_m), intercept = FALSE)
B_plus  <- bs(tau_star_m + eps, knots = knots_internal, degree = 3,
              Boundary.knots = c(0, tau_star_m), intercept = FALSE)
B_minus <- bs(tau_star_m - eps, knots = knots_internal, degree = 3,
              Boundary.knots = c(0, tau_star_m), intercept = FALSE)

dB_star <- (B_plus - B_minus) / (2 * eps)

C_grid <- B_grid
for (j in seq_len(ncol(B_grid))) {
  C_grid[, j] <- B_grid[, j] - B_star[1, j] - dB_star[1, j] * (tau_grid - tau_star_m)
}
C_grid[tau_grid >= tau_star_m, ] <- 0

C_df <- as.data.frame(C_grid)
C_df$tau_m <- tau_grid

library(tidyr)
library(ggplot2)
C_long <- pivot_longer(C_df, cols = -tau_m, names_to = "basis", values_to = "value")

ggplot(C_long, aes(x = tau_m, y = value, color = basis)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = tau_star_m, linetype = 2) +
  theme_bw() +
  labs(
    x = "Months since immunization",
    y = "Basis value",
    title = "Tail-constrained cubic spline basis functions"
  )


 ################# 2 ###################3

beta <- coef(fit_obj$fit_ns)
beta <- beta[!is.na(beta)]

# columnas del modelo que corresponden a la base
basis_names <- names(beta)

# asegurarnos de que coincidan con las columnas de C_grid
colnames(C_grid) <- basis_names

# contribuciones individuales
contr_grid <- sweep(C_grid, 2, beta, `*`)
contr_df <- as.data.frame(contr_grid)
contr_df$tau_m <- tau_grid

contr_long <- pivot_longer(contr_df, cols = -tau_m, names_to = "term", values_to = "contribution")

ggplot(contr_long, aes(x = tau_m, y = contribution, color = term)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = tau_star_m, linetype = 2) +
  theme_bw() +
  labs(
    x = "Months since immunization",
    y = "Contribution to log-HR",
    title = "Term-by-term contributions to log-HR"
  )
###################### 3 sumados #######################


eta_grid <- rowSums(contr_grid)

plot_df <- data.frame(
  tau_m = tau_grid,
  logHR = eta_grid,
  HR = exp(eta_grid),
  VE = 1 - exp(eta_grid)
)

ggplot(plot_df, aes(x = tau_m, y = logHR)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = tau_star_m, linetype = 2) +
  theme_bw() +
  labs(
    x = "Months since immunization",
    y = "log-HR(tau)",
    title = "Estimated log-HR curve"
  )


# ============================================================
# BOOTSTRAP HELPERS
# ============================================================

make_boot_sample_by_group <- function(df) {
  groups <- unique(as.character(df$Group))
  sampled_groups <- sample(groups, size = length(groups), replace = TRUE)

  out_list <- vector("list", length(sampled_groups))

  for (i in seq_along(sampled_groups)) {
    g <- sampled_groups[i]

    tmp <- df[df$Group == g, , drop = FALSE]

    # Renombrar Group para que cada réplica bootstrap sea un estrato distinto
    tmp$Group <- paste0("bootG_", i)

    # Renombrar RUN para que si el mismo grupo sale 2 veces,
    # los individuos no colapsen en el cluster robusto
    tmp$RUN <- paste0(tmp$RUN, "_boot", i)

    out_list[[i]] <- tmp
  }

  dplyr::bind_rows(out_list)
}

bootstrap_waning_model <- function(df_prepared,
                                   split_width_days,
                                   spline_df,
                                   tau_grid_by = 0.05,
                                   tau_star_m = 15,
                                   use_tail_constraint = TRUE,
                                   R = 200,
                                   seed = 123,
                                   verbose = TRUE) {

  set.seed(seed)

  # Ajuste original
  fit0 <- fit_waning_model(
    df_prepared = df_prepared,
    split_width_days = split_width_days,
    spline_df = spline_df,
    tau_grid_by = tau_grid_by,
    tau_star_m = tau_star_m,
    use_tail_constraint = use_tail_constraint
  )

  tau_grid <- fit0$ve_curve$tau_m
  B <- matrix(NA_real_, nrow = R, ncol = length(tau_grid))

  for (r in seq_len(R)) {
    if (verbose && (r %% 10 == 0 || r == 1)) {
      message("Bootstrap replicate ", r, " / ", R)
    }

    df_b <- make_boot_sample_by_group(df_prepared)

    fit_b <- tryCatch(
      fit_waning_model(
        df_prepared = df_b,
        split_width_days = split_width_days,
        spline_df = spline_df,
        tau_grid_by = tau_grid_by,
        tau_star_m = tau_star_m,
        use_tail_constraint = use_tail_constraint
      ),
      error = function(e) NULL
    )

    if (!is.null(fit_b)) {
      # guardamos logHR
      B[r, ] <- fit_b$ve_curve$logHR
    }
  }

  # cuántas réplicas funcionaron
  ok <- apply(B, 1, function(x) all(is.finite(x)))
  if (verbose) {
    message("Successful bootstrap replicates: ", sum(ok), " / ", R)
  }

  B_ok <- B[ok, , drop = FALSE]

  if (nrow(B_ok) < 20) {
    stop("Muy pocas réplicas bootstrap exitosas para construir IC.")
  }

  q_logHR <- apply(B_ok, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

  out_curve <- fit0$ve_curve
  out_curve$logHR_boot_lo <- q_logHR[1, ]
  out_curve$logHR_boot_hi <- q_logHR[2, ]

  out_curve$HR_boot_lo <- exp(out_curve$logHR_boot_lo)
  out_curve$HR_boot_hi <- exp(out_curve$logHR_boot_hi)

  out_curve$VE_boot_lo <- 1 - out_curve$HR_boot_hi
  out_curve$VE_boot_hi <- 1 - out_curve$HR_boot_lo

  list(
    fit_original = fit0,
    boot_logHR = B_ok,
    ve_curve_boot = out_curve,
    n_success = nrow(B_ok),
    n_total = R
  )
}


boot_res <- bootstrap_waning_model(
  df_prepared = prepared_by_width[["w14"]],
  split_width_days = 14,
  spline_df = 5,
  tau_grid_by = 0.05,
  tau_star_m = 15,
  use_tail_constraint = TRUE,
  R = 100,
  seed = 123,
  verbose = TRUE
)


curve_boot <- boot_res$ve_curve_boot

p_boot <- ggplot(curve_boot, aes(x = tau_m, y = VE)) +
  geom_ribbon(aes(ymin = VE_boot_lo, ymax = VE_boot_hi), alpha = 0.20) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(
    x = "Months since immunization",
    y = "VE(t) = 1 - HR(t)",
    title = paste0(
      "Tail-constrained Cox waning curve with bootstrap CI (",
      boot_res$n_success, "/", boot_res$n_total, " successful replicates)"
    )
  ) +
  theme_bw()

print(p_boot)
