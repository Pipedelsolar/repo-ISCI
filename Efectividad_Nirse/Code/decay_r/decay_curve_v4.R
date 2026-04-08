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
 # ── Predicción en grilla — versión corregida ──
    tau_grid_m  <- seq(0, tau_star_m, by = tau_grid_by)
    B_grid_nl   <- eval_nonlinear_basis(tau_grid_m, knots, spline_df,
                                        tau_star_m, proj_coef)
    active_grid <- as.numeric(tau_grid_m < tau_star_m)

    Xg_orig <- cbind(
    beta0 = active_grid,
    beta1 = active_grid * tau_grid_m,
    active_grid * B_grid_nl
    )
    colnames(Xg_orig) <- colnames(X_orig)

    # Punto estimado (igual que antes)
    eta <- as.vector(Xg_orig %*% theta_hat)

    # ── SE corregido: calcular en espacio alpha (null-space) ──
    # En vez de propagar V_theta (que tiene rango deficiente cerca de tau_0)
    # usamos directamente V_alpha en el espacio reparametrizado
    Xg_z   <- Xg_orig %*% N_keep                         # mapear a espacio alpha
    se_eta <- sqrt(pmax(0, rowSums((Xg_z %*% V_alpha) * Xg_z)))

    # Forzar se = 0 exactamente en tau_0 (por construcción)
    se_eta[tau_grid_m >= tau_star_m] <- 0

    # También forzar se pequeño muy cerca de tau_0
    # (la restricción elimina incertidumbre en ese punto)
    near_tau0 <- tau_grid_m >= (tau_star_m - 0.5)
    se_eta[near_tau0] <- se_eta[near_tau0] *
    (1 - (tau_grid_m[near_tau0] - (tau_star_m - 0.5)) / 0.5)^2

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

print(all_model_summary)

# ============================================================
# 6) FIGURA SENSIBILIDAD (p_tau_grid)
# ============================================================

labels_01 <- function(x) {
  x[abs(x) < 1e-8] <- 0
  sprintf("%.1f", x)
}

p_tau_grid <- all_curves %>%
  mutate(
    split_width_days = factor(split_width_days, levels = c(14)),
    spline_df = factor(spline_df, levels = c(3, 4, 5, 6),
                       labels = c("ns(df = 3)", "ns(df = 4)",
                                  "ns(df = 5)", "ns(df = 6)")),
    tau_0 = factor(tau_0, levels = c(13, 14, 15),
                   labels = c(expression(tau[0] == 13),
                              expression(tau[0] == 14),
                              expression(tau[0] == 15)))
  ) %>%
  ggplot(aes(x = tau_m, y = VE,
             color = split_width_days, fill = split_width_days)) +
  geom_ribbon(aes(ymin = VE_lo, ymax = VE_hi), alpha = 0.08, linewidth = 0) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_grid(rows = vars(spline_df), cols = vars(tau_0),
             scales = "free_y", labeller = label_parsed) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = labels_01) +
  scale_x_continuous(breaks = seq(0, 15, by = 3)) +
  coord_cartesian(xlim = c(0, 15)) +
  labs(
    x     = "Months since immunization",
    y     = "VE(t) = 1 - HR(t)",
    color = "Split width (days)",
    fill  = "Split width (days)",
    title = "Sensitivity to spline flexibility and tail constraint"
  ) +
  theme_bw(base_size = 11) +
  theme(aspect.ratio = 0.55, legend.position = "bottom")

print(p_tau_grid)


# ============================================================
# 7) FIGURA PRINCIPAL: VE(t) paper-ready
# ============================================================

key_sel  <- "w14_df3_tau_15"
res_sel  <- results[[key_sel]]
curve    <- res_sel$ve_curve
support  <- res_sel$support_summary

support_hist <- support %>%
  mutate(rel_pt = person_time_days / max(person_time_days, na.rm = TRUE) * 0.18)

col_line <- "#1B4F72"
col_band <- "#2E86C1"
col_hist <- "grey82"

p_main <- ggplot() +
  geom_rect(
    data = support_hist,
    aes(xmin = month_bin, xmax = month_bin + 0.9,
        ymin = -0.04, ymax = -0.04 + rel_pt),
    fill = col_hist, color = NA, alpha = 0.9
  ) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.45) +
  geom_ribbon(data = curve,
              aes(x = tau_m, ymin = VE_lo, ymax = VE_hi),
              fill = col_band, alpha = 0.15) +
  geom_line(data = curve,
            aes(x = tau_m, y = VE),
            color = col_line, linewidth = 1.1) +
  geom_vline(xintercept = 15, linetype = "dotted",
             color = "grey50", linewidth = 0.5) +
  annotate("text", x = 15.15, y = 0.92,
           label = expression(tau[0] == 15),
           hjust = 0, size = 3.2, color = "grey40", family = "serif") +
  annotate("text", x = 0.1, y = -0.01,
           label = "Person-time\n(relative)",
           hjust = 0, vjust = 1, size = 2.6,
           color = "grey50", family = "serif") +
  scale_x_continuous(breaks = 0:15, expand = c(0.01, 0)) +
  scale_y_continuous(
    breaks = c(-0.25, 0, 0.25, 0.5, 0.75, 1.0),
    labels = scales::percent_format(accuracy = 1),
    limits = c(-0.27, 1.02)
  ) +
  coord_cartesian(xlim = c(0, 15.6)) +
  labs(
    x        = "Months since immunization",
    y        = "Vaccine effectiveness  VE(t) = 1 − HR(t)",
    title    = "Waning of vaccine effectiveness over time",
    subtitle = "Cox model with constrained cubic spline  ·  split = 14 d, df = 3, τ₀ = 15 mo",
    caption  = "Shaded area: pointwise 95% CI.  Bars: relative person-time at risk (immunized arm)."
  ) +
  theme_classic(base_size = 12, base_family = "serif") +
  theme(
    plot.title         = element_text(face = "bold", size = 20),
    plot.subtitle      = element_text(color = "grey35", size = 9.5),
    plot.caption       = element_text(color = "grey45", size = 8, hjust = 0),
    axis.title         = element_text(size = 10.5),
    axis.text          = element_text(size = 9.5),
    panel.grid.major.y = element_line(color = "grey93", linewidth = 0.35),
    panel.grid.major.x = element_blank(),
    plot.margin        = margin(12, 16, 8, 8)
  )

print(p_main)
ggsave("ve_curve_paper.pdf", p_main, width = 6.5, height = 4.2, device = cairo_pdf)
ggsave("ve_curve_paper.png", p_main, width = 6.5, height = 4.2, dpi = 320)

# ============================================================
# 8) FUNCIÓN GENÉRICA (familias alternativas)
# ============================================================

make_f_spline <- function(tau, tau_star_m, spline_df, knots = NULL) {
  obj <- make_bs_basis(tau, spline_df = spline_df,
                       tau_star_m = tau_star_m, knots = knots)
  list(B = obj$B, knots = obj$knots)
}

make_f_log <- function(tau, tau_star_m, ...) {
  list(B = matrix(log(1 + tau), ncol = 1), knots = NULL)
}

make_f_weibull <- function(tau, tau_star_m, gamma = 2, ...) {
  list(B = matrix(tau^gamma, ncol = 1), knots = NULL)
}

make_f_exponential <- function(tau, tau_star_m, lambda = 0.3, ...) {
  list(B = matrix(exp(-lambda * tau), ncol = 1), knots = NULL)
}

make_f_power <- function(tau, tau_star_m, p = 0.5, ...) {
  list(B = matrix(tau^p, ncol = 1), knots = NULL)
}

fit_waning_generic <- function(df_prepared,
                                split_width_days,
                                tau_star_m   = 15,
                                tau_grid_by  = 0.05,
                                f_basis_fn   = make_f_spline,
                                f_basis_args = list(spline_df = 3),
                                family_label = "spline") {

  df            <- df_prepared
  tau_eval_clip <- pmin(df$tau_eval_m, tau_star_m)
  active        <- as.numeric(df$inmunizado == 1 & df$tau_eval_m < tau_star_m)

  tau_imm <- tau_eval_clip[df$inmunizado == 1 & df$tau_eval_m < tau_star_m]
  tau_imm <- tau_imm[is.finite(tau_imm)]

  f_support <- do.call(f_basis_fn, c(list(tau = tau_imm, tau_star_m = tau_star_m),
                                      f_basis_args))
  knots_out <- f_support$knots

  f_all   <- do.call(f_basis_fn, c(list(tau = tau_eval_clip, tau_star_m = tau_star_m),
                                    c(f_basis_args, list(knots = knots_out))))$B
  n_f     <- ncol(f_all)
  f_names <- paste0("f", seq_len(n_f))

  X_orig <- cbind(
    beta0 = active,
    beta1 = active * tau_eval_clip,
    sweep(f_all, 1, active, `*`)
  )
  colnames(X_orig) <- c("beta0", "beta1", f_names)

  eps         <- 1e-6
  tau_at_star <- tau_star_m - eps

  eval_f_at <- function(t) {
    do.call(f_basis_fn, c(list(tau = t, tau_star_m = tau_star_m),
                           c(f_basis_args, list(knots = knots_out))))$B
  }

  f_star  <- eval_f_at(tau_at_star)
  df_star <- (eval_f_at(tau_at_star + eps) - eval_f_at(tau_at_star - eps)) / (2 * eps)

  A <- rbind(
    c(1, tau_at_star, as.numeric(f_star)),
    c(0, 1,           as.numeric(df_star))
  )
  colnames(A) <- colnames(X_orig)
  rownames(A) <- c("value_at_tau0", "slope_at_tau0")

  Nmat <- null_space(A)
  Z    <- X_orig %*% Nmat
  colnames(Z) <- paste0("z", seq_len(ncol(Z)))
  df_fit <- cbind(df, as.data.frame(Z))

  f_ns <- as.formula(paste0(
    "Surv(start, stop, event_vrs) ~ ",
    paste(c(colnames(Z), "strata(Group)", "cluster(RUN)"), collapse = " + ")
  ))

  fit_ns    <- coxph(f_ns, data = df_fit, ties = "efron", robust = TRUE,
                     x = TRUE, model = TRUE)
  fit_const <- coxph(
    Surv(start, stop, event_vrs) ~ inmunizado + strata(Group) + cluster(RUN),
    data = df, ties = "efron", robust = TRUE
  )

  alpha_full <- coef(fit_ns)
  keep_alpha <- !is.na(alpha_full)
  alpha      <- alpha_full[keep_alpha]
  V_alpha    <- vcov(fit_ns)[keep_alpha, keep_alpha, drop = FALSE]
  N_keep     <- Nmat[, keep_alpha, drop = FALSE]

  theta_hat <- as.vector(N_keep %*% alpha)
  names(theta_hat) <- colnames(X_orig)
  V_theta <- N_keep %*% V_alpha %*% t(N_keep)
  rownames(V_theta) <- colnames(V_theta) <- colnames(X_orig)

  # ── CORRECCIÓN: grilla acotada a tau_star_m ──
  tau_grid_m <- seq(0, tau_star_m, by = tau_grid_by)

  f_grid     <- do.call(f_basis_fn, c(list(tau = tau_grid_m, tau_star_m = tau_star_m),
                                       c(f_basis_args, list(knots = knots_out))))$B
  active_g   <- as.numeric(tau_grid_m < tau_star_m)

  Xg <- cbind(
    beta0 = active_g,
    beta1 = active_g * tau_grid_m,
    sweep(f_grid, 1, active_g, `*`)
  )
  colnames(Xg) <- colnames(X_orig)

  eta    <- as.vector(Xg %*% theta_hat)
  se_eta <- sqrt(rowSums((Xg %*% V_theta) * Xg))

  ve_curve <- data.frame(
    family = family_label,
    tau_m  = tau_grid_m,
    logHR  = eta, se = se_eta,
    HR     = exp(eta),
    HR_lo  = exp(eta - 1.96 * se_eta),
    HR_hi  = exp(eta + 1.96 * se_eta),
    VE     = 1 - exp(eta),
    VE_lo  = 1 - exp(eta + 1.96 * se_eta),
    VE_hi  = 1 - exp(eta - 1.96 * se_eta)
  )

  model_summary <- data.frame(
    family        = family_label,
    tau_0         = tau_star_m,
    n_params_free = sum(keep_alpha),
    loglik_ns     = fit_ns$loglik[2],
    loglik_const  = fit_const$loglik[2],
    LRT_stat      = 2 * (fit_ns$loglik[2] - fit_const$loglik[2]),
    AIC           = -2 * fit_ns$loglik[2] + 2 * sum(keep_alpha),
    p_value       = pchisq(
      2 * (fit_ns$loglik[2] - fit_const$loglik[2]),
      df         = sum(keep_alpha) - sum(!is.na(coef(fit_const))),
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
# 9) COMPARAR FAMILIAS
# ============================================================

data_w14 <- prepared_by_width[["w14"]]

res_families <- list(
  spline = fit_waning_generic(
    data_w14, split_width_days = 14, tau_star_m = 15,
    f_basis_fn = make_f_spline, f_basis_args = list(spline_df = 3),
    family_label = "Spline (df=3)"
  ),
  log = fit_waning_generic(
    data_w14, split_width_days = 14, tau_star_m = 15,
    f_basis_fn = make_f_log, f_basis_args = list(),
    family_label = "Logarítmica"
  ),
  weibull = fit_waning_generic(
    data_w14, split_width_days = 14, tau_star_m = 15,
    f_basis_fn = make_f_weibull, f_basis_args = list(gamma = 2),
    family_label = "Weibull (γ=2)"
  ),
  exponential = fit_waning_generic(
    data_w14, split_width_days = 14, tau_star_m = 15,
    f_basis_fn = make_f_exponential, f_basis_args = list(lambda = 0.3),
    family_label = "Exponencial (λ=0.3)"
  )
)

comp_table <- bind_rows(lapply(res_families, `[[`, "model_summary")) %>%
  mutate(delta_AIC = AIC - min(AIC)) %>%
  arrange(AIC)

print(comp_table)

df_curves <- bind_rows(lapply(res_families, `[[`, "ve_curve"))

familia_colors <- c(
  "Spline (df=3)"       = "#1B4F72",
  "Logarítmica"         = "#C0392B",
  "Weibull (γ=2)"       = "#1E8449",
  "Exponencial (λ=0.3)" = "#7D3C98"
)

# ── Overlay sin IC ──
p_overlay <- df_curves %>%
  ggplot(aes(tau_m, VE, color = family)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey60", linewidth = 0.4) +
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
    plot.title         = element_text(face = "bold", size = 20),
    plot.subtitle      = element_text(color = "grey35", size = 15),
    legend.position    = "bottom",
    legend.text        = element_text(size = 12),
    legend.key.width   = unit(1.1, "cm"),
    panel.grid.major.y = element_line(color = "grey93", linewidth = 0.3),
    plot.margin        = margin(12, 16, 8, 10),
    aspect.ratio       = 0.55
  )

print(p_overlay)

# ── Facetas con IC ──
library(ggplot2)
library(scales)
library(grid)

p_facet <- df_curves %>%
  ggplot(aes(x = tau_m, y = VE, color = family, fill = family)) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.5
  ) +
  geom_ribbon(
    aes(ymin = VE_lo, ymax = VE_hi),
    alpha = 0.18,
    linewidth = 0
  ) +
  geom_line(linewidth = 1.15) +
  facet_wrap(~ family, ncol = 2) +
  scale_color_manual(values = familia_colors) +
  scale_fill_manual(values = familia_colors) +
  scale_y_continuous(
    labels = label_percent(accuracy = 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1)
  ) +
  scale_x_continuous(
    breaks = seq(0, 15, by = 3)
  ) +
  coord_cartesian(xlim = c(0, 15), ylim = c(0, 1.1)) +
  labs(
    x = "Months since immunization",
    y = expression(VE(tau) == 1 - HR(tau)),
    title = "Waning of vaccine effectiveness",
    subtitle = "Functional-form comparison. Shaded areas show pointwise 95% confidence intervals."
  ) +
  theme_bw(base_size = 25, base_family = "serif") +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0),
    plot.subtitle = element_text(size = 13, color = "grey30", margin = margin(b = 10)),
    axis.title.x = element_text(size = 15, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 15, face = "bold", margin = margin(r = 10)),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none",
    strip.background = element_rect(fill = "grey94", color = "grey70"),
    strip.text = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.35),
    panel.spacing = unit(1.1, "lines"),
    plot.margin = margin(14, 18, 12, 12)
  )

print(p_facet)
ggsave("ve_overlay.pdf",  p_overlay, width = 7.5, height = 5.0, device = cairo_pdf)
ggsave("ve_overlay.png",  p_overlay, width = 7.5, height = 5.0, dpi = 320)
ggsave("ve_facet.pdf",    p_facet,   width = 7.5, height = 6.0, device = cairo_pdf)
ggsave("ve_facet.png",    p_facet,   width = 7.5, height = 6.0, dpi = 320)

# ============================================================
# 10) TABLA DE COEFICIENTES
# ============================================================
library(kableExtra)

extract_coef_block <- function(res, family_label) {
  th   <- res$theta_hat
  se_t <- sqrt(diag(res$V_theta))
  z_t  <- th / se_t
  p_t  <- 2 * pnorm(abs(z_t), lower.tail = FALSE)
  stars <- case_when(
    p_t < 0.001 ~ "***", p_t < 0.01 ~ "**", p_t < 0.05 ~ "*", TRUE ~ ""
  )
  cell <- paste0(sprintf("%.3f", th), stars, "\n(", sprintf("%.3f", se_t), ")")
  data.frame(param = names(th), value = cell, stringsAsFactors = FALSE) %>%
    setNames(c("param", family_label))
}

tabs <- lapply(names(res_families), function(nm) {
  extract_coef_block(res_families[[nm]], res_families[[nm]]$model_summary$family)
})

tab_coef <- Reduce(function(a, b) full_join(a, b, by = "param"), tabs)

param_order  <- c("beta0", "beta1", paste0("f", 1:10))
param_labels <- c(beta0 = "β₀  (intercept)", beta1 = "β₁  (linear slope)",
                  f1 = "f₁(τ)", f2 = "f₂(τ)", f3 = "f₃(τ)")
tab_coef <- tab_coef %>%
  mutate(param = factor(param, levels = param_order)) %>%
  arrange(param) %>%
  mutate(param = as.character(param))
tab_coef$param <- dplyr::recode(tab_coef$param, !!!param_labels)

metrics <- bind_rows(lapply(res_families, `[[`, "model_summary")) %>%
  mutate(
    AIC    = round(-2 * loglik_ns + 2 * n_params_free, 2),
    loglik = round(loglik_ns, 2),
    LRT_p  = ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))
  )

make_metric_row <- function(label, values) {
  setNames(as.data.frame(t(c(label, as.character(values))),
                         stringsAsFactors = FALSE),
           c("param", colnames(tab_coef)[-1]))
}

rows_metrics <- bind_rows(
  make_metric_row("N (free params)", metrics$n_params_free),
  make_metric_row("Log-likelihood",  metrics$loglik),
  make_metric_row("AIC",             metrics$AIC),
  make_metric_row("LRT p-value",     metrics$LRT_p)
)

tab_coef[is.na(tab_coef)] <- "—"
tab_final  <- bind_rows(tab_coef, rows_metrics)
n_coef     <- nrow(tab_coef)
n_met      <- nrow(rows_metrics)
col_names  <- c("", colnames(tab_final)[-1])

html_out <- tab_final %>%
  kable("html", col.names = col_names,
        align   = c("l", rep("c", ncol(tab_final) - 1)),
        caption = "Model coefficients by functional form  ·  β₀ + β₁τ + β₂f(τ)  ·  τ₀ = 15 mo",
        escape  = FALSE) %>%
  kable_styling(full_width = FALSE,
                bootstrap_options = c("striped","hover","condensed"),
                font_size = 13) %>%
  row_spec(n_coef, extra_css = "border-bottom: 2px solid #555;") %>%
  row_spec(seq(n_coef + 1, n_coef + n_met),
           background = "#f5f5f5", italic = TRUE) %>%
  row_spec(0, bold = TRUE, background = "#2C3E50", color = "white") %>%
  column_spec(1, bold = TRUE, width = "12em") %>%
  as.character()

writeLines(paste0(
  '<!DOCTYPE html><html><head><meta charset="utf-8">',
  '<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css">',
  '<style>body{padding:30px;font-family:Georgia,serif;}',
  'td{white-space:pre-line;vertical-align:middle;}',
  'caption{font-weight:bold;font-size:14px;margin-bottom:8px;}</style>',
  '</head><body>', html_out, '</body></html>'
), "coef_comparison_table.html")

cat("✓ Tabla guardada en coef_comparison_table.html\n")

cat("\nNota: estimate (SE). Significancia: * p<0.05  ** p<0.01  *** p<0.001\n")