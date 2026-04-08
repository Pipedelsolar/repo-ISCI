# ============================================================
# FPM WANING MODEL — Jennings et al. (2025) framework
# Script limpio y consolidado
# ============================================================
 
library(survival)
library(splines2)
library(dplyr)
library(ggplot2)
library(readr)
library(patchwork)
library(scales)
library(knitr)
 
# ============================================================
# 0) DATOS
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
 
df_raw$RUN        <- as.character(df_raw$RUN)
df_raw$Group      <- as.factor(df_raw$Group)
df_raw$start      <- as.numeric(df_raw$start)
df_raw$stop       <- as.numeric(df_raw$stop)
df_raw$event_vrs  <- as.integer(df_raw$event_vrs)
df_raw$inmunizado <- as.integer(df_raw$inmunizado)
 
if (!"t_inm" %in% names(df_raw)) {
  ref_date        <- as.Date("2024-04-01")
  df_raw$fechaInm <- as.Date(df_raw$fechaInm)
  df_raw$t_inm    <- as.numeric(df_raw$fechaInm - ref_date)
}
 
df_raw <- df_raw %>% filter(stop > start)
 
stopifnot(all(df_raw$event_vrs  %in% c(0, 1)))
stopifnot(all(df_raw$inmunizado %in% c(0, 1)))
 
# ============================================================
# 1) PREPARAR DATOS FPM
# Cada persona puede tener hasta 2 segmentos:
#   inmunizado=0 (antes de vacunarse)
#   inmunizado=1 (después de vacunarse)
# El tiempo del FPM es tau = tiempo desde vacunación (meses)
# ============================================================
 
df_fpm <- df_raw %>%
  arrange(RUN, start) %>%
  group_by(RUN, inmunizado) %>%
  summarise(
    t_entry        = min(start),
    t_exit         = max(stop),
    event          = as.integer(last(event_vrs) == 1),
    t_inm          = first(t_inm),
    Group          = first(Group),
    .groups        = "drop"
  ) %>%
  mutate(
    t_segment_days = t_exit - t_entry,
    tau_start_days = ifelse(inmunizado == 1, pmax(0, t_entry - t_inm), NA_real_),
    tau_stop_days  = ifelse(inmunizado == 1, pmax(0, t_exit  - t_inm), NA_real_),
    tau_start_m    = tau_start_days / 30.4375,
    tau_stop_m     = tau_stop_days  / 30.4375
  ) %>%
  filter(t_segment_days > 0)
 
cat("── Resumen df_fpm ──\n")
df_fpm %>%
  group_by(inmunizado) %>%
  summarise(
    n_segmentos  = n(),
    n_personas   = n_distinct(RUN),
    n_eventos    = sum(event),
    t_medio_dias = round(mean(t_segment_days), 1),
    .groups = "drop"
  ) %>%
  print()
 
# ============================================================
# 2) FUNCIÓN PRINCIPAL: fit_fpm_waning
# ============================================================
 
fit_fpm_waning <- function(df_fpm,
                            df_trt      = 4,
                            tau_star_m  = 15,
                            tau_grid_by = 0.05,
                            add_p95     = TRUE,
                            n_quad      = 200,   # puntos grilla integración
                            verbose     = TRUE) {
 
  # ── Datos brazo inmunizado ──
  df_imm <- df_fpm %>%
    filter(inmunizado == 1) %>%
    mutate(
      t  = pmax(tau_stop_m, 0.001),
      ev = as.integer(event)
    ) %>%
    filter(is.finite(t) & t > 0.001) %>%
    as.data.frame()
 
  if (verbose) {
    cat(sprintf("N inmunizados: %d   Eventos: %d\n",
                nrow(df_imm), sum(df_imm$ev)))
  }
 
  t_events_imm <- df_imm$t[df_imm$ev == 1]
  t_events_imm <- t_events_imm[is.finite(t_events_imm) & t_events_imm > 0]
 
  # ── Boundary knots ──
  bk <- c(0.001, tau_star_m)
 
  # ── Knots internos ──
  n_internal <- max(df_trt - 1L, 0L)
 
  if (n_internal > 0 && length(t_events_imm) >= 4) {
    t_ev_clip <- t_events_imm[t_events_imm > bk[1] & t_events_imm < bk[2]]
    probs     <- seq(0, 1, length.out = n_internal + 2L)[-c(1L, n_internal + 2L)]
    knots_    <- quantile(t_ev_clip, probs = probs, na.rm = TRUE)
 
    if (add_p95) {
      p95           <- quantile(t_ev_clip, 0.95, na.rm = TRUE)
      knots_[length(knots_)] <- p95
    }
    knots_ <- sort(unique(round(knots_, 6)))
    knots_ <- knots_[knots_ > bk[1] & knots_ < bk[2]]
  } else {
    knots_ <- NULL
  }
 
  if (verbose) {
    cat(sprintf("Boundary: [%.3f, %.3f]   Internal knots: [%s]\n",
                bk[1], bk[2],
                paste(round(knots_, 3), collapse = ", ")))
  }
 
  # ── Base NCS en datos ──
  B_imm <- as.matrix(naturalSpline(
    df_imm$t,
    knots          = knots_,
    Boundary.knots = bk,
    intercept      = TRUE
  ))
  n_basis     <- ncol(B_imm)
  col_names_B <- paste0("ncs", seq_len(n_basis))
  colnames(B_imm) <- col_names_B
 
  if (verbose) cat(sprintf("n_basis = %d\n", n_basis))
 
  # ── Base NCS en grilla común (pre-calculada para integración vectorizada) ──
  t_max    <- max(df_imm$t)
  s_common <- seq(bk[1], t_max, length.out = n_quad)
 
  B_common <- as.matrix(naturalSpline(
    s_common,
    knots          = knots_,
    Boundary.knots = bk,
    intercept      = TRUE
  ))
 
  # ── Punto de partida: Weibull simple + coxph ──
  wb_init <- tryCatch(
    survreg(Surv(t, ev) ~ 1, data = df_imm, dist = "weibull"),
    error = function(e) NULL
  )
 
  if (!is.null(wb_init)) {
    log_shape_init <- -log(wb_init$scale)
    log_scale_init <- log(abs(exp(coef(wb_init)[1])) + 0.1)
  } else {
    log_shape_init <- 0
    log_scale_init <- log(mean(df_imm$t))
  }
 
  df_cox   <- cbind(df_imm[, c("t", "ev")], as.data.frame(B_imm))
  f_cox    <- as.formula(paste0("Surv(t, ev) ~ ",
                                paste(col_names_B, collapse = " + ")))
  cox_init <- tryCatch(
    coxph(f_cox, data = df_cox, ties = "efron"),
    error = function(e) NULL
  )
 
  gamma_init <- if (!is.null(cox_init)) {
    cx <- coef(cox_init)
    cx[is.na(cx)] <- 0
    as.numeric(cx)
  } else {
    rep(0, n_basis)
  }
 
  init_par <- c(log_shape_init, log_scale_init, gamma_init)
 
  # ── Log-likelihood vectorizada (Weibull + NCS aditivo) ──
  # log h(t) = log h0(t) + NCS(t)' gamma
  # H(t)     = integral_0^t h(s) ds  →  trapecio en grilla común
  make_neg_ll <- function(fixed_idx = integer(0),
                           fixed_val = numeric(0)) {
    function(par) {
      log_shape <- par[1]
      log_scale <- par[2]
      gamma     <- par[3:(2 + n_basis)]
 
      if (length(fixed_idx) > 0) gamma[fixed_idx] <- fixed_val
 
      shape <- exp(log_shape)
      scale <- exp(log_scale)
 
      if (!is.finite(shape) || !is.finite(scale) ||
          shape <= 0 || scale <= 0) return(1e10)
 
      # ── log-hazard en t_i (vectorizado) ──
      log_h_t <- log(shape) - log(scale) +
        (shape - 1) * log(df_imm$t / scale) +
        as.numeric(B_imm %*% gamma)
 
      # ── Hazard acumulado via grilla común ──
      log_h_s  <- log(shape) - log(scale) +
        (shape - 1) * log(s_common / scale) +
        as.numeric(B_common %*% gamma)
      h_s      <- exp(log_h_s)
 
      # Integral acumulada (trapecio)
      ds       <- diff(s_common)
      dH_steps <- ds * (head(h_s, -1) + tail(h_s, -1)) / 2
      H_common <- c(0, cumsum(dH_steps))
 
      # Interpolar H(t_i)
      H_t <- approx(s_common, H_common, xout = df_imm$t, rule = 2)$y
 
      ll <- sum(df_imm$ev * log_h_t) - sum(H_t)
 
      if (!is.finite(ll)) return(1e10)
      -ll
    }
  }
 
  neg_ll_unc <- make_neg_ll()
  neg_ll_con <- make_neg_ll(fixed_idx = c(1L, as.integer(n_basis)),
                             fixed_val = c(0, 0))
 
  if (verbose) cat("Optimizando unconstrained...\n")
  opt_unc <- optim(init_par, neg_ll_unc, method = "BFGS",
                   control = list(maxit = 3000, reltol = 1e-10),
                   hessian = TRUE)
 
  if (verbose) cat("Optimizando constrained...\n")
  opt_con <- optim(init_par, neg_ll_con, method = "BFGS",
                   control = list(maxit = 3000, reltol = 1e-10),
                   hessian = TRUE)
 
  loglik_unc <- -opt_unc$value
  loglik_con <- -opt_con$value
  lrt_stat   <- 2 * max(loglik_unc - loglik_con, 0)
 
  if (verbose) {
    cat(sprintf("  LL unc=%.3f   LL con=%.3f\n", loglik_unc, loglik_con))
    cat(sprintf("  LRT=%.3f   p=%.4g\n",
                lrt_stat, pchisq(lrt_stat, df = 2, lower.tail = FALSE)))
  }
 
  # ── Predicción en grilla ──
  tau_grid_m <- seq(0.01, ceiling(t_max), by = tau_grid_by)
 
  gamma_con           <- opt_con$par[3:(2 + n_basis)]
  gamma_con[c(1L, as.integer(n_basis))] <- 0   # forzar constraint exacto
 
  B_grid <- as.matrix(naturalSpline(
    tau_grid_m,
    knots          = knots_,
    Boundary.knots = bk,
    intercept      = TRUE
  ))
 
  log_HR                           <- as.numeric(B_grid %*% gamma_con)
  log_HR[tau_grid_m >= tau_star_m] <- 0   # HR = 1 después de tau*
 
  # SE via delta method (solo parámetros gamma libres)
  free_idx      <- setdiff(seq_len(n_basis), c(1L, as.integer(n_basis)))
  gamma_par_idx <- 2L + free_idx
 
  vcov_con <- tryCatch(solve(opt_con$hessian), error = function(e) NULL)
 
  if (!is.null(vcov_con) && max(gamma_par_idx) <= nrow(vcov_con)) {
    V_free   <- vcov_con[gamma_par_idx, gamma_par_idx, drop = FALSE]
    G        <- B_grid[, free_idx, drop = FALSE]
    se_logHR <- sqrt(pmax(0, rowSums((G %*% V_free) * G)))
  } else {
    se_logHR <- rep(NA_real_, length(tau_grid_m))
  }
  se_logHR[tau_grid_m >= tau_star_m] <- 0
 
  HR <- exp(log_HR)
 
  ve_curve <- data.frame(
    df_trt   = df_trt,
    tau_star = tau_star_m,
    p95_knot = add_p95,
    tau_m    = tau_grid_m,
    logHR    = log_HR,
    se       = se_logHR,
    HR       = HR,
    HR_lo    = exp(log_HR - 1.96 * se_logHR),
    HR_hi    = exp(log_HR + 1.96 * se_logHR),
    VE       = 1 - HR,
    VE_lo    = 1 - exp(log_HR + 1.96 * se_logHR),
    VE_hi    = 1 - exp(log_HR - 1.96 * se_logHR)
  )
 
  model_summary <- data.frame(
    df_trt     = df_trt,
    tau_star   = tau_star_m,
    p95_knot   = add_p95,
    n_imm      = nrow(df_imm),
    n_events   = sum(df_imm$ev),
    loglik_unc = round(loglik_unc, 3),
    loglik_con = round(loglik_con, 3),
    LRT_stat   = round(lrt_stat,   3),
    LRT_p      = pchisq(lrt_stat, df = 2, lower.tail = FALSE),
    AIC_unc    = round(-2 * loglik_unc + 2 * (2 + n_basis), 2),
    AIC_con    = round(-2 * loglik_con + 2 * (2 + length(free_idx)), 2)
  )
 
  list(
    opt_unc       = opt_unc,
    opt_con       = opt_con,
    ve_curve      = ve_curve,
    model_summary = model_summary,
    knots         = knots_,
    bk            = bk,
    n_basis       = n_basis,
    free_idx      = free_idx,
    df_imm        = df_imm,
    tau_star_m    = tau_star_m,
    B_imm         = B_imm,
    s_common      = s_common,
    B_common      = B_common
  )
}
 
# ============================================================
# 3) GRILLA DE MODELOS
# ============================================================
 
df_grid <- expand.grid(
  df_trt     = c(3, 4, 5),
  tau_star_m = c(13, 14, 15),
  add_p95    = c(FALSE, TRUE),
  stringsAsFactors = FALSE
)
 
results_fpm <- list()
 
for (i in seq_len(nrow(df_grid))) {
  key <- paste0(
    "df",   df_grid$df_trt[i],
    "_tau", df_grid$tau_star_m[i],
    "_p95", as.integer(df_grid$add_p95[i])
  )
 
  message(sprintf("[%d/%d] df=%d  tau*=%g  p95=%s",
                  i, nrow(df_grid),
                  df_grid$df_trt[i],
                  df_grid$tau_star_m[i],
                  df_grid$add_p95[i]))
 
  results_fpm[[key]] <- tryCatch(
    fit_fpm_waning(
      df_fpm      = df_fpm,
      df_trt      = df_grid$df_trt[i],
      tau_star_m  = df_grid$tau_star_m[i],
      add_p95     = df_grid$add_p95[i],
      n_quad      = 200,
      verbose     = FALSE
    ),
    error = function(e) {
      message("  ERROR: ", conditionMessage(e))
      NULL
    }
  )
}
 
# ============================================================
# 4) CONSOLIDAR RESULTADOS
# ============================================================
 
all_summaries <- bind_rows(
  lapply(names(results_fpm), function(k) {
    x <- results_fpm[[k]]
    if (is.null(x)) return(NULL)
    x$model_summary %>% mutate(model_key = k)
  })
) %>%
  mutate(
    delta_AIC_con = AIC_con - min(AIC_con, na.rm = TRUE),
    p95_label     = ifelse(p95_knot, "With p95 knot", "Without p95 knot")
  ) %>%
  arrange(AIC_con)
 
all_curves <- bind_rows(
  lapply(names(results_fpm), function(k) {
    x <- results_fpm[[k]]
    if (is.null(x)) return(NULL)
    x$ve_curve %>% mutate(model_key = k)
  })
)
 
cat("\n── Resumen de modelos ──\n")
print(all_summaries)
 
# ============================================================
# 5) FIGURA 1 — Sensitivity grid
# ============================================================
 
p_sensitivity <- all_curves %>%
  mutate(
    df_label  = factor(paste0("df = ", df_trt),
                       levels = paste0("df = ", c(3, 4, 5))),
    tau_label = factor(paste0("τ₀ = ", tau_star, " mo"),
                       levels = paste0("τ₀ = ", c(13, 14, 15), " mo")),
    p95_label = ifelse(p95_knot, "With p95 knot", "Without p95 knot")
  ) %>%
  filter(tau_m <= 15) %>%
  ggplot(aes(tau_m, VE,
             color    = p95_label,
             fill     = p95_label,
             linetype = p95_label)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.35) +
  geom_ribbon(aes(ymin = VE_lo, ymax = VE_hi),
              alpha = 0.10, linewidth = 0) +
  geom_line(linewidth = 0.85) +
  facet_grid(rows = vars(df_label), cols = vars(tau_label)) +
  scale_color_manual(values = c("#1B4F72", "#C0392B")) +
  scale_fill_manual(values  = c("#1B4F72", "#C0392B")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     breaks = c(0, 0.5, 1)) +
  scale_x_continuous(breaks = seq(0, 15, by = 5)) +
  coord_cartesian(xlim = c(0, 15), ylim = c(-0.2, 1.1)) +
  labs(
    x        = "Months since immunization",
    y        = "VE(τ) = 1 − HR(τ)",
    color    = NULL, fill = NULL, linetype = NULL,
    title    = "FPM waning model — sensitivity to df and τ₀",
    subtitle = "NCS with/without p95 knot  ·  Jennings et al. (2025) framework"
  ) +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(color = "grey35", size = 9),
    legend.position  = "bottom",
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    strip.text       = element_text(face = "bold", size = 9),
    panel.grid.minor = element_blank()
  )
 
print(p_sensitivity)
ggsave("fpm_sensitivity.pdf", p_sensitivity,
       width = 8, height = 7, device = cairo_pdf)
 
# ============================================================
# 6) FIGURA 2 — Mejor modelo (paper-ready)
# ============================================================
 
best_key   <- all_summaries$model_key[1]
res_best   <- results_fpm[[best_key]]
curve_best <- res_best$ve_curve
 
support_hist <- res_best$df_imm %>%
  mutate(month_bin = floor(t)) %>%
  group_by(month_bin) %>%
  summarise(person_time = n(), events = sum(ev), .groups = "drop") %>%
  mutate(rel_pt = person_time / max(person_time) * 0.12)
 
p_best <- ggplot() +
  geom_rect(
    data = support_hist,
    aes(xmin = month_bin, xmax = month_bin + 0.85,
        ymin = -0.05, ymax = -0.05 + rel_pt),
    fill = "grey82", color = NA, alpha = 0.9
  ) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  geom_ribbon(data = curve_best,
              aes(x = tau_m, ymin = VE_lo, ymax = VE_hi),
              fill = "#2E86C1", alpha = 0.15) +
  geom_line(data = curve_best,
            aes(x = tau_m, y = VE),
            color = "#1B4F72", linewidth = 1.1) +
  geom_vline(xintercept = res_best$tau_star_m,
             linetype = "dotted", color = "grey45", linewidth = 0.5) +
  annotate("text",
           x = res_best$tau_star_m + 0.2, y = 0.90,
           label = sprintf("τ₀ = %g mo", res_best$tau_star_m),
           hjust = 0, size = 3.3, color = "grey40", family = "serif") +
  annotate("text", x = 0.2, y = -0.02,
           label = "Person-time (relative)",
           hjust = 0, vjust = 1, size = 2.6,
           color = "grey55", family = "serif") +
  scale_x_continuous(breaks = seq(0, 15, by = 3), expand = c(0.01, 0)) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     breaks = c(-0.25, 0, 0.25, 0.5, 0.75, 1)) +
  coord_cartesian(xlim = c(0, res_best$tau_star_m + 0.5),
                  ylim = c(-0.12, 1.05)) +
  labs(
    x        = "Months since immunization",
    y        = "Vaccine effectiveness  VE(τ) = 1 − HR(τ)",
    title    = "Waning of vaccine effectiveness — FPM",
    subtitle = sprintf("NCS df=%d · τ₀=%g mo · p95=%s · 95%% CI via delta method",
                       res_best$model_summary$df_trt,
                       res_best$tau_star_m,
                       ifelse(res_best$model_summary$p95_knot, "yes", "no")),
    caption  = "Bars: relative person-time (immunized arm). Shaded: pointwise 95% CI."
  ) +
  theme_classic(base_size = 12, base_family = "serif") +
  theme(
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(color = "grey35", size = 9.5),
    plot.caption       = element_text(color = "grey45", size = 8, hjust = 0),
    panel.grid.major.y = element_line(color = "grey93", linewidth = 0.3),
    aspect.ratio       = 0.55,
    plot.margin        = margin(12, 16, 8, 10)
  )
 
print(p_best)
ggsave("fpm_ve_best.pdf", p_best,
       width = 7, height = 4.5, device = cairo_pdf)
 
# ============================================================
# 7) FIGURA 3 — HR unconstrained vs constrained
# ============================================================
 
tau_g  <- seq(0.01, res_best$tau_star_m + 1, by = 0.05)
 
B_g <- as.matrix(naturalSpline(
  tau_g,
  knots          = res_best$knots,
  Boundary.knots = res_best$bk,
  intercept      = TRUE
))
 
par_unc      <- res_best$opt_unc$par
par_con      <- res_best$opt_con$par
 
gamma_unc    <- par_unc[3:(2 + res_best$n_basis)]
gamma_con    <- par_con[3:(2 + res_best$n_basis)]
gamma_con[c(1L, as.integer(res_best$n_basis))] <- 0
 
HR_unc <- exp(as.numeric(B_g %*% gamma_unc))
HR_con <- exp(as.numeric(B_g %*% gamma_con))
HR_con[tau_g >= res_best$tau_star_m] <- 1
 
df_HR <- data.frame(
  tau    = rep(tau_g, 2),
  HR     = c(HR_unc, HR_con),
  model  = rep(c("Unconstrained", "Constrained"), each = length(tau_g))
)
 
p_HR <- ggplot(df_HR, aes(tau, HR, color = model, linetype = model)) +
  geom_hline(yintercept = 1, linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = res_best$tau_star_m, linetype = "dotted",
             color = "grey45", linewidth = 0.5) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = c("#1B4F72", "#C0392B")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  coord_cartesian(ylim = c(0, 1.5),
                  xlim = c(0, res_best$tau_star_m + 1)) +
  labs(x = "Months since immunization", y = "HR(τ)",
       title = "Unconstrained vs constrained HR — FPM",
       color = NULL, linetype = NULL) +
  theme_classic(base_size = 11, base_family = "serif") +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom",
        panel.grid.major.y = element_line(color = "grey93", linewidth = 0.3))
 
print(p_HR)
ggsave("fpm_HR_comparison.pdf", p_HR,
       width = 6, height = 4, device = cairo_pdf)
 
# ============================================================
# 8) TABLA COMPARATIVA
# ============================================================
 
all_summaries %>%
  select(df_trt, tau_star, p95_label, n_events,
         loglik_unc, loglik_con, LRT_stat, LRT_p,
         AIC_con, delta_AIC_con) %>%
  mutate(
    LRT_p         = ifelse(LRT_p < 0.001, "<0.001",
                           sprintf("%.3f", LRT_p)),
    delta_AIC_con = round(delta_AIC_con, 2)
  ) %>%
  kable("simple",
        col.names = c("df", "τ₀", "p95 knot", "Events",
                      "LL uncons.", "LL cons.",
                      "LRT χ²", "p-value",
                      "AIC (cons.)", "ΔAIC"),
        digits = 3)