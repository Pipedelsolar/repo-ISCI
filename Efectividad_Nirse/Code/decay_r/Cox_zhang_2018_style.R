library(survival)
library(splines)
library(ggplot2)
library(dplyr)
library(readr)
library(scales)
library(patchwork)
# install.packages("future")
# # Paso 1: forzar instalación binaria de parallelly actualizado
# install.packages("parallelly", type = "binary")

# # Paso 2: luego future y sus dependencias en binario
# install.packages("future",       type = "binary")
# install.packages("future.apply", type = "binary")

# # Paso 3: lava en binario
# install.packages("lava", type = "binary")

# # Paso 4: timereg
# install.packages("timereg", type = "binary")

# # Verificar
library(timereg)
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

# Tiempo de inmunización
if (!"t_inm" %in% names(df_raw)) {
  ref_date        <- as.Date("2024-04-01")
  df_raw$fechaInm <- as.Date(df_raw$fechaInm)
  df_raw$t_inm    <- as.numeric(df_raw$fechaInm - ref_date)
}

# tau = meses desde vacunación (0 para no inmunizados)
df_raw <- df_raw %>%
  mutate(
    tau_start_days = ifelse(inmunizado == 1, pmax(0, start - t_inm), 0),
    tau_stop_days  = ifelse(inmunizado == 1, pmax(0, stop  - t_inm), 0),
    tau_eval_days  = ifelse(inmunizado == 1,
                            0.5 * (tau_start_days + tau_stop_days), 0),
    tau_eval_m     = pmax(0, tau_eval_days / 30.4375)
  )

cat("── Resumen general ──\n")
cat("Personas únicas:", n_distinct(df_raw$RUN), "\n")
cat("Eventos totales:", sum(df_raw$event_vrs), "\n")
cat("Eventos inmunizados:", sum(df_raw$event_vrs[df_raw$inmunizado == 1]), "\n")
cat("Eventos controles:", sum(df_raw$event_vrs[df_raw$inmunizado == 0]), "\n")

# ============================================================
# 1) MODELO BASE — Cox con inmunizado como covariate constante
#    (verifica que hay efecto antes de modelar el waning)
# ============================================================

fit_base <- coxph(
  Surv(start, stop, event_vrs) ~ inmunizado + strata(Group) + cluster(RUN),
  data   = df_raw,
  ties   = "efron",
  robust = TRUE
)

cat("\n── Modelo base (efecto constante) ──\n")
print(summary(fit_base))

# ============================================================
# 2) TEST DE PROPORCIONALIDAD (Zhang et al., sección Schoenfeld)
#    ¿El efecto de inmunizado varía con el tiempo?
# ============================================================

zph_base <- cox.zph(fit_base, transform = "km")

cat("\n── Test de proporcionalidad (Schoenfeld) ──\n")
print(zph_base)

# ============================================================
# 3) FIGURA 1: Schoenfeld residuals — replica Figure 1 del paper
# ============================================================


# Extraer datos del plot para ggplot
zph_data <- data.frame(
  time  = zph_base$x,
  beta  = zph_base$y[, "inmunizado"],
#   se    = sqrt(zph_base$var[, "inmunizado"]) # si disponible
#   se    = sqrt(zph_base$var[, 1]) # si disponible
)

hr_const <- coef(fit_base)["inmunizado"]

p_schoenfeld <- ggplot(zph_data, aes(x = time, y = beta)) +
  geom_point(alpha = 0.25, size = 1.8, color = "grey40") +
  geom_smooth(method = "loess", span = 0.5,
              color = "#1B4F72", linewidth = 1.5,
              fill  = "#2E86C1", alpha = 0.15) +
  geom_hline(yintercept = 0,        linetype = "dotted",
             color = "black",  linewidth = 0.8) +
  geom_hline(yintercept = hr_const, linetype = "dashed",
             color = "#1E8449", linewidth = 0.8) +
  coord_cartesian(
    ylim = c(-3, 0.5)    # fuerza el rango Y directamente
  ) +
  labs(
    x        = "Time (KM-transformed scale)",
    y        = "β(t) for inmunizado",
    title    = "Time-varying coefficient — Schoenfeld residuals",
    subtitle = sprintf("Global PH test p = %.4f  ·  Dashed = average HR",
                       zph_base$table["inmunizado", "p"]),
  ) +
  annotate("text", x = min(zph_data_plot$time), y = 0.15,
           label = "Null effect", hjust = 0, size = 9,
           color = "black", family = "serif") +
  annotate("text", x = min(zph_data_plot$time), y = hr_const + 0.15,
           label = "Average HR", hjust = 0, size = 9,
           color = "#1E8449", family = "serif") +
  theme_classic(base_size = 20, base_family = "serif") +
  theme(
    plot.title    = element_text(face = "bold", size = 32),
    plot.subtitle = element_text(color = "grey35", size = 24),
    plot.caption  = element_text(color = "grey45", size = 20, hjust = 0),
    axis.title    = element_text(size = 28),
    axis.text     = element_text(size = 26),
    aspect.ratio  = 0.6    # más ancho que alto — reduce el alargamiento
  )

print(p_schoenfeld)
ggsave("schoenfeld_residuals.pdf", p_schoenfeld,
       width = 12, height = 6, device = cairo_pdf)   # dimensiones más anchas
ggsave("schoenfeld_residuals.png", p_schoenfeld,
       width = 12, height = 6, dpi = 200)


# ============================================================
# 4) MÉTODO A — Step function (Zhang et al., sección "Step function")
#    Replica exactamente lo del paper con survSplit()
#    Cortes cada 3 meses post-vacunación
# ============================================================

# Preparar dataset solo con brazo inmunizado para estimar VE(tau)
# El truco: usamos tau como eje de tiempo para el inmunizado
# y añadimos la interacción inmunizado:strata(tau_bin)

cut_points_m <- c(3, 6, 9, 12, 15)   # meses — ajustar según datos
cut_points_d <- cut_points_m * 30.4375

# Crear variable de bin de tau en df_raw
df_raw <- df_raw %>%
  mutate(
    tau_bin = cut(
      tau_eval_m,
      breaks = c(-Inf, cut_points_m, Inf),
      labels = c(
        paste0("0-", cut_points_m[1], "m"),
        paste0(cut_points_m[-length(cut_points_m)], "-",
               cut_points_m[-1], "m"),
        paste0(">", tail(cut_points_m, 1), "m")
      ),
      right  = FALSE
    ),
    # Para no inmunizados: bin de referencia
    tau_bin = as.character(tau_bin),
    tau_bin = ifelse(inmunizado == 0, "Control", tau_bin),
    tau_bin = factor(tau_bin)
  )

cat("\n── Distribución de intervalos tau ──\n")
df_raw %>%
  filter(inmunizado == 1) %>%
  group_by(tau_bin) %>%
  summarise(
    n_filas  = n(),
    eventos  = sum(event_vrs),
    pt_dias  = round(sum(stop - start)),
    .groups  = "drop"
  ) %>%
  print()

# Modelo step function: interacción inmunizado × tau_bin
# log-HR(tau) = sum_k beta_k * I(tau in bin_k) * inmunizado
fit_step <- coxph(
  Surv(start, stop, event_vrs) ~
    inmunizado:tau_bin + strata(Group) + cluster(RUN),
  data   = df_raw,
  ties   = "efron",
  robust = TRUE
)

cat("\n── Modelo step function ──\n")
print(summary(fit_step))

# Extraer HR y VE por bin
coef_step <- coef(fit_step)
se_step   <- sqrt(diag(vcov(fit_step)))

idx_bins  <- grep("inmunizado:tau_bin", names(coef_step))
bin_names <- gsub("inmunizado:tau_bin", "", names(coef_step)[idx_bins])

df_step_ve <- data.frame(
  tau_bin  = bin_names,
  logHR    = coef_step[idx_bins],
  se       = se_step[idx_bins]
) %>%
  # Filtrar Control y >15m que no tienen tau_mid numérico
  filter(!tau_bin %in% c("Control", ">15m")) %>%
  mutate(
    # Extraer límites directamente del nombre del bin
    lo      = as.numeric(gsub("^(\\d+\\.?\\d*)-.*", "\\1", tau_bin)),
    hi      = as.numeric(gsub("^\\d+\\.?\\d*-(\\d+\\.?\\d*)m$", "\\1", tau_bin)),
    tau_mid = (lo + hi) / 2,
    HR      = exp(logHR),
    HR_lo   = exp(logHR - 1.96 * se),
    HR_hi   = exp(logHR + 1.96 * se),
    VE      = 1 - HR,
    VE_lo   = 1 - HR_hi,
    VE_hi   = 1 - HR_lo
  ) %>%
  arrange(tau_mid) %>%
  mutate(tau_bin = factor(tau_bin, levels = tau_bin))

cat("\n── VE por bin ordenado ──\n")
print(df_step_ve %>% select(tau_bin, tau_mid, VE, VE_lo, VE_hi))

# ============================================================
# 5) FIGURA 2: VE step function — replica Figure 2 del paper
# ============================================================

df_step_ve <- df_step_ve %>%
  arrange(tau_mid)

p_step <- ggplot(df_step_ve, aes(x = tau_mid, y = VE)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.6) +
  geom_errorbar(aes(ymin = VE_lo, ymax = VE_hi),
                width = 0.3, color = "#1B4F72", linewidth = 1.0) +
  geom_point(size = 5, color = "#1B4F72", fill = "white",
             shape = 21, stroke = 2) +
  geom_step(aes(x = tau_mid, y = VE),
            color = "#1B4F72", linewidth = 1.0,
            direction = "mid") +
  scale_x_continuous(
    breaks = df_step_ve$tau_mid,          # solo los puntos del df
    labels = df_step_ve$tau_bin           # etiquetas exactas del df
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    breaks = seq(-0.25, 1, by = 0.25)
  ) +
  coord_cartesian(
    xlim = c(min(df_step_ve$tau_mid) - 0.5,
             max(df_step_ve$tau_mid) + 0.5),
    ylim = c(-0.3, 1.1)
  ) +
  labs(
    x        = "Months since immunization",
    y        = "Vaccine effectiveness  VE = 1 − HR",
    title    = "Vaccine effectiveness — step function",
    subtitle = "HR estimated per time bin  ·  95% CI",
  ) +
  theme_classic(base_size = 28, base_family = "serif") +
  theme(
    plot.title    = element_text(face = "bold", size = 30),
    plot.subtitle = element_text(color = "grey35", size = 22),
    plot.caption  = element_text(color = "grey45", size = 18, hjust = 0),
    axis.title    = element_text(size = 26),
    axis.text.x   = element_text(angle = 30, hjust = 1, size = 22),
    axis.text.y   = element_text(size = 22),
    aspect.ratio  = 0.45    # ancho, no alargado
  )

print(p_step)
ggsave("ve_step_function.pdf", p_step,
       width = 14, height = 7, device = cairo_pdf)
ggsave("ve_step_function.png", p_step,
       width = 14, height = 7, dpi = 200)

# ============================================================
# 6) MÉTODO B — Función continua con tt()
#    (Zhang et al., sección "Continuous function")
#    log-HR(tau) = beta * g(tau) donde g es una función de tau
# ============================================================

# Para usar tt() necesitamos un dataset donde el tiempo sea tau
# Usamos solo el brazo inmunizado con tau como eje de tiempo
df_tt <- df_raw %>%
  mutate(
    # tt() trabaja con tiempo de análisis = stop
    # Usamos tau_eval_m como covariable time-varying
    # La interacción inmunizado × g(tau) se modela via tt()
  )

# Modelo B1: g(tau) = log(tau + 1) — efecto logarítmico
fit_tt_log <- coxph(
  Surv(start, stop, event_vrs) ~
    inmunizado + tt(inmunizado) + strata(Group) + cluster(RUN),
  data   = df_raw,
  ties   = "efron",
  robust = TRUE,
  tt     = function(x, t, ...) {
    # tau en meses desde inmunización
    tau_m <- pmax(0, (t - df_raw$t_inm[match(
      df_raw$RUN, df_raw$RUN)]) / 30.4375)
    x * log(tau_m + 1)
  }
)

# Versión más simple y robusta: crear la covariable explícitamente
# beta(tau) = beta0 + beta1 * log(tau + 1)
df_raw <- df_raw %>%
  mutate(
    inm_x_logtau = inmunizado * log(tau_eval_m + 1),
    inm_x_tau    = inmunizado * tau_eval_m,
    inm_x_tau2   = inmunizado * tau_eval_m^2,
    inm_x_sqrttau = inmunizado * sqrt(tau_eval_m)
  )

# Modelo B1: log-HR(tau) = beta0 + beta1 * log(tau+1)
fit_log <- coxph(
  Surv(start, stop, event_vrs) ~
    inmunizado + inm_x_logtau + strata(Group) + cluster(RUN),
  data   = df_raw,
  ties   = "efron",
  robust = TRUE
)

# Modelo B2: log-HR(tau) = beta0 + beta1 * tau  (lineal)
fit_linear <- coxph(
  Surv(start, stop, event_vrs) ~
    inmunizado + inm_x_tau + strata(Group) + cluster(RUN),
  data   = df_raw,
  ties   = "efron",
  robust = TRUE
)

# Modelo B3: log-HR(tau) = beta0 + beta1 * tau + beta2 * tau^2
fit_quad <- coxph(
  Surv(start, stop, event_vrs) ~
    inmunizado + inm_x_tau + inm_x_tau2 + strata(Group) + cluster(RUN),
  data   = df_raw,
  ties   = "efron",
  robust = TRUE
)

# Modelo B4: log-HR(tau) = beta0 + beta1 * sqrt(tau)
fit_sqrt <- coxph(
  Surv(start, stop, event_vrs) ~
    inmunizado + inm_x_sqrttau + strata(Group) + cluster(RUN),
  data   = df_raw,
  ties   = "efron",
  robust = TRUE
)

# Modelo B5: log-HR(tau) = beta0 + beta1*tau + spline(tau)
# Spline cúbico residualizado con restricciones en tau_star

# Parámetros
spline_df  <- 3
tau_star_m <- 15

# Preparar base spline (reutilizando funciones del script anterior)
tau_eval_clip <- pmin(df_raw$tau_eval_m, tau_star_m)
active_idx    <- df_raw$inmunizado == 1 & df_raw$tau_eval_m < tau_star_m
tau_imm       <- tau_eval_clip[active_idx]
tau_imm       <- tau_imm[is.finite(tau_imm)]

bs_support_obj <- make_bs_basis(tau_imm, spline_df, tau_star_m, knots = NULL)
knots_spl      <- bs_support_obj$knots
B_support      <- bs_support_obj$B

B_all <- make_bs_basis(tau_eval_clip, spline_df, tau_star_m, knots_spl)$B

rez       <- residualize_basis_against_linear(B_support, tau_imm, B_all, tau_eval_clip)
B_nl_all  <- rez$B_nl_all
proj_coef <- rez$proj_coef

# Añadir columnas spline al df_raw
spline_cols <- paste0("inm_x_spl", seq_len(ncol(B_nl_all)))
for (j in seq_len(ncol(B_nl_all))) {
  df_raw[[spline_cols[j]]] <- as.numeric(df_raw$inmunizado == 1) * B_nl_all[, j]
}

# Ajustar modelo
fit_spline <- coxph(
  as.formula(paste0(
    "Surv(start, stop, event_vrs) ~ inmunizado + inm_x_tau + ",
    paste(spline_cols, collapse = " + "),
    " + strata(Group) + cluster(RUN)"
  )),
  data   = df_raw,
  ties   = "efron",
  robust = TRUE
)

cat("\nB5 — Spline cúbico (df=3):\n")
print(summary(fit_spline)$coefficients)

cat("\n── Modelos función continua ──\n")
cat("\nB1 — Logarítmico:\n");   print(summary(fit_log)$coefficients)
cat("\nB2 — Lineal:\n");         print(summary(fit_linear)$coefficients)
cat("\nB3 — Cuadrático:\n");     print(summary(fit_quad)$coefficients)
cat("\nB4 — Raíz cuadrada:\n");  print(summary(fit_sqrt)$coefficients)

# Comparación AIC
aic_table <- data.frame(
  modelo = c("Constante", "Log(tau+1)", "Lineal", "Cuadrático", "Sqrt(tau)"),
  loglik = c(
    fit_base$loglik[2],
    fit_log$loglik[2],
    fit_linear$loglik[2],
    fit_quad$loglik[2],
    fit_sqrt$loglik[2]
  ),
  df = c(1, 2, 2, 3, 2)
) %>%
  mutate(
    AIC      = -2 * loglik + 2 * df,
    delta_AIC = AIC - min(AIC)
  ) %>%
  arrange(AIC)

cat("\n── Comparación AIC ──\n")
print(aic_table)

# ============================================================
# 7) PREDICCIÓN DE VE(tau) PARA CADA MODELO CONTINUO
# ============================================================

tau_grid <- seq(0, 15, by = 0.1)

# Función auxiliar: predice log-HR en grilla dado modelo y coefs
predict_ve <- function(fit, tau_grid, model_name) {

  co  <- coef(fit)
  vc  <- vcov(fit)

  # Identificar coeficientes relevantes
  nm  <- names(co)

  # Construir X en grilla
  if (model_name == "Log(tau+1)") {
    X <- cbind(1, log(tau_grid + 1))
    idx <- match(c("inmunizado", "inm_x_logtau"), nm)
  } else if (model_name == "Lineal") {
    X <- cbind(1, tau_grid)
    idx <- match(c("inmunizado", "inm_x_tau"), nm)
  } else if (model_name == "Cuadrático") {
    X <- cbind(1, tau_grid, tau_grid^2)
    idx <- match(c("inmunizado", "inm_x_tau", "inm_x_tau2"), nm)
  } else if (model_name == "Sqrt(tau)") {
    X <- cbind(1, sqrt(tau_grid))
    idx <- match(c("inmunizado", "inm_x_sqrttau"), nm)
  } else if (model_name == "Spline (df=3)") {
    X   <- cbind(1, tau_grid, eval_nonlinear_basis(
        tau_grid, knots_spl, spline_df, tau_star_m, proj_coef))
    idx <- match(c("inmunizado", "inm_x_tau", spline_cols), names(coef(fit_spline)))
    idx <- idx[!is.na(idx)]
    b   <- coef(fit_spline)[idx]
    V   <- vcov(fit_spline)[idx, idx, drop = FALSE]
    X   <- X[, seq_len(length(idx)), drop = FALSE]
  } else {  # Constante
    X <- matrix(1, nrow = length(tau_grid), ncol = 1)
    idx <- match("inmunizado", nm)
  } 

  idx <- idx[!is.na(idx)]
  b   <- co[idx]
  V   <- vc[idx, idx, drop = FALSE]
  X   <- X[, seq_len(length(idx)), drop = FALSE]

  eta    <- as.numeric(X %*% b)
  se_eta <- sqrt(rowSums((X %*% V) * X))

  data.frame(
    modelo = model_name,
    tau_m  = tau_grid,
    logHR  = eta,
    se     = se_eta,
    HR     = exp(eta),
    HR_lo  = exp(eta - 1.96 * se_eta),
    HR_hi  = exp(eta + 1.96 * se_eta),
    VE     = 1 - exp(eta),
    VE_lo  = 1 - exp(eta + 1.96 * se_eta),
    VE_hi  = 1 - exp(eta - 1.96 * se_eta)
  )
}

df_ve_curves <- bind_rows(
  predict_ve(fit_base,   tau_grid, "Constante"),
  predict_ve(fit_log,    tau_grid, "Log(tau+1)"),
  predict_ve(fit_linear, tau_grid, "Lineal"),
  predict_ve(fit_quad,   tau_grid, "Cuadrático"),
  predict_ve(fit_sqrt,   tau_grid, "Sqrt(tau)")
)

df_ve_curves <- bind_rows(
  df_ve_curves,
  predict_ve(fit_spline, tau_grid, "Spline (df=3)")
)

# ============================================================
# 8) FIGURA 3: VE continua — todas las formas funcionales
# ============================================================

model_colors <- c(
  "Constante"   = "grey50",
  "Log(tau+1)"  = "#1B4F72",
  "Lineal"      = "#C0392B",
  "Cuadrático"  = "#1E8449",
  "Sqrt(tau)"   = "#7D3C98",
  "Spline (df=3)" = "#E67E22"
)

# Panel A: sin IC, todas juntas
p_cont_overlay <- df_ve_curves %>%
  filter(modelo != "Constante") %>%
  ggplot(aes(tau_m, VE, color = modelo)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  # Línea del modelo constante como referencia
  geom_hline(
    yintercept = 1 - exp(coef(fit_base)["inmunizado"]),
    linetype   = "dotted", color = "grey50", linewidth = 0.5
  ) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = model_colors) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     breaks = seq(-0.5, 1, by = 0.25)) +
  scale_x_continuous(breaks = seq(0, 15, by = 3)) +
  coord_cartesian(xlim = c(0, 15), ylim = c(-0.3, 1.1)) +
  labs(
    x        = "Months since immunization",
    y        = "VE(τ) = 1 − HR(τ)",
    color    = "Functional form",
    title    = "Vaccine effectiveness — continuous time-varying coefficient",
    subtitle = "Dotted line = constant VE (PH model)"
  ) +
  theme_classic(base_size = 12, base_family = "serif") +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "grey35", size = 9.5),
    legend.position = "bottom",
    aspect.ratio  = 0.55
  )

# Panel B: facetas con IC
p_cont_facet <- df_ve_curves %>%
  filter(modelo != "Constante") %>%
  ggplot(aes(tau_m, VE, color = modelo, fill = modelo)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  geom_ribbon(aes(ymin = VE_lo, ymax = VE_hi),
              alpha = 0.12, linewidth = 0) +
  geom_line(linewidth = 0.95) +
  facet_wrap(~ modelo, ncol = 2) +
  scale_color_manual(values = model_colors) +
  scale_fill_manual(values  = model_colors) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     breaks = seq(-0.5, 1, by = 0.25)) +
  scale_x_continuous(breaks = seq(0, 15, by = 5)) +
  coord_cartesian(xlim = c(0, 15), ylim = c(-0.3, 1.1)) +
  labs(
    x        = "Months since immunization",
    y        = "VE(τ) = 1 − HR(τ)",
    title    = "Vaccine effectiveness — continuous forms with 95% CI",
    subtitle = "Shaded area: pointwise 95% CI · strata(Group) + cluster(RUN)"
  ) +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(color = "grey40", size = 9),
    legend.position  = "none",
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    strip.text       = element_text(face = "bold", size = 9.5),
    panel.grid.minor = element_blank()
  )

print(p_cont_overlay)
print(p_cont_facet)

ggsave("ve_continuous_overlay.pdf", p_cont_overlay,
       width = 7, height = 4.5, device = cairo_pdf)
ggsave("ve_continuous_facet.pdf",   p_cont_facet,
       width = 7.5, height = 6,   device = cairo_pdf)

# ============================================================
# 9) FIGURA 4: Combinar step + continuo en un panel
#    (figura de resumen, paper-ready)
# ============================================================

# Añadir step function como escalones al mismo gráfico
df_step_rect <- df_step_ve %>%
  filter(!is.na(tau_mid)) %>%
  mutate(
    x_lo = c(0, cut_points_m)[seq_len(n())],
    x_hi = c(cut_points_m, 15)[seq_len(n())]
  )

p_summary <- ggplot() +
  # Step function como rectángulos sombreados
  geom_rect(
    data = df_step_rect,
    aes(xmin = x_lo, xmax = x_hi, ymin = VE_lo, ymax = VE_hi),
    fill = "grey85", alpha = 0.6, color = NA
  ) +
  geom_segment(
    data = df_step_rect,
    aes(x = x_lo, xend = x_hi, y = VE, yend = VE),
    color = "grey50", linewidth = 1.0
  ) +
  # Curvas continuas encima
  geom_line(
    data = df_ve_curves %>%
      filter(modelo %in% c("Log(tau+1)", "Sqrt(tau)")),
    aes(tau_m, VE, color = modelo),
    linewidth = 1.0
  ) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  scale_color_manual(
    values = c("Log(tau+1)" = "#1B4F72", "Sqrt(tau)" = "#7D3C98"),
    labels = c("Log(tau+1)" = "Logarítmica", "Sqrt(tau)" = "Raíz cuadrada")
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     breaks = seq(-0.25, 1, by = 0.25)) +
  scale_x_continuous(breaks = seq(0, 15, by = 3)) +
  coord_cartesian(xlim = c(0, 15), ylim = c(-0.2, 1.05)) +
  labs(
    x        = "Months since immunization",
    y        = "Vaccine effectiveness  VE(τ) = 1 − HR(τ)",
    color    = "Continuous model",
    title    = "Waning of vaccine effectiveness",
    subtitle = "Grey: step function estimate ± 95% CI  ·  Lines: continuous models",
    caption  = "Step function: bins of 3 months. Continuous: time-varying coefficient Cox model."
  ) +
  theme_classic(base_size = 12, base_family = "serif") +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "grey35", size = 9.5),
    plot.caption  = element_text(color = "grey45", size = 8, hjust = 0),
    legend.position = "bottom",
    legend.text     = element_text(size = 10),
    aspect.ratio    = 0.55
  )

print(p_summary)
ggsave("ve_summary_paper.pdf", p_summary,
       width = 7, height = 4.5, device = cairo_pdf)
ggsave("ve_summary_paper.png", p_summary,
       width = 7, height = 4.5, dpi = 320)

# ============================================================
# 11) TABLA RESUMEN — todos los modelos
# ============================================================

tab_models <- data.frame(
  Modelo    = c("PH constante", "Log(τ+1)", "Lineal", "Cuadrático", "Sqrt(τ)"),
  log_lik   = round(c(
    fit_base$loglik[2],
    fit_log$loglik[2],
    fit_linear$loglik[2],
    fit_quad$loglik[2],
    fit_sqrt$loglik[2]
  ), 2),
  df        = c(1, 2, 2, 3, 2),
  LRT_vs_PH = round(2 * c(
    0,
    fit_log$loglik[2]    - fit_base$loglik[2],
    fit_linear$loglik[2] - fit_base$loglik[2],
    fit_quad$loglik[2]   - fit_base$loglik[2],
    fit_sqrt$loglik[2]   - fit_base$loglik[2]
  ), 3)
) %>%
  mutate(
    AIC      = round(-2 * log_lik + 2 * df, 2),
    delta_AIC = round(AIC - min(AIC), 2),
    p_vs_PH  = c(NA, round(pchisq(LRT_vs_PH[-1], df = 1,
                                   lower.tail = FALSE), 4))
  )

cat("\n══ Tabla comparativa de modelos ══\n")
print(tab_models)

# Guardar como HTML
library(knitr)
library(kableExtra)

html_tab <- tab_models %>%
  kable("html",
        caption = "Model comparison — time-varying VE(τ)",
        col.names = c("Model", "Log-lik", "df",
                      "LRT vs PH", "AIC", "ΔAIC", "p vs PH"),
        align = c("l", rep("r", 6))) %>%
  kable_styling(full_width = FALSE,
                bootstrap_options = c("striped","hover","condensed"),
                font_size = 13) %>%
  row_spec(which.min(tab_models$AIC),
           bold = TRUE, background = "#D6EAF8") %>%
  row_spec(0, bold = TRUE, background = "#2C3E50", color = "white") %>%
  as.character()

writeLines(paste0(
  '<!DOCTYPE html><html><head><meta charset="utf-8">',
  '<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css">',
  '<style>body{padding:30px;font-family:Georgia,serif;}</style>',
  '</head><body>', html_tab, '</body></html>'
), "model_comparison_table.html")

cat("✓ Tabla guardada en model_comparison_table.html\n")