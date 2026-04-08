

# ANTES DE CORRER ESTE ARCHIVO SE DEBE CORRER DECAY_CURVE.R #############

library(survival)
library(splines)
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)

setwd("C:/Users/ntrig/Desktop/ISCI/Proyectos/Efectividad_Nirse/Code/decay_r")

# ============================================================
# A) ELEGIR KNOTS INTERNOS PARA EL EFECTO RESTRINGIDO
#    - mediana
#    - percentil 95 (estilo Jennings)
# ============================================================
choose_hr1_knots <- function(tau_obs, tau_star_m, probs = c(0.25, 0.50, 0.75, 0.95)) {
  x <- tau_obs[is.finite(tau_obs) & tau_obs > 0 & tau_obs < tau_star_m]

  if (length(x) < 20) {
    stop("Muy poco soporte antes de tau_star para construir knots internos.")
  }

  k <- as.numeric(quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 7))
  k <- sort(unique(round(k, 6)))
  k <- k[k > 0 & k < tau_star_m]

  if (length(k) < 2) {
    # fallback mínimo
    k <- as.numeric(quantile(x, probs = c(0.33, 0.67), na.rm = TRUE, names = FALSE, type = 7))
    k <- sort(unique(round(k, 6)))
    k <- k[k > 0 & k < tau_star_m]
  }

  if (length(k) == 0) {
    stop("No se pudieron construir knots internos válidos.")
  }

  k
}

# ============================================================
# B) OBJETO DE CONFIGURACIÓN DE LA BASE RESTRINGIDA
# ============================================================
make_hr1_basis_setup <- function(tau_obs, tau_star_m, knot_probs = c(0.25, 0.50, 0.75, 0.95)) {
  knots <- choose_hr1_knots(
    tau_obs     = tau_obs,
    tau_star_m  = tau_star_m,
    probs       = knot_probs
  )

  list(
    tau_star_m = tau_star_m,
    knots      = knots,
    bknots     = c(0, tau_star_m),
    h          = max(1e-6, tau_star_m * 1e-6)
  )
}

# ============================================================
# C) EVALUAR BASE RESTRINGIDA:
#    C_j(tau) = B_j(tau) - B_j(tau*) - B_j'(tau*)(tau-tau*)
#    => C_j(tau)=0 para tau >= tau*
# ============================================================
eval_hr1_basis <- function(x, setup) {
  x <- as.numeric(x)

  basis_fun <- function(z) {
    ns(
      z,
      knots = setup$knots,
      Boundary.knots = setup$bknots,
      intercept = FALSE
    )
  }

  B      <- basis_fun(x)
  B_star <- basis_fun(setup$tau_star_m)
  B_plus <- basis_fun(setup$tau_star_m + setup$h)
  B_minus<- basis_fun(setup$tau_star_m - setup$h)

  dB_star <- (B_plus - B_minus) / (2 * setup$h)

  C <- B
  for (j in seq_len(ncol(B))) {
    C[, j] <- B[, j] - B_star[1, j] - dB_star[1, j] * (x - setup$tau_star_m)
  }

  # Por estabilidad numérica: exactamente 0 desde tau_star en adelante
  C[x >= setup$tau_star_m, ] <- 0

  colnames(C) <- paste0("imm_hr1_", seq_len(ncol(C)))
  attr(C, "B_star")  <- B_star
  attr(C, "dB_star") <- dB_star
  attr(C, "setup")   <- setup
  C
}

get_full_rank_keep <- function(M, tol = 1e-10) {
  q <- qr(M, tol = tol)
  keep <- sort(q$pivot[seq_len(q$rank)])
  keep
}

# ============================================================
# D) AJUSTE DEL COX RESTRINGIDO HR(t)=1 DESDE tau_star
# ============================================================
fit_waning_model_hr1 <- function(df_prepared, tau_star_m = 15, tau_grid_by = 0.05) {
  df <- df_prepared

  tau_imm <- df$tau_eval_m[df$inmunizado == 1]
  tau_imm <- tau_imm[is.finite(tau_imm)]

  setup <- make_hr1_basis_setup(
    tau_obs    = tau_imm,
    tau_star_m = tau_star_m,
    knot_probs = c(0.25, 0.50, 0.75, 0.95)
  )

  # Base restringida "raw"
  C_all_raw <- eval_hr1_basis(df$tau_eval_m, setup = setup)

  # Quedarse sólo con columnas linealmente independientes
  rank_rows <- df$inmunizado == 1 & is.finite(df$tau_eval_m) & df$tau_eval_m < tau_star_m

  keep_cols <- get_full_rank_keep(C_all_raw[rank_rows, , drop = FALSE])

  C_all <- C_all_raw[, keep_cols, drop = FALSE]
  colnames(C_all) <- paste0("imm_hr1_", seq_len(ncol(C_all)))

  # Agregar interacción con inmunización
  for (j in seq_len(ncol(C_all))) {
    df[[colnames(C_all)[j]]] <- df$inmunizado * C_all[, j]
  }

  rhs <- paste(
    c(
      colnames(C_all),
      "strata(Group)",
      "cluster(RUN)"
    ),
    collapse = " + "
  )

  f_hr1 <- as.formula(
    paste0("Surv(start, stop, event_vrs) ~ ", rhs)
  )

  fit_hr1 <- coxph(
    formula = f_hr1,
    data    = df,
    ties    = "efron",
    robust  = TRUE,
    x       = TRUE,
    model   = TRUE
  )

  tau_grid_m <- seq(
    0,
    ceiling(max(df$tau_eval_m[df$inmunizado == 1], na.rm = TRUE)),
    by = tau_grid_by
  )

  C_grid_raw <- eval_hr1_basis(tau_grid_m, setup = setup)
  C_grid <- C_grid_raw[, keep_cols, drop = FALSE]
  colnames(C_grid) <- colnames(C_all)

  Xg <- C_grid

  # Seguridad extra: dropear coeficientes NA si coxph dejó alguno
  beta_full <- coef(fit_hr1)
  keep_beta <- !is.na(beta_full)

  beta <- beta_full[keep_beta]
  V    <- vcov(fit_hr1)[keep_beta, keep_beta, drop = FALSE]
  Xg   <- Xg[, names(beta), drop = FALSE]

  eta    <- as.vector(Xg %*% beta)
  se_eta <- sqrt(rowSums((Xg %*% V) * Xg))

  HR    <- exp(eta)
  HR_lo <- exp(eta - 1.96 * se_eta)
  HR_hi <- exp(eta + 1.96 * se_eta)

  VE    <- 1 - HR
  VE_lo <- 1 - HR_hi
  VE_hi <- 1 - HR_lo

  # Forzar exactitud numérica desde tau_star
  idx_tail <- tau_grid_m >= tau_star_m

  HR[idx_tail]    <- 1
  HR_lo[idx_tail] <- 1
  HR_hi[idx_tail] <- 1

  VE[idx_tail]    <- 0
  VE_lo[idx_tail] <- 0
  VE_hi[idx_tail] <- 0

  ve_curve <- data.frame(
    model   = paste0("Constrained HR->1 at ", tau_star_m, " months"),
    tau_m   = tau_grid_m,
    logHR   = eta,
    se      = se_eta,
    HR      = HR,
    HR_lo   = HR_lo,
    HR_hi   = HR_hi,
    VE      = VE,
    VE_lo   = VE_lo,
    VE_hi   = VE_hi
  )

  list(
    df_model = df,
    fit_hr1  = fit_hr1,
    ve_curve = ve_curve,
    setup    = setup,
    keep_cols = keep_cols
  )
}

# ============================================================
# E) AJUSTE CON RESTRICCIÓN HR=1 DESDE tau_star
# ============================================================
df14 <- prepared_by_width[["w14"]]

# Elige tau_star.
# Te recomiendo partir con 15 meses como sensibilidad estructurada.
tau_star_m <- 15

res_hr1 <- fit_waning_model_hr1(
  df_prepared = df14,
  tau_star_m  = tau_star_m,
  tau_grid_by = 0.05
)

summary(res_hr1$fit_hr1)

# Ver los knots usados
print(res_hr1$setup)


free_curve <- results[["w14_df3"]]$ve_curve %>%
  transmute(
    model = "Free spline (df=3, split=14)",
    tau_m, VE, VE_lo, VE_hi
  )

constrained_curve <- res_hr1$ve_curve %>%
  transmute(
    model, tau_m, VE, VE_lo, VE_hi
  )

compare_curves <- bind_rows(free_curve, constrained_curve)

p_compare <- ggplot(compare_curves, aes(x = tau_m, y = VE, color = model, fill = model)) +
  geom_ribbon(aes(ymin = VE_lo, ymax = VE_hi), alpha = 0.12, linewidth = 0) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = tau_star_m, linetype = 3) +
  labs(
    x = "Months since immunization",
    y = "VE(t) = 1 - HR(t)",
    title = paste0("Free vs constrained Cox waning curve (HR=1 from ", tau_star_m, " months)")
  ) +
  theme_bw()

print(p_compare)