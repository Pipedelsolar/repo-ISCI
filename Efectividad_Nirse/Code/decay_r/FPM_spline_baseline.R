
# ANTES DE CORRER ESTE ARCHIVO SE DEBE CORRER DECAY_CURVE.R #############

# install.packages(c("survPen", "splines2", "survival"))
setwd("C:/Users/ntrig/Desktop/ISCI/Proyectos/Efectividad_Nirse/Code")


library(survival)
library(survPen)
library(splines2)
library(dplyr)
library(ggplot2)
library(boot)

# ============================================================
# 0) DESCARGAR Y CARGAR EL CÓDIGO MODIFICADO DE JENNINGS
# ============================================================
constrain_url <- "https://raw.githubusercontent.com/angusjennings/spline-model/main/ConstrainNR.R"
constrain_file <- file.path(tempdir(), "ConstrainNR.R")
download.file(constrain_url, constrain_file, mode = "wb")
source(constrain_file)

# ============================================================
# 1) FUNCIÓN naturalSpline2 "backwards" DEL EJEMPLO DE JENNINGS
# ============================================================
naturalSpline2 <- function(x, log = FALSE, df = NULL, knots = NULL,
                           intercept = FALSE, Boundary.knots = NULL,
                           derivs = 0L, integral = FALSE, ...) {

  if (log) {
    x <- log(x)
    if (!is.null(knots)) knots <- log(knots)
    if (!is.null(Boundary.knots)) Boundary.knots <- log(Boundary.knots)
  }

  x <- -x
  if (!is.null(knots)) knots <- -knots
  if (!is.null(Boundary.knots)) Boundary.knots <- -Boundary.knots

  naturalSpline(
    x,
    df = df,
    knots = knots,
    intercept = intercept,
    Boundary.knots = Boundary.knots,
    derivs = derivs,
    integral = integral,
    ...
  )
}

# ============================================================
# 2) HELPERS
# ============================================================
safe_quantiles <- function(x, probs) {
  q <- as.numeric(quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 7))
  q <- sort(q)
  for (i in 2:length(q)) {
    if (q[i] <= q[i - 1]) q[i] <- q[i - 1] + 1e-6
  }
  q
}

add_basis_cols <- function(df, x, all_knots, prefix) {
  B <- as.matrix(
    naturalSpline2(
      x = x,
      Boundary.knots = all_knots[c(1, length(all_knots))],
      knots = all_knots[-c(1, length(all_knots))],
      intercept = FALSE
    )
  )
  colnames(B) <- paste0(prefix, seq_len(ncol(B)))
  bind_cols(df, as.data.frame(B))
}

add_treatment_basis <- function(df, tau, all_knots, treat_var = "inmunizado", prefix = "trt_") {
  B <- as.matrix(
    naturalSpline2(
      x = tau,
      Boundary.knots = all_knots[c(1, length(all_knots))],
      knots = all_knots[-c(1, length(all_knots))],
      intercept = FALSE
    )
  )
  out_names <- paste0(prefix, seq_len(ncol(B)))
  for (j in seq_len(ncol(B))) {
    df[[out_names[j]]] <- df[[treat_var]] * B[, j]
  }
  attr(df, "basis_names") <- out_names
  df
}

make_formula_from_names <- function(base_names, trt_names, extra_names = NULL, trt_main = "inmunizado") {
  rhs <- c(base_names, trt_main, trt_names, extra_names)
  as.formula(paste("~", paste(rhs, collapse = " + ")))
}

predict_hr_curve_survpen <- function(fit, tau_grid_m, t_base_m, k_base, k_trt, trt_prefix, trt_main = "inmunizado") {

  nd <- data.frame(
    t1_m = rep(t_base_m, length(tau_grid_m)),
    tau_m = tau_grid_m,
    inmunizado = 1
  )

  nd <- add_basis_cols(nd, x = nd$t1_m, all_knots = k_base, prefix = "bh_")
  nd <- add_treatment_basis(nd, tau = nd$tau_m, all_knots = k_trt, treat_var = "inmunizado", prefix = trt_prefix)

  nd_ref <- nd
  nd_ref[[trt_main]] <- 0
  trt_cols <- grep(paste0("^", trt_prefix), names(nd_ref), value = TRUE)
  nd_ref[, trt_cols] <- 0

  pr <- predict(fit, newdata = nd, newdata.ref = nd_ref, type = "HR")
  pr <- as.data.frame(pr)

  data.frame(
    tau_m = tau_grid_m,
    HR    = pr$HR,
    HR_lo = pr$HR.inf,
    HR_hi = pr$HR.sup,
    VE    = 1 - pr$HR,
    VE_lo = 1 - pr$HR.sup,
    VE_hi = 1 - pr$HR.inf
  )
}

make_all_knots <- function(x, lower, upper, internal_probs = c(0.33, 0.67, 0.95), eps = 1e-6) {
  x <- x[is.finite(x) & x > lower & x < upper]

  if (length(x) < 10) {
    stop("Muy poco soporte dentro de los boundary knots para construir knots internos.")
  }

  k_int <- as.numeric(quantile(x, probs = internal_probs, na.rm = TRUE, names = FALSE, type = 7))

  # Forzar que queden estrictamente dentro de (lower, upper)
  k_int <- pmin(pmax(k_int, lower + eps), upper - eps)

  # Forzar orden estricto
  k_int <- sort(k_int)
  if (length(k_int) >= 2) {
    for (i in 2:length(k_int)) {
      if (k_int[i] <= k_int[i - 1]) {
        k_int[i] <- min(upper - eps, k_int[i - 1] + eps)
      }
    }
  }

  k_int <- unique(round(k_int, 6))

  if (length(k_int) == 0) {
    stop("No se pudieron construir knots internos válidos.")
  }

  c(lower, k_int, upper)
}

# ============================================================
# 3) BASE DEL ANÁLISIS
# ============================================================
df_j <- prepared_by_width[["w30"]] %>%
  filter(stop > start) %>%
  mutate(
    t0_m    = start / 30, #.4375,
    t1_m    = stop  / 30, #.4375,
    t_inm_m = t_inm / 30, #.4375,
    tau_m   = tau_eval_m
  )

# si quieres limitar a seguimiento positivo real
df_j <- df_j %>% filter(t1_m > t0_m)

# punto donde fuerzas HR=1
tau_star_m <- 5

# número total de knots (incluye boundary knots), siguiendo el ejemplo
deg_base <- 5
deg_trt  <- 5

# ============================================================
# 4) KNOTS
# ============================================================
# baseline: cuantiles de tiempos de evento
event_t <- df_j$t1_m[df_j$event_vrs == 1]
k_base <- safe_quantiles(event_t, seq(0, 1, length.out = deg_base))

# ============================================================
# 4) KNOTS CORREGIDOS
# ============================================================

# Baseline: boundary = min/max tiempo de evento, internos en cuantiles intermedios
event_t <- df_j$t1_m[df_j$event_vrs == 1 & is.finite(df_j$t1_m)]

k_base <- make_all_knots(
  x = event_t,
  lower = min(event_t, na.rm = TRUE),
  upper = max(event_t, na.rm = TRUE),
  internal_probs = c(0.25, 0.50, 0.75)
)

# tau para inmunizados
tau_event <- df_j$tau_m[df_j$event_vrs == 1 & df_j$inmunizado == 1 & is.finite(df_j$tau_m)]
if (length(tau_event) < 10) {
  tau_event <- df_j$tau_m[df_j$inmunizado == 1 & is.finite(df_j$tau_m)]
}

# Unconstrained treatment spline
k_trt_un <- make_all_knots(
  x = tau_event,
  lower = 0,
  upper = max(tau_event, na.rm = TRUE),
  internal_probs = c(0.25, 0.50, 0.75)
)

# k_trt_un <- c(0, 1, 3, 6, 12, max(df_j$tau_m[df_j$inmunizado == 1], na.rm = TRUE))

# Constrained treatment spline:
# boundaries fijos en 0 y tau_star_m
# internos dentro de (0, tau_star_m), con knot tardío tipo p95
k_trt_cons <- make_all_knots(
  x = tau_event,
  lower = 0,
  upper = tau_star_m,
  internal_probs = c(0.33, 0.67, 0.95)
)

print(k_base)
print(k_trt_un)
print(k_trt_cons)

stopifnot(all(diff(k_base) > 0))
stopifnot(all(diff(k_trt_un) > 0))
stopifnot(all(diff(k_trt_cons) > 0))

# ============================================================
# 5) MODELO UNCONSTRAINED
# ============================================================
df_un <- df_j %>%
  add_basis_cols(x = .$t1_m, all_knots = k_base, prefix = "bh_") %>%
  add_treatment_basis(tau = .$tau_m, all_knots = k_trt_un, treat_var = "inmunizado", prefix = "tu_")

base_names_un <- grep("^bh_", names(df_un), value = TRUE)
trt_names_un  <- grep("^tu_", names(df_un), value = TRUE)

form_un <- make_formula_from_names(
  base_names = base_names_un,
  trt_names  = trt_names_un,
  extra_names = NULL,
  trt_main = "inmunizado"
)

fit_un <- survPen_cons(
  formula = form_un,
  data    = df_un,
  t1      = t1_m,
  t0      = t0_m,
  event   = event_vrs == 1
)

print(fit_un$coefficients)


# ============================================================
# 6) MODELO CONSTRAINED
# ============================================================
df_co <- df_j %>%
  add_basis_cols(x = .$t1_m, all_knots = k_base, prefix = "bh_") %>%
  add_treatment_basis(tau = .$tau_m, all_knots = k_trt_cons, treat_var = "inmunizado", prefix = "tc_")

base_names_co <- grep("^bh_", names(df_co), value = TRUE)
trt_names_co  <- grep("^tc_", names(df_co), value = TRUE)

form_co <- make_formula_from_names(
  base_names = base_names_co,
  trt_names  = trt_names_co,
  extra_names = NULL,
  trt_main = "inmunizado"
)

# nombres esperados de coeficientes en el modelo constrained
coef_names_co <- colnames(model.matrix(form_co, data = df_co))

beta_ini <- setNames(rep(0, length(coef_names_co)), coef_names_co)

# baseline fijo = intercepto + spline baseline
baseline_fixed_names <- c("(Intercept)", base_names_co)

# copiar desde el unconstrained
beta_un <- fit_un$coefficients
beta_ini[baseline_fixed_names] <- beta_un[baseline_fixed_names]

# fijar a 0 el PH main effect del tratamiento y la forward extrapolation variable
zero_names <- c("inmunizado", "tc_1")
beta_ini[zero_names] <- 1e-12

cons_idx <- match(c(baseline_fixed_names, zero_names), coef_names_co)
cons_idx <- cons_idx[!is.na(cons_idx)]

fit_co <- survPen_cons(
  formula  = form_co,
  data     = df_co,
  t1       = t1_m,
  t0       = t0_m,
  event    = event_vrs == 1,
  beta.ini = unname(beta_ini),
  cons     = cons_idx
)

print(fit_co$coefficients)


# ============================================================
# 7) CURVAS HR/VE
# ============================================================
tau_grid_m <- seq(0, ceiling(max(df_j$tau_m[df_j$inmunizado == 1], na.rm = TRUE)), by = 0.05)

# elegir un tiempo baseline fijo; como el HR cancela baseline, da lo mismo mientras sea igual en ambos grupos
t_base_pred <- median(df_j$t1_m[df_j$event_vrs == 1], na.rm = TRUE)

curve_un <- predict_hr_curve_survpen(
  fit       = fit_un,
  tau_grid_m = tau_grid_m,
  t_base_m  = t_base_pred,
  k_base    = k_base,
  k_trt     = k_trt_un,
  trt_prefix = "tu_"
) %>%
  mutate(model = "Jennings unconstrained")

curve_co <- predict_hr_curve_survpen(
  fit       = fit_co,
  tau_grid_m = tau_grid_m,
  t_base_m  = t_base_pred,
  k_base    = k_base,
  k_trt     = k_trt_cons,
  trt_prefix = "tc_"
) %>%
  mutate(model = paste0("Jennings constrained (HR=1 from ", tau_star_m, " m)"))

# por estabilidad visual, forzar exactitud desde tau_star
curve_co <- curve_co %>%
  mutate(
    HR    = ifelse(tau_m >= tau_star_m, 1, HR),
    HR_lo = ifelse(tau_m >= tau_star_m, 1, HR_lo),
    HR_hi = ifelse(tau_m >= tau_star_m, 1, HR_hi),
    VE    = ifelse(tau_m >= tau_star_m, 0, VE),
    VE_lo = ifelse(tau_m >= tau_star_m, 0, VE_lo),
    VE_hi = ifelse(tau_m >= tau_star_m, 0, VE_hi)
  )

curves_j <- bind_rows(curve_un, curve_co)

ggplot(curves_j, aes(x = tau_m, y = VE, color = model, fill = model)) +
  geom_ribbon(aes(ymin = VE_lo, ymax = VE_hi), alpha = 0.12, linewidth = 0) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = tau_star_m, linetype = 3) +
  theme_bw() +
  labs(
    x = "Months since immunization",
    y = "VE(t) = 1 - HR(t)",
    title = "Jennings-survPen adaptation"
  )


  # ============================================================
# 8) BOOTSTRAP POR GROUP
# ============================================================

# strict_interior_knots <- function(k, lower, upper, eps = 1e-4) {
#   k <- sort(k)
#   k <- pmin(pmax(k, lower + eps), upper - eps)

#   if (length(k) >= 2) {
#     for (i in 2:length(k)) {
#       if (k[i] <= k[i - 1]) {
#         k[i] <- k[i - 1] + eps
#       }
#     }
#   }

#   # si por los ajustes alguno quedó fuera por arriba, redistribuir
#   if (length(k) > 0 && max(k) >= upper) {
#     k <- seq(lower + eps, upper - eps, length.out = length(k) + 2)[-c(1, length(k) + 2)]
#   }

#   k
# }

# make_all_knots_boot <- function(x, lower, upper, n_internal = 3,
#                                 probs = NULL, eps = 1e-4) {

#   if (!is.finite(lower) || !is.finite(upper) || upper <= lower + 10 * eps) {
#     stop("Boundary knots inválidos.")
#   }

#   x <- x[is.finite(x) & x > lower + eps & x < upper - eps]

#   if (is.null(probs)) {
#     probs <- seq_len(n_internal) / (n_internal + 1)
#   }

#   if (length(x) >= max(10, n_internal + 2) && length(unique(x)) >= 2) {
#     k_int <- as.numeric(
#       quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 7)
#     )
#   } else {
#     # fallback determinístico si la réplica viene pobre
#     k_int <- seq(lower, upper, length.out = n_internal + 2)[-c(1, n_internal + 2)]
#   }

#   k_int <- strict_interior_knots(k_int, lower = lower, upper = upper, eps = eps)

#   c(lower, k_int, upper)
# }

# fit_jennings_once <- function(dat, tau_star_m = 4, deg_base = 5, deg_trt = 5, tau_grid_m) {

#   dat <- dat %>%
#     filter(stop > start) %>%
#     mutate(
#       t0_m    = start / 30.4375,
#       t1_m    = stop  / 30.4375,
#       t_inm_m = t_inm / 30.4375,
#       tau_m   = tau_eval_m
#     )

#   n_int_base <- deg_base - 2
#   n_int_trt  <- deg_trt  - 2

#   # -------------------------
#   # Baseline knots
#   # -------------------------
#   event_t <- dat$t1_m[dat$event_vrs == 1 & is.finite(dat$t1_m)]

#   if (length(event_t) < 10) {
#     stop("Muy pocos eventos para construir knots baseline en esta réplica.")
#   }

#   k_base <- make_all_knots_boot(
#     x = event_t,
#     lower = min(event_t, na.rm = TRUE),
#     upper = max(event_t, na.rm = TRUE),
#     n_internal = n_int_base,
#     probs = seq_len(n_int_base) / (n_int_base + 1)
#   )

#   # -------------------------
#   # Treatment knots: usar soporte inmunizado, no solo eventos
#   # -------------------------
#   tau_imm_all <- dat$tau_m[dat$inmunizado == 1 & is.finite(dat$tau_m)]
#   tau_imm_pre <- dat$tau_m[
#     dat$inmunizado == 1 &
#       is.finite(dat$tau_m) &
#       dat$tau_m > 0 &
#       dat$tau_m < tau_star_m
#   ]

#   if (length(tau_imm_all) < 10) {
#     stop("Muy poco soporte inmunizado en esta réplica.")
#   }

#   # unconstrained: boundary superior en máximo observado
#   upper_un <- max(tau_imm_all, na.rm = TRUE)
#   if (!is.finite(upper_un) || upper_un <= 0) {
#     stop("Boundary superior inválido para spline unconstrained.")
#   }

#   k_trt_un <- make_all_knots_boot(
#     x = tau_imm_all,
#     lower = 0,
#     upper = upper_un,
#     n_internal = n_int_trt,
#     probs = seq_len(n_int_trt) / (n_int_trt + 1)
#   )

#   # constrained: boundary superior fijo en tau_star_m
#   # usar p95 como knot tardío si hay 3 internos
#   probs_cons <- if (n_int_trt == 3) c(0.33, 0.67, 0.95) else seq_len(n_int_trt) / (n_int_trt + 1)

#   x_cons <- if (length(tau_imm_pre) >= 10) tau_imm_pre else tau_imm_all[tau_imm_all < tau_star_m]

#   if (length(x_cons) < 5) {
#     stop("Muy poco soporte antes de tau_star en esta réplica.")
#   }

#   k_trt_cons <- make_all_knots_boot(
#     x = x_cons,
#     lower = 0,
#     upper = tau_star_m,
#     n_internal = n_int_trt,
#     probs = probs_cons
#   )

#   # -------------------------
#   # Unconstrained fit
#   # -------------------------
#   df_un <- dat %>%
#     add_basis_cols(x = .$t1_m, all_knots = k_base, prefix = "bh_") %>%
#     add_treatment_basis(tau = .$tau_m, all_knots = k_trt_un, treat_var = "inmunizado", prefix = "tu_")

#   base_names_un <- grep("^bh_", names(df_un), value = TRUE)
#   trt_names_un  <- grep("^tu_", names(df_un), value = TRUE)

#   form_un <- make_formula_from_names(
#     base_names = base_names_un,
#     trt_names  = trt_names_un,
#     trt_main   = "inmunizado"
#   )

#   fit_un <- survPen_cons(
#     formula = form_un,
#     data    = df_un,
#     t1      = t1_m,
#     t0      = t0_m,
#     event   = event_vrs == 1
#   )

#   # -------------------------
#   # Constrained fit
#   # -------------------------
#   df_co <- dat %>%
#     add_basis_cols(x = .$t1_m, all_knots = k_base, prefix = "bh_") %>%
#     add_treatment_basis(tau = .$tau_m, all_knots = k_trt_cons, treat_var = "inmunizado", prefix = "tc_")

#   base_names_co <- grep("^bh_", names(df_co), value = TRUE)
#   trt_names_co  <- grep("^tc_", names(df_co), value = TRUE)

#   form_co <- make_formula_from_names(
#     base_names = base_names_co,
#     trt_names  = trt_names_co,
#     trt_main   = "inmunizado"
#   )

#   coef_names_co <- colnames(model.matrix(form_co, data = df_co))

#   beta_ini <- setNames(rep(0, length(coef_names_co)), coef_names_co)

#   baseline_fixed_names <- c("(Intercept)", base_names_co)
#   beta_ini[baseline_fixed_names] <- fit_un$coefficients[baseline_fixed_names]

#   beta_ini[c("inmunizado", "tc_1")] <- 1e-12

#   cons_idx <- match(c(baseline_fixed_names, "inmunizado", "tc_1"), coef_names_co)
#   cons_idx <- cons_idx[!is.na(cons_idx)]

#   fit_co <- survPen_cons(
#     formula  = form_co,
#     data     = df_co,
#     t1       = t1_m,
#     t0       = t0_m,
#     event    = event_vrs == 1,
#     beta.ini = unname(beta_ini),
#     cons     = cons_idx
#   )

#   # -------------------------
#   # Prediction
#   # -------------------------
#   t_base_pred <- median(dat$t1_m[dat$event_vrs == 1], na.rm = TRUE)

#   pred <- predict_hr_curve_survpen(
#     fit        = fit_co,
#     tau_grid_m = tau_grid_m,
#     t_base_m   = t_base_pred,
#     k_base     = k_base,
#     k_trt      = k_trt_cons,
#     trt_prefix = "tc_"
#   )

#   pred$HR[pred$tau_m >= tau_star_m] <- 1

#   log(pred$HR)
# }

# boot_fun_group <- function(groups, indices, data_full, tau_star_m, deg_base, deg_trt, tau_grid_m) {
#   sampled_groups <- groups[indices]

#   dat_b <- bind_rows(lapply(seq_along(sampled_groups), function(i) {
#     out <- data_full %>% filter(Group == sampled_groups[i])
#     out$boot_group_id <- i
#     out
#   }))

#   out <- tryCatch(
#     fit_jennings_once(
#       dat = dat_b,
#       tau_star_m = tau_star_m,
#       deg_base = deg_base,
#       deg_trt = deg_trt,
#       tau_grid_m = tau_grid_m
#     ),
#     error = function(e) {
#       message("Bootstrap replicate failed: ", conditionMessage(e))
#       rep(NA_real_, length(tau_grid_m))
#     }
#   )

#   out
# }

# set.seed(123)

# groups <- unique(df_j$Group)

# bt <- boot(
#   data = groups,
#   statistic = function(g, idx) boot_fun_group(
#     groups = g,
#     indices = idx,
#     data_full = prepared_by_width[["w14"]],
#     tau_star_m = tau_star_m,
#     deg_base = deg_base,
#     deg_trt = deg_trt,
#     tau_grid_m = tau_grid_m
#   ),
#   R = 1
# )
# sum(apply(bt$t, 1, function(x) all(is.na(x))))

# # percentiles bootstrap
# ci_mat <- t(apply(bt$t, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE))
# HR_boot_lo <- exp(ci_mat[, 1])
# HR_boot_hi <- exp(ci_mat[, 2])

# curve_co$HR_lo_boot <- HR_boot_lo
# curve_co$HR_hi_boot <- HR_boot_hi
# curve_co$VE_lo_boot <- 1 - HR_boot_hi
# curve_co$VE_hi_boot <- 1 - HR_boot_lo

# ggplot(curves_j, aes(x = tau_m, y = VE, color = model, fill = model)) +
#   geom_ribbon(aes(ymin = VE_lo_boot, ymax = VE_hi_boot), alpha = 0.12, linewidth = 0) +
#   geom_line(linewidth = 1) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   geom_vline(xintercept = tau_star_m, linetype = 3) +
#   theme_bw() +
#   labs(
#     x = "Months since immunization",
#     y = "VE(t) = 1 - HR(t)",
#     title = "Jennings-survPen adaptation"
#   )