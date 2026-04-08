# ============================================================
# JENNINGS-STYLE FPM FOR RSV WANING
# faithful to:
# - backwards/reversed restricted cubic spline basis
# - 2-step constrained fit
# - log-time baseline
# - log(1+tau) treatment-effect spline
#
# NOTES:
# 1) This script DOES NOT force the early 0-5 month shape.
#    It only gives more flexibility there via early treatment knots.
# 2) It DOES force HR=1 from tau_star_m onward in the constrained fit.
# 3) It does NOT include bootstrap yet.
# ============================================================

rm(list = ls())

library(readr)
library(dplyr)
library(tidyr)
library(survival)
library(survPen)
library(ggplot2)

# ============================================================
# 0) USER OPTIONS
# ============================================================

setwd("C:/Users/ntrig/Desktop/ISCI/Proyectos/Efectividad_Nirse/Code")

data_file <- "../Data/df_tv_matched_11_03_2026.csv"

# split width for immunized intervals
split_width_days <- 14

# sensitivity parameter:
# choose the month from which HR is forced to 1
# start with 8, 10, or 12 as sensible sensitivity values.
tau_star_m <- 15

# number of total baseline knots (includes boundaries)
baseline_nk <- 3

# treatment knots for unconstrained model, on raw month scale.
# these are EARLY by design, to let the model flex in 0-5 months.
tau_internal_un_m <- c(1, 3, 5, 8)

# treatment knots for constrained model, on raw month scale.
# these are clipped to be < tau_star_m automatically.
tau_internal_co_m <- c(1, 5, 10)

# ============================================================
# 1) DOWNLOAD/SOURCE JENNINGS CONSTRAINED survPen CODE
# ============================================================

# If you already saved ConstrainNR.R locally, replace these 3 lines with:
# source("C:/path/to/ConstrainNR.R")

constrain_url  <- "https://raw.githubusercontent.com/angusjennings/spline-model/main/ConstrainNR.R"
constrain_file <- file.path(tempdir(), "ConstrainNR.R")
download.file(constrain_url, constrain_file, mode = "wb")
source(constrain_file)

# ============================================================
# 2) FUNCTIONS FROM / ALIGNED WITH THE PAPER
# ============================================================

# x if x>0 or 0 if x<=0
pos <- function(x) {
  sapply(x, function(z) if (z > 0) z else 0)
}

# reversed restricted cubic spline
# same spirit as the paper's R code
RevResCubicSpline <- function(x, k = NULL, nk, int = FALSE, pre = NULL, drop = FALSE) {
  if (is.null(k)) {
    k <- quantile(x, seq(0, 1, length.out = nk), na.rm = TRUE)
  }
  if (!is.null(pre)) {
    pre <- paste0(pre, ".")
  }

  K <- length(k)

  s <- data.frame(s1 = x)

  if (int) s <- cbind(1, s)

  if (K > 1) {
    for (j in 1:(K - 2)) {
      l <- (k[K - j] - k[1]) / (k[K] - k[1])

      s[[paste0("s", j)]] <-
        pos((k[K - j] - x)^3) -
        l * pos((k[K] - x)^3) -
        (1 - l) * pos((k[1] - x)^3)
    }

    s[[paste0("s", K - 1)]] <- x
  }

  s <- as.matrix(s)
  colnames(s) <- paste0(pre, colnames(s))

  if (drop) {
    return(s[, -ncol(s), drop = FALSE])
  } else {
    return(s)
  }
}

# ============================================================
# 3) DATA PREP
# ============================================================

df_tv_matched <- read_csv(
  data_file,
  col_types = cols(
    DIAG9  = col_character(),
    DIAG10 = col_character(),
    DIAG11 = col_character()
  ),
  show_col_types = FALSE
)

df_raw <- df_tv_matched

need <- c("RUN", "start", "stop", "event_vrs", "Group", "inmunizado",'year_nac','mes_nac_name','sexo','region','SEMANAS')
stopifnot(all(need %in% names(df_raw)))

df_raw$RUN        <- as.character(df_raw$RUN)
df_raw$Group      <- as.factor(df_raw$Group)
df_raw$start      <- as.numeric(df_raw$start)
df_raw$stop       <- as.numeric(df_raw$stop)
df_raw$event_vrs  <- as.integer(df_raw$event_vrs)
df_raw$inmunizado <- as.integer(df_raw$inmunizado)

df_raw <- df_raw %>% filter(stop > start)

stopifnot(all(df_raw$event_vrs %in% c(0, 1)))
stopifnot(all(df_raw$inmunizado %in% c(0, 1)))

# same time origin as your previous work
if (!"t_inm" %in% names(df_raw)) {
  ref_date <- as.Date("2024-04-01")

  if (!"fechaInm" %in% names(df_raw)) {
    stop("Need either t_inm or fechaInm in the dataset.")
  }

  df_raw$fechaInm <- as.Date(df_raw$fechaInm)
  df_raw$t_inm    <- as.numeric(df_raw$fechaInm - ref_date)
}

df_raw$t_inm <- as.numeric(df_raw$t_inm)

if (any(df_raw$inmunizado == 1 & is.na(df_raw$t_inm))) {
  stop("There are rows with inmunizado=1 but missing t_inm/fechaInm.")
}

# ============================================================
# 4) INTERVAL REBUILDING
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
      Group         = first(Group),
      inmunizado    = first(inmunizado),
      t_inm         = first(t_inm),
      year_nac      = first(year_nac),
      mes_nac_name  = first(mes_nac_name),
      sexo          = first(sexo),
      SEMANAS       = first(SEMANAS),
      region        = first(region),
      start         = first(start),
      stop          = last(stop),
      event_vrs     = as.integer(last(event_vrs)),
      .groups       = "drop"
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
  ) %>% filter(stop > start)

  out$event_vrs <- 0L
  if (nrow(out) > 0) out$event_vrs[nrow(out)] <- as.integer(event)

  out
}

rebuild_intervals <- function(data, split_width_days = 14) {
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

        tmp$RUN          <- row_i$RUN
        tmp$Group        <- row_i$Group
        tmp$inmunizado   <- row_i$inmunizado
        tmp$t_inm        <- row_i$t_inm
        tmp$year_nac     <- row_i$year_nac
        tmp$mes_nac_name <- row_i$mes_nac_name
        tmp$sexo         <- row_i$sexo
        tmp$SEMANAS      <- row_i$SEMANAS
        tmp$region       <- row_i$region

        split_list[[i]] <- tmp[, c(
            "RUN", "Group", "inmunizado", "t_inm",
            "year_nac", "mes_nac_name", "sexo", "SEMANAS", "region",
            "start", "stop", "event_vrs"
        )]
        } else {
            split_list[[i]] <- row_i[, c(
                "RUN", "Group", "inmunizado", "t_inm",
                "year_nac", "mes_nac_name", "sexo", "SEMANAS", "region",
                "start", "stop", "event_vrs"
                )]
                }
                
  }

  bind_rows(split_list) %>%
    arrange(RUN, start, stop) %>%
    filter(stop >= start)
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

identify_forward_var <- function(knots_logtau, tau_star_m, prefix = "tc") {
  x_tail <- seq(log1p(tau_star_m) + 1e-6, log1p(tau_star_m + 6), length.out = 200)

  B_tail <- as.data.frame(
    RevResCubicSpline(
      x = x_tail,
      k = knots_logtau,
      nk = length(knots_logtau),
      int = FALSE,
      pre = prefix,
      drop = FALSE
    )
  )

  sds <- sapply(B_tail, sd)
  print(sds)

  names(which.max(sds))
}

df_j <- rebuild_intervals(df_raw, split_width_days = split_width_days) %>%
  add_tau_variables() %>%
  mutate(
    tau_m = tau_eval_m,
    t0_m  = start / 30.4375,
    t1_m  = stop  / 30.4375
  ) %>%
  filter(t1_m >= t0_m)

df_j <- df_j %>%
  mutate(
    year_nac     = as.numeric(year_nac),
    mes_nac_name = as.factor(mes_nac_name),
    sexo         = as.factor(sexo),
    SEMANAS      = as.numeric(SEMANAS),
    region       = as.factor(region)
  )

# ============================================================
# 5) KNOT HELPERS
# ============================================================

strict_inside <- function(x, lower, upper, eps = 1e-6) {
  x <- sort(unique(x))
  x <- x[x > lower + eps & x < upper - eps]
  x
}

make_baseline_knots_logt <- function(t_event_m, nk = 5) {
  x <- log(pmax(t_event_m, 1e-6))
  as.numeric(quantile(x, probs = seq(0, 1, length.out = nk), na.rm = TRUE, type = 7))
}

make_tau_knots_un <- function(tau_m, internal_months = c(1, 3, 5, 8)) {
  upper <- max(tau_m, na.rm = TRUE)
  raw <- c(0, strict_inside(internal_months, 0, upper), upper)
  log1p(raw)
}

make_tau_knots_co <- function(tau_m, tau_star_m, internal_months = c(1, 3, 5)) {
  tau_pre <- tau_m[is.finite(tau_m) & tau_m > 0 & tau_m < tau_star_m]

  if (length(tau_pre) >= 10) {
    p95 <- as.numeric(quantile(tau_pre, 0.95, na.rm = TRUE, type = 7))
    ints <- c(internal_months, p95)
  } else {
    ints <- internal_months
  }

  raw <- c(0, strict_inside(ints, 0, tau_star_m), tau_star_m)
  log1p(raw)
}

# ============================================================
# 6) DESIGN MATRIX HELPERS
# ============================================================

add_rev_spline_cols <- function(df, x, knots, prefix) {
  B <- RevResCubicSpline(
    x = x,
    k = knots,
    nk = length(knots),
    int = FALSE,
    pre = prefix,
    drop = FALSE
  )
  bind_cols(df, as.data.frame(B))
}

add_trt_rev_spline_cols <- function(df, x, knots, treat_var = "inmunizado", prefix = "trt") {
  B <- RevResCubicSpline(
    x = x,
    k = knots,
    nk = length(knots),
    int = FALSE,
    pre = prefix,
    drop = FALSE
  )
  B <- as.data.frame(B)
  for (nm in names(B)) {
    df[[nm]] <- df[[treat_var]] * B[[nm]]
  }
  df
}

make_formula <- function(base_names, trt_names, adj_terms = NULL, trt_main = "inmunizado") {
  rhs <- c(base_names, adj_terms, trt_main, trt_names)
  as.formula(paste("~", paste(rhs, collapse = " + ")))
}

mode_value <- function(x) {
  ux <- unique(x[!is.na(x)])
  ux[which.max(tabulate(match(x, ux)))]
}

build_reference_profile <- function(df) {
  list(
    year_nac     = median(df$year_nac, na.rm = TRUE),
    mes_nac_name = mode_value(df$mes_nac_name),
    sexo         = mode_value(df$sexo),
    SEMANAS      = median(df$SEMANAS, na.rm = TRUE),
    region       = mode_value(df$region)
  )
}

get_adjustment_terms <- function(df) {
  terms <- c()

  if ("year_nac" %in% names(df) && sum(is.finite(df$year_nac)) > 1 &&
      sd(df$year_nac, na.rm = TRUE) > 0) {
    terms <- c(terms, "year_nac")
  }

  if ("mes_nac_name" %in% names(df) && nlevels(droplevels(df$mes_nac_name)) > 1) {
    terms <- c(terms, "mes_nac_name")
  }

  if ("sexo" %in% names(df) && nlevels(droplevels(df$sexo)) > 1) {
    terms <- c(terms, "sexo")
  }

  if ("SEMANAS" %in% names(df) && sum(is.finite(df$SEMANAS)) > 1 &&
      sd(df$SEMANAS, na.rm = TRUE) > 0) {
    terms <- c(terms, "SEMANAS")
  }

  if ("region" %in% names(df) && nlevels(droplevels(df$region)) > 1) {
    terms <- c(terms, "region")
  }

  terms
}

# ============================================================
# 7) PREDICTION HELPER
# ============================================================

predict_hr_curve_survpen_rev <- function(
    fit,
    tau_grid_m,
    t_pred_m,
    k_base_logt,
    k_trt_logtau,
    trt_prefix,
    ref_profile,
    trt_main = "inmunizado",
    constrained = FALSE,
    tau_star_m = NULL
) {

    nd <- data.frame(
        t1_m = rep(t_pred_m, length(tau_grid_m)),
        tau_m = tau_grid_m,
        inmunizado = 1,
        year_nac = ref_profile$year_nac,
        mes_nac_name = factor(ref_profile$mes_nac_name, levels = levels(df_j$mes_nac_name)),
        sexo = factor(ref_profile$sexo, levels = levels(df_j$sexo)),
        SEMANAS = ref_profile$SEMANAS,
        region = factor(ref_profile$region, levels = levels(df_j$region))
        )

  nd <- add_rev_spline_cols(
    nd,
    x = log(pmax(nd$t1_m, 1e-6)),
    knots = k_base_logt,
    prefix = "bh"
  )

  nd <- add_trt_rev_spline_cols(
    nd,
    x = log1p(nd$tau_m),
    knots = k_trt_logtau,
    treat_var = "inmunizado",
    prefix = trt_prefix
  )

  nd_ref <- nd
  nd_ref[[trt_main]] <- 0
  trt_cols <- grep(paste0("^", trt_prefix, "\\."), names(nd_ref), value = TRUE)
  nd_ref[, trt_cols] <- 0

  pr <- predict(fit, newdata = nd, newdata.ref = nd_ref, type = "HR")
  pr <- as.data.frame(pr)

  out <- data.frame(
    tau_m = tau_grid_m,
    HR    = pr$HR,
    HR_lo = pr$HR.inf,
    HR_hi = pr$HR.sup,
    VE    = 1 - pr$HR,
    VE_lo = 1 - pr$HR.sup,
    VE_hi = 1 - pr$HR.inf
  )

  if (isTRUE(constrained) && !is.null(tau_star_m)) {
    idx <- out$tau_m >= tau_star_m
    out$HR[idx]    <- 1
    out$HR_lo[idx] <- 1
    out$HR_hi[idx] <- 1
    out$VE[idx]    <- 0
    out$VE_lo[idx] <- 0
    out$VE_hi[idx] <- 0
  }

  out
}

# ============================================================
# 8) CHOOSE KNOTS
# ============================================================

event_t_m <- df_j$t1_m[df_j$event_vrs == 1]
tau_imm_m <- df_j$tau_eval_m[df_j$inmunizado == 1]

k_base_logt <- make_baseline_knots_logt(event_t_m, nk = baseline_nk)
k_trt_un    <- make_tau_knots_un(tau_imm_m, internal_months = tau_internal_un_m)
k_trt_co    <- make_tau_knots_co(tau_imm_m, tau_star_m = tau_star_m, internal_months = tau_internal_co_m)

forward_var_co <- identify_forward_var(
  knots_logtau = k_trt_co,
  tau_star_m   = tau_star_m,
  prefix       = "tc"
)

cat("\nForward extrapolation variable detected:\n")
print(forward_var_co)

cat("\nBaseline knots (log t):\n")
print(k_base_logt)

cat("\nTreatment unconstrained knots (log(1+tau)):\n")
print(k_trt_un)

cat("\nTreatment constrained knots (log(1+tau)):\n")
print(k_trt_co)

# ============================================================
# 9) STEP 1: UNCONSTRAINED JENNINGS-STYLE MODEL
# ============================================================

adj_terms <- get_adjustment_terms(df_j)
print(adj_terms)

ref_profile <- build_reference_profile(df_j)
print(ref_profile)

df_un <- df_j %>%
  add_rev_spline_cols(
    x = log(pmax(.$t1_m, 1e-6)),
    knots = k_base_logt,
    prefix = "bh"
  ) %>%
  add_trt_rev_spline_cols(
    x = log1p(.$tau_m),
    knots = k_trt_un,
    treat_var = "inmunizado",
    prefix = "tu"
  )

base_names_un <- grep("^bh\\.", names(df_un), value = TRUE)
trt_names_un  <- grep("^tu\\.", names(df_un), value = TRUE)



form_un <- make_formula(
  base_names = base_names_un,
  trt_names  = trt_names_un,
  adj_terms  = adj_terms,
  trt_main   = "inmunizado"
)

fit_un <- survPen_cons(
  formula = form_un,
  data    = df_un,
  t1      = t1_m,
  t0      = t0_m,
  event   = event_vrs == 1
)

cat("\nUnconstrained coefficients:\n")
print(fit_un$coefficients)

# ============================================================
# 10) STEP 2: CONSTRAINED MODEL
#     Fix:
#     - baseline coefficients from step 1
#     - treatment PH main effect = 0
#     - forward extrapolation spline variable = 0
# ============================================================

df_co <- df_j %>%
  add_rev_spline_cols(
    x = log(pmax(.$t1_m, 1e-6)),
    knots = k_base_logt,
    prefix = "bh"
  ) %>%
  add_trt_rev_spline_cols(
    x = log1p(.$tau_m),
    knots = k_trt_co,
    treat_var = "inmunizado",
    prefix = "tc"
  )

base_names_co <- grep("^bh\\.", names(df_co), value = TRUE)
trt_names_co  <- grep("^tc\\.", names(df_co), value = TRUE)

form_co <- make_formula(
  base_names = base_names_co,
  trt_names  = trt_names_co,
  adj_terms  = adj_terms,
  trt_main   = "inmunizado"
)

coef_names_co <- colnames(model.matrix(form_co, data = df_co))
beta_ini <- setNames(rep(0, length(coef_names_co)), coef_names_co)

baseline_fixed_names <- c("(Intercept)", base_names_co)
beta_ini[baseline_fixed_names] <- fit_un$coefficients[baseline_fixed_names]

# inicializar covariables comunes desde el unconstrained
common_names <- intersect(names(fit_un$coefficients), coef_names_co)
beta_ini[common_names] <- fit_un$coefficients[common_names]

# luego sobreescribir las constraints de tratamiento
zero_names <- c("inmunizado", forward_var_co)
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

cat("\nConstrained coefficients:\n")
print(fit_co$coefficients)

# ============================================================
# 11) PREDICT CURVES
# ============================================================

tau_grid_m <- seq(0, ceiling(max(df_j$tau_m[df_j$inmunizado == 1], na.rm = TRUE)), by = 0.05)

# choose a representative analysis time for conditional HR prediction
t_pred_m <- median(df_j$t1_m[df_j$event_vrs == 1], na.rm = TRUE)

curve_un <- predict_hr_curve_survpen_rev(
  fit          = fit_un,
  tau_grid_m   = tau_grid_m,
  t_pred_m     = t_pred_m,
  k_base_logt  = k_base_logt,
  k_trt_logtau = k_trt_un,
  trt_prefix   = "tu",
  ref_profile  = ref_profile,
  constrained  = FALSE
) %>%
  mutate(model = "Jennings unconstrained")

curve_co <- predict_hr_curve_survpen_rev(
  fit          = fit_co,
  tau_grid_m   = tau_grid_m,
  t_pred_m     = t_pred_m,
  k_base_logt  = k_base_logt,
  k_trt_logtau = k_trt_co,
  trt_prefix   = "tc",
  ref_profile  = ref_profile,
  constrained  = TRUE,
  tau_star_m   = tau_star_m
) %>%
  mutate(model = paste0("Jennings constrained (HR=1 from ", tau_star_m, " m)"))

curves_all <- bind_rows(curve_un, curve_co)

# ============================================================
# 12) DIAGNOSTIC PLOTS
# ============================================================

p_compare <- ggplot(curves_all, aes(x = tau_m, y = VE, color = model, fill = model)) +
  geom_ribbon(aes(ymin = VE_lo, ymax = VE_hi), alpha = 0.12, linewidth = 0) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey30") +
  geom_vline(xintercept = tau_star_m, linetype = 3) +
  theme_bw() +
  labs(
    x = "Months since immunization",
    y = "VE(t) = 1 - HR(t)",
    title = "Jennings-style FPM with reversed restricted cubic splines"
  )

print(p_compare)

# ============================================================
# 13) OPTIONAL: OVERLAY WITH YOUR COX FREE CURVE
#     only if you already created results[['w14_df3']]
# ============================================================

if (exists("results")) {
  if ("w14_df3" %in% names(results)) {
    cox_curve <- results[["w14_df3"]]$ve_curve %>%
      transmute(
        tau_m, VE, VE_lo, VE_hi,
        model = "Cox free spline (df=3, split=14)"
      )

    p_overlay <- bind_rows(
      cox_curve,
      curve_un %>% select(tau_m, VE, VE_lo, VE_hi, model),
      curve_co %>% select(tau_m, VE, VE_lo, VE_hi, model)
    ) %>%
      ggplot(aes(x = tau_m, y = VE, color = model, fill = model)) +
      geom_ribbon(aes(ymin = VE_lo, ymax = VE_hi), alpha = 0.08, linewidth = 0) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 0, linetype = 2, color = "grey30") +
      geom_vline(xintercept = tau_star_m, linetype = 3) +
      theme_bw() +
      labs(
        x = "Months since immunization",
        y = "VE(t) = 1 - HR(t)",
        title = "Cox vs Jennings-style FPM"
      )

    print(p_overlay)
  }
}

# ============================================================
# 14) EXPORT
# ============================================================

dir.create("../Output/jennings_rev_rcs", showWarnings = FALSE, recursive = TRUE)

write.csv(curve_un, file = "../Output/jennings_rev_rcs/curve_unconstrained.csv", row.names = FALSE)
write.csv(curve_co, file = "../Output/jennings_rev_rcs/curve_constrained.csv", row.names = FALSE)

ggsave(
  filename = "../Output/jennings_rev_rcs/jennings_compare.png",
  plot = p_compare,
  width = 9,
  height = 6,
  dpi = 300
)

cat("\nDone.\n")
cat("tau_star_m =", tau_star_m, "\n")
cat("Suggested next sensitivity values: 8, 10, 12 months.\n")