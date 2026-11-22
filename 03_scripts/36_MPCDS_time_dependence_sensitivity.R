# ============================================================
# 36_MPCDS_time_dependence_sensitivity.R
# Reproduces time-dependence analyses for the 504-patient LUAD
# MPCDS survival cohort with proliferation and negative-ESTIMATE
# relative purity proxy covariates.
#
# Run AFTER 35_proliferation_purity_sensitivity.R has produced:
#   LUAD_MPCDS_adjustment_patient_level.csv
# Writes CSV results and an original B/W PNG to 04_results;
# writes coloured PNG/PDF to 05_figures/manuscript/individual/mpcds
# Exploratory apparent/development-cohort analysis; not CV.
# ============================================================

if (!requireNamespace("survival", quietly = TRUE)) {
  stop("Please install the survival package before running this script.")
}

# survSplit() parses the formula LHS and expects the unqualified
# call Surv(...), not survival::Surv(...).
library(survival)

project_root <- "C:/Users/ibrah/R-Lang/Masters_Thesis/MSc_Lung_Cancer_PCD"
output_dir <- file.path(project_root, "04_results", "proliferation_purity_sensitivity")
input_file <- file.path(output_dir, "LUAD_MPCDS_adjustment_patient_level.csv")
if (!file.exists(input_file)) stop("Required input missing: ", input_file)

dat <- read.csv(input_file, stringsAsFactors = FALSE)
required <- c("patient_id", "OS_time", "OS_status", "MPCDS_z",
              "proliferation_z", "relative_purity_proxy_z")
stopifnot(all(required %in% names(dat)), !anyDuplicated(dat$patient_id))
stopifnot(nrow(dat) == 504L, sum(dat$OS_status) == 182L,
          all(is.finite(dat$OS_time)), all(dat$OS_time > 0),
          all(dat$OS_status %in% c(0, 1)),
          all(is.finite(as.matrix(dat[c("MPCDS_z", "proliferation_z", "relative_purity_proxy_z")]))))

# 1. Fully adjusted proportional-hazards model and residual test
fit_ph <- survival::coxph(
  survival::Surv(OS_time, OS_status) ~
    MPCDS_z + proliferation_z + relative_purity_proxy_z,
  data = dat, ties = "efron", x = TRUE
)
ph <- survival::cox.zph(fit_ph, transform = "km")
ph_table <- data.frame(term = rownames(ph$table), ph$table,
                       row.names = NULL, check.names = FALSE)
write.csv(ph_table, file.path(output_dir, "LUAD_MPCDS_PH_diagnostics_km.csv"), row.names = FALSE)

# Preserve the original black-and-white diagnostic with the model outputs.
png(file.path(output_dir, "LUAD_MPCDS_Schoenfeld_diagnostic.png"),
    width = 1800, height = 1200, res = 200)
plot(ph, var = "MPCDS_z", resid = TRUE, se = TRUE,
     main = "LUAD MPCDS: Schoenfeld Residual Diagnostic")
abline(h = coef(fit_ph)["MPCDS_z"], lty = 2)
dev.off()

# Create separate, colour-coded figure files in the EXISTING figure structure.
# This is an exploratory diagnostic, not an assigned manuscript panel or
# supplementary figure number until the authors approve its inclusion.
figure_dir <- file.path(
  project_root, "05_figures", "manuscript", "individual", "mpcds"
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

plot_colored_schoenfeld <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mar = c(5.0, 5.2, 3.3, 1.2), family = "sans")

  plot(ph, var = "MPCDS_z", resid = TRUE, se = TRUE,
       col = c("#006D77", "#78B8BE"),
       lwd = c(3, 1.6), lty = c(1, 2),
       pch = 1, cex = 0.75,
       xlab = "Follow-up time (days)",
       ylab = "Estimated log hazard ratio for MPCDS",
       main = "LUAD MPCDS: Schoenfeld residual diagnostic")
  abline(h = unname(coef(fit_ph)["MPCDS_z"]),
         col = "#D97706", lty = 3, lwd = 2)
  legend("topright",
         legend = c("Smoothed MPCDS coefficient",
                    "Approximate +/- 2 SE",
                    "Constant-coefficient Cox model"),
         col = c("#006D77", "#78B8BE", "#D97706"),
         lwd = c(3, 1.6, 2), lty = c(1, 2, 3),
         bty = "n", cex = 0.75)
}

figure_png <- file.path(
  figure_dir, "LUAD_MPCDS_Schoenfeld_diagnostic_COLOUR.png"
)
figure_pdf <- file.path(
  figure_dir, "LUAD_MPCDS_Schoenfeld_diagnostic_COLOUR.pdf"
)

png(figure_png, width = 2000, height = 1400, res = 220, bg = "white")
plot_colored_schoenfeld()
dev.off()

pdf(figure_pdf, width = 9, height = 6.3, useDingbats = FALSE)
plot_colored_schoenfeld()
dev.off()

# 2. Exploratory log-time interaction model, time in days
fit_tv <- survival::coxph(
  survival::Surv(OS_time, OS_status) ~
    MPCDS_z + proliferation_z + relative_purity_proxy_z + tt(MPCDS_z),
  data = dat, ties = "efron",
  tt = function(x, t, ...) x * log(t / 365.25)
)
tv_coef <- summary(fit_tv)$coefficients
tv_results <- data.frame(term = rownames(tv_coef), tv_coef,
                         row.names = NULL, check.names = FALSE)
write.csv(tv_results, file.path(output_dir, "LUAD_MPCDS_log_time_interaction.csv"), row.names = FALSE)

# 3. Split follow-up into 0-1, 1-3 and >3 year intervals
split_dat <- survival::survSplit(
  Surv(OS_time, OS_status) ~ .,
  data = dat, cut = c(365.25, 1095.75),
  start = "tstart", episode = "period"
)
stopifnot(sum(split_dat$OS_status) == 182L,
          all(split_dat$OS_time > split_dat$tstart),
          setequal(unique(split_dat$period), 1:3))

split_dat$MPCDS_0_1 <- ifelse(split_dat$period == 1L, split_dat$MPCDS_z, 0)
split_dat$MPCDS_1_3 <- ifelse(split_dat$period == 2L, split_dat$MPCDS_z, 0)
split_dat$MPCDS_after_3 <- ifelse(split_dat$period == 3L, split_dat$MPCDS_z, 0)

fit_period <- survival::coxph(
  survival::Surv(tstart, OS_time, OS_status) ~
    MPCDS_0_1 + MPCDS_1_3 + MPCDS_after_3 +
    proliferation_z + relative_purity_proxy_z,
  data = split_dat, ties = "efron"
)
fit_constant <- survival::coxph(
  survival::Surv(tstart, OS_time, OS_status) ~
    MPCDS_z + proliferation_z + relative_purity_proxy_z,
  data = split_dat, ties = "efron"
)

terms <- c("MPCDS_0_1", "MPCDS_1_3", "MPCDS_after_3")
ci <- stats::confint(fit_period)
period_results <- data.frame(
  period = c("0-1 year", "1-3 years", "Beyond 3 years"),
  events = as.integer(table(factor(split_dat$period[split_dat$OS_status == 1L],
                                   levels = 1:3))),
  HR = as.numeric(exp(coef(fit_period)[terms])),
  lower_CI = as.numeric(exp(ci[terms, 1])),
  upper_CI = as.numeric(exp(ci[terms, 2])),
  p_value = as.numeric(summary(fit_period)$coefficients[terms, "Pr(>|z|)"])
)
write.csv(period_results, file.path(output_dir, "LUAD_MPCDS_time_period_HRs.csv"), row.names = FALSE)

ll1 <- stats::logLik(fit_constant)
ll2 <- stats::logLik(fit_period)
lr_stat <- 2 * (as.numeric(ll2) - as.numeric(ll1))
lr_df <- attr(ll2, "df") - attr(ll1, "df")
comparison <- data.frame(
  chi_square = lr_stat, df = lr_df,
  p_value = stats::pchisq(lr_stat, df = lr_df, lower.tail = FALSE)
)
write.csv(comparison, file.path(output_dir, "LUAD_MPCDS_time_model_comparison.csv"), row.names = FALSE)

cat("\n=== PH DIAGNOSTICS ===\n")
print(ph)
cat("\n=== LOG-TIME MPCDS INTERACTION ===\n")
print(tv_coef["tt(MPCDS_z)", , drop = FALSE])
cat("\n=== EVENTS BY PERIOD ===\n")
print(table(split_dat$period, split_dat$OS_status))
cat("\n=== TIME-PERIOD MPCDS HAZARD RATIOS ===\n")
print(period_results, row.names = FALSE)
cat("\n=== MODEL COMPARISON ===\n")
print(comparison, row.names = FALSE)
cat("\nAnalysis results written to:\n", output_dir, "\n")
cat("Coloured figures written to:\n", figure_dir, "\n")
