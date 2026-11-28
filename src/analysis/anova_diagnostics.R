# Run from the project root after model_results.csv has been generated.
# Rscript --vanilla src/analysis/anova_diagnostics.R

input <- "results/tables/model_results.csv"
tables_dir <- "results/tables"
figures_dir <- "results/figures"
d <- read.csv(input, stringsAsFactors = FALSE)

required <- c("dataset_id", "model", "relationship", "noise", "noise_free_rmse", "test_rmse")
if (!all(required %in% names(d))) stop("The results table is missing columns.")
if (anyNA(d[required])) stop("The analysis columns contain missing values.")
for (response in c("noise_free_rmse", "test_rmse")) {
  if (!is.numeric(d[[response]]) || any(!is.finite(d[[response]])) ||
      any(d[[response]] <= 0)) stop(paste(response, "must be finite and positive."))
}

model_levels <- c("OLS", "kNN", "random forest")
if (!setequal(d$model, model_levels) ||
    !setequal(d$relationship, c("linear", "quadratic")) ||
    !setequal(d$noise, c("low", "high"))) {
  stop("Unexpected factor levels in the results table.")
}
d$dataset_id <- factor(d$dataset_id)
d$model <- factor(d$model, levels = model_levels)
d$relationship <- factor(d$relationship, levels = c("linear", "quadratic"))
d$noise <- factor(d$noise, levels = c("low", "high"))

if (nrow(d) != 120L || nlevels(d$dataset_id) != 40L ||
    any(table(d$dataset_id, d$model) != 1L)) {
  stop("Expected 40 datasets with one evaluation of each of the three models.")
}
scenario_counts <- tapply(
  paste(d$relationship, d$noise), d$dataset_id,
  function(x) length(unique(x))
)
if (any(scenario_counts != 1L)) stop("A dataset has inconsistent scenario labels.")
datasets <- unique(d[c("dataset_id", "relationship", "noise")])
if (any(table(datasets$relationship, datasets$noise) != 10L)) {
  stop("Expected ten independent datasets per relationship-noise combination.")
}

d <- d[order(d$dataset_id, d$model), ]
rownames(d) <- NULL
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Relationship and noise vary between datasets. Model varies within datasets.
blocked_anova <- aov(
  log(noise_free_rmse) ~ model * relationship * noise + Error(dataset_id),
  data = d
)
d$log_noise_free_rmse <- log(d$noise_free_rmse)

# Check each error level separately. These lm fits are for plots;
# all significance tests come from the aov model above.
dataset_means <- aggregate(log_noise_free_rmse ~ dataset_id + relationship + noise,
                           d, mean)
between_fit <- lm(log_noise_free_rmse ~ relationship * noise, data = dataset_means)
within_fit <- lm(
  log_noise_free_rmse ~ dataset_id + model + model:relationship +
    model:noise + model:relationship:noise,
  data = d
)

means <- aggregate(cbind(noise_free_rmse, test_rmse, log_noise_free_rmse) ~
                     relationship + noise + model, d, mean)
deviations <- aggregate(cbind(noise_free_rmse, test_rmse) ~
                          relationship + noise + model, d, sd)
names(means)[4:6] <- c("mean_noise_free_rmse", "mean_test_rmse", "mean_log_noise_free_rmse")
names(deviations)[4:5] <- c("sd_noise_free_rmse", "sd_test_rmse")
descriptive <- merge(means, deviations, by = c("relationship", "noise", "model"))
descriptive$geometric_mean_noise_free_rmse <- exp(descriptive$mean_log_noise_free_rmse)
descriptive$n <- 10L
descriptive <- descriptive[order(descriptive$relationship, descriptive$noise,
                               descriptive$model), ]
write.csv(descriptive, file.path(tables_dir, "rmse_summary.csv"), row.names = FALSE)

initial_summary <- summary(blocked_anova)
anova_rows <- lapply(names(initial_summary), function(stratum) {
  tab <- initial_summary[[stratum]][[1]]
  data.frame(
    stratum = sub("Error: ", "", stratum, fixed = TRUE),
    response = "log(noise_free_rmse)",
    term = trimws(rownames(tab)), df = tab$Df,
    sum_sq = tab$`Sum Sq`, mean_sq = tab$`Mean Sq`,
    f_value = tab$`F value`, p_value = tab$`Pr(>F)`,
    row.names = NULL
  )
})
write.csv(do.call(rbind, anova_rows),
          file.path(tables_dir, "blocked_anova.csv"), row.names = FALSE)

png(file.path(figures_dir, "anova_diagnostics.png"),
    width = 2100, height = 1550, res = 180)
tryCatch({
  par(mfrow = c(2, 2), mar = c(4.4, 4.5, 3.2, 1),
      oma = c(2.5, 0, 2, 0), las = 1, cex = 0.95)
  plot(fitted(between_fit), residuals(between_fit), pch = 16,
       col = "gray35",
       xlab = "Fitted dataset mean of log(RMSE)", ylab = "Residual",
       main = "A. Between-dataset residuals")
  abline(h = 0, lty = 2, col = "gray50")
  qqnorm(residuals(between_fit), pch = 16, col = "gray35",
         main = "B. Between-dataset Q-Q plot")
  qqline(residuals(between_fit), col = "black", lwd = 2, lty = 2)
  plot(fitted(within_fit), residuals(within_fit), pch = 16,
       col = "gray35",
       xlab = "Fitted log(noise-free RMSE)", ylab = "Residual",
       main = "C. Within-dataset residuals")
  abline(h = 0, lty = 2, col = "gray50")
  qqnorm(residuals(within_fit), pch = 16, col = "gray35",
         main = "D. Within-dataset Q-Q plot")
  qqline(residuals(within_fit), col = "black", lwd = 2, lty = 2)
  mtext("ANOVA diagnostics: log(noise-free RMSE)", outer = TRUE, side = 3, line = .4,
        font = 2, cex = 1.25)
  mtext("Dataset means: 40 datasets. Within-dataset residuals: 120 paired model evaluations.",
        outer = TRUE, side = 1, line = .2, cex = .85)
}, finally = dev.off())

png(file.path(figures_dir, "model_performance.png"),
    width = 2100, height = 1450, res = 180)
tryCatch({
  par(mfrow = c(2, 2), mar = c(4, 5.5, 3, 1),
      oma = c(2.8, 0, 2.4, 0), las = 1)
  for (noise_level in levels(d$noise)) {
    for (relation in levels(d$relationship)) {
      part <- d[d$relationship == relation & d$noise == noise_level, ]
      m <- exp(tapply(part$log_noise_free_rmse, part$model, mean))
      plot(1:3, m, type = "n", log = "y", ylim = range(d$noise_free_rmse) * c(.75, 1.4),
           xlim = c(.8, 3.2), xaxt = "n", yaxt = "n", xlab = "",
           ylab = "Noise-free RMSE (log scale)",
           main = paste(tools::toTitleCase(relation), "-", noise_level, "noise"))
      axis(1, at = 1:3, labels = c("OLS", "kNN", "Random forest"))
      ticks <- c(.005, .02, .1, .5, 2, 10, 30)
      axis(2, at = ticks, labels = c("0.005", "0.02", "0.1", "0.5", "2", "10", "30"))
      abline(h = ticks, col = "gray90")
      for (id in unique(part$dataset_id)) {
        pair <- part[part$dataset_id == id, ]
        lines(as.integer(pair$model), pair$noise_free_rmse,
              type = "o", pch = 16, cex = .5, col = "gray75")
      }
      lines(1:3, m, type = "o", pch = 16, lwd = 2.5, col = "black")
      text(1:3, m, labels = sprintf("%.3f", m), pos = 3,
           offset = .6, col = "black", cex = .95)
    }
  }
  mtext("Model performance across data conditions", outer = TRUE,
        side = 3, line = .6, font = 2, cex = 1.25)
  mtext("Black: geometric mean. Gray: paired results from each of ten datasets.",
        outer = TRUE, side = 1, line = .2, cex = .85)
  mtext("Lower is better. ANOVA uses log(noise-free RMSE); test RMSE is retained in the summary table.",
        outer = TRUE, side = 1, line = 1.4, cex = .8)
}, finally = dev.off())

print(descriptive, row.names = FALSE)
print(initial_summary)
