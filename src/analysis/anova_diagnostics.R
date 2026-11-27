# Run from the project root after model_results.csv has been generated.
# Rscript --vanilla src/analysis/anova_diagnostics.R
library(nlme)

input <- "results/tables/model_results.csv"
tables_dir <- "results/tables"
figures_dir <- "results/figures"
d <- read.csv(input, stringsAsFactors = FALSE)

required <- c("dataset_id", "model", "relationship", "noise", "test_rmse")
if (!all(required %in% names(d))) stop("The results table is missing columns.")
if (anyNA(d[required])) stop("The analysis columns contain missing values.")
if (!is.numeric(d$test_rmse) || any(!is.finite(d$test_rmse)) ||
    any(d$test_rmse <= 0)) stop("Test RMSE must be finite and positive.")

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
  test_rmse ~ relationship * noise * model + Error(dataset_id / model),
  data = d
)
initial_fit <- lm(
  test_rmse ~ dataset_id + model + model:relationship +
    model:noise + model:relationship:noise,
  data = d
)

# Treatment-specific variances and a shared unstructured correlation matrix.
d$model_order <- as.integer(d$model)
d$treatment <- interaction(d$relationship, d$noise, d$model)
contrasts(d$relationship) <- contr.sum(2)
contrasts(d$noise) <- contr.sum(2)
contrasts(d$model) <- contr.sum(3)
gls_fit <- nlme::gls(
  test_rmse ~ relationship * noise * model,
  data = d,
  correlation = nlme::corSymm(form = ~ model_order | dataset_id),
  weights = nlme::varIdent(form = ~ 1 | treatment),
  method = "REML",
  control = nlme::glsControl(msMaxIter = 300)
)

means <- aggregate(test_rmse ~ relationship + noise + model, d, mean)
deviations <- aggregate(test_rmse ~ relationship + noise + model, d, sd)
names(means)[4] <- "mean_test_rmse"
names(deviations)[4] <- "sd_test_rmse"
descriptive <- merge(means, deviations, by = c("relationship", "noise", "model"))
descriptive$n <- 10L
descriptive <- descriptive[order(descriptive$relationship, descriptive$noise,
                               descriptive$model), ]
write.csv(descriptive, file.path(tables_dir, "rmse_summary.csv"), row.names = FALSE)

initial_summary <- summary(blocked_anova)
anova_rows <- lapply(names(initial_summary), function(stratum) {
  tab <- initial_summary[[stratum]][[1]]
  data.frame(
    stratum = sub("Error: ", "", stratum, fixed = TRUE),
    term = trimws(rownames(tab)), df = tab$Df,
    sum_sq = tab$`Sum Sq`, mean_sq = tab$`Mean Sq`,
    f_value = tab$`F value`, p_value = tab$`Pr(>F)`,
    row.names = NULL
  )
})
write.csv(do.call(rbind, anova_rows),
          file.path(tables_dir, "blocked_anova.csv"), row.names = FALSE)

gls_tests <- anova(gls_fit, type = "marginal")
gls_table <- data.frame(
  term = rownames(gls_tests), numerator_df = gls_tests$numDF,
  denominator_df = gls_fit$dims$N - gls_fit$dims$p,
  f_value = gls_tests$`F-value`, p_value = gls_tests$`p-value`,
  row.names = NULL
)
write.csv(gls_table, file.path(tables_dir, "gls_anova.csv"), row.names = FALSE)

normalized <- as.numeric(residuals(gls_fit, type = "normalized"))
png(file.path(figures_dir, "anova_diagnostics.png"),
    width = 2100, height = 1550, res = 180)
tryCatch({
  par(mfrow = c(2, 2), mar = c(4.4, 4.5, 3.2, 1),
      oma = c(2.5, 0, 2, 0), las = 1, cex = 0.95)
  plot(fitted(initial_fit), residuals(initial_fit), pch = 16,
       col = "gray35",
       xlab = "Fitted test RMSE", ylab = "Residual",
       main = "A. Initial blocked model")
  abline(h = 0, lty = 2, col = "gray50")
  qqnorm(residuals(initial_fit), pch = 16, col = "gray35",
         main = "B. Initial residual Q-Q plot")
  qqline(residuals(initial_fit), col = "black", lwd = 2, lty = 2)
  plot(fitted(gls_fit), normalized, pch = 16,
       col = "gray35",
       xlab = "Fitted test RMSE", ylab = "Normalized residual",
       main = "C. Adjusted GLS model")
  abline(h = 0, lty = 2, col = "gray50")
  qqnorm(normalized, pch = 16, col = "gray35",
         main = "D. GLS residual Q-Q plot")
  qqline(normalized, col = "black", lwd = 2, lty = 2)
  mtext("Residual diagnostics", outer = TRUE, side = 3, line = .4,
        font = 2, cex = 1.25)
  mtext("GLS allows treatment-specific variances and correlation within datasets.",
        outer = TRUE, side = 1, line = .2, cex = .85)
}, finally = dev.off())

png(file.path(figures_dir, "model_performance.png"),
    width = 2100, height = 1450, res = 180)
tryCatch({
  par(mfrow = c(2, 2), mar = c(4, 4.8, 3, 1),
      oma = c(2.8, 0, 2.4, 0), las = 1)
  for (noise_level in levels(d$noise)) {
    for (relation in levels(d$relationship)) {
      part <- d[d$relationship == relation & d$noise == noise_level, ]
      m <- tapply(part$test_rmse, part$model, mean)
      plot(1:3, m, type = "n", log = "y", ylim = c(.35, 45),
           xlim = c(.8, 3.2), xaxt = "n", yaxt = "n", xlab = "",
           ylab = "Test RMSE (log display scale)",
           main = paste(tools::toTitleCase(relation), "-", noise_level, "noise"))
      axis(1, at = 1:3, labels = c("OLS", "kNN", "Random forest"))
      axis(2, at = c(.5, 1, 2, 5, 10, 30), labels = c("0.5", "1", "2", "5", "10", "30"))
      abline(h = c(.5, 1, 2, 5, 10, 30), col = "gray90")
      for (id in unique(part$dataset_id)) {
        pair <- part[part$dataset_id == id, ]
        lines(as.integer(pair$model), pair$test_rmse,
              type = "o", pch = 16, cex = .5, col = "gray75")
      }
      lines(1:3, m, type = "o", pch = 16, lwd = 2.5, col = "black")
      text(1:3, m, labels = sprintf("%.2f", m), pos = 3,
           offset = .6, col = "black", cex = .95)
    }
  }
  mtext("Model performance across data conditions", outer = TRUE,
        side = 3, line = .6, font = 2, cex = 1.25)
  mtext("Black: arithmetic mean. Gray: paired results from each of ten datasets.",
        outer = TRUE, side = 1, line = .2, cex = .85)
  mtext("Lower is better. The logarithmic display does not transform the analysis response.",
        outer = TRUE, side = 1, line = 1.4, cex = .8)
}, finally = dev.off())

print(descriptive, row.names = FALSE)
print(initial_summary)
print(gls_tests)
