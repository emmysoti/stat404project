source("src/data_generation/data_generation.R")
library(FNN)
library(ranger)
set.seed(979)

linear_data <- c(linear_high, linear_low)
quadratic_data <- c(quadratic_high, quadratic_low)
linear_noise <- c(rep("high", length(linear_high)), rep("low", length(linear_low)))
quadratic_noise <- c(rep("high", length(quadratic_high)), rep("low", length(quadratic_low)))
linear_splits <- list()
quadratic_splits <- list()

for (i in seq_along(linear_data)) {
  dataset <- linear_data[[i]]

  train_rows <- sample(seq_len(nrow(dataset)), size = 0.8 * nrow(dataset))
  test_rows <- setdiff(seq_len(nrow(dataset)), train_rows)

  linear_splits[[i]] <- list(train = train_rows, test = test_rows)
}

for (i in seq_along(quadratic_data)) {
  dataset <- quadratic_data[[i]]
  train_rows <- sample(seq_len(nrow(dataset)), size = 0.8 * nrow(dataset))
  test_rows <- setdiff(seq_len(nrow(dataset)), train_rows)
  quadratic_splits[[i]] <- list(train = train_rows, test = test_rows)
}

# Create shared CV folds and a random model order for each dataset.
set.seed(979)
linear_folds <- list()
quadratic_folds <- list()
linear_order <- list()
quadratic_order <- list()

for (i in seq_along(linear_data)) {
  linear_folds[[i]] <- sample(rep(1:5, length.out = length(linear_splits[[i]]$train)))
  linear_order[[i]] <- sample(c("OLS", "kNN", "random forest"))
}
for (i in seq_along(quadratic_data)) {
  quadratic_folds[[i]] <- sample(rep(1:5, length.out = length(quadratic_splits[[i]]$train)))
  quadratic_order[[i]] <- sample(c("OLS", "kNN", "random forest"))
}

k_values <- c(1, 3, 5, 10, 20, 40, 80, 160)
depths <- c(2, 4, 6, 8, 12)

fit_knn <- function(train_data, test_data, folds, true_y) {
  cv_rmse <- numeric(length(k_values))

  for (i in seq_along(k_values)) {
    fold_errors <- numeric(5)

    for (fold in 1:5) {
      cv_train <- train_data[folds != fold, ]
      cv_valid <- train_data[folds == fold, ]

      knn_cv <- knn.reg(
        train = as.matrix(cv_train["x"]),
        test = as.matrix(cv_valid["x"]),
        y = cv_train$y,
        k = k_values[i]
      )
      fold_errors[fold] <- mean((cv_valid$y - knn_cv$pred)^2)
    }

    cv_rmse[i] <- sqrt(mean(fold_errors))
  }

  best_k <- k_values[which.min(cv_rmse)]
  knn <- knn.reg(
    train = as.matrix(train_data["x"]),
    test = as.matrix(test_data["x"]),
    y = train_data$y,
    k = best_k
  )

  data.frame(
    best_k = best_k,
    noise_free_rmse = sqrt(mean((true_y - knn$pred)^2)),
    test_rmse = sqrt(mean((test_data$y - knn$pred)^2))
  )
}


fit_rf <- function(train_data, test_data, folds, true_y) {
  cv_rmse <- numeric(length(depths))

  for (i in seq_along(depths)) {
    fold_errors <- numeric(5)

    for (fold in 1:5) {
      cv_train <- train_data[folds != fold, ]
      cv_valid <- train_data[folds == fold, ]
      rf_cv <- ranger::ranger(
        y ~ x,
        data = cv_train,
        mtry = 1,
        num.trees = 500,
        max.depth = depths[i]
      )
      predictions <- predict(rf_cv, data = cv_valid)$predictions
      fold_errors[fold] <- mean((cv_valid$y - predictions)^2)
    }

    cv_rmse[i] <- sqrt(mean(fold_errors))
  }

  best_depth <- depths[which.min(cv_rmse)]
  rf <- ranger::ranger(
    y ~ x,
    data = train_data,
    mtry = 1,
    num.trees = 500,
    max.depth = best_depth
  )
  predictions <- predict(rf, data = test_data)$predictions

  data.frame(
    best_depth = best_depth,
    noise_free_rmse = sqrt(mean((true_y - predictions)^2)),
    test_rmse = sqrt(mean((test_data$y - predictions)^2))
  )
}


# Evaluate all models on each dataset in the randomized order.
results_list <- list()
result_number <- 0
set.seed(979)

for (relationship in c("linear", "quadratic")) {
  if (relationship == "linear") {
    datasets <- linear_data
    splits <- linear_splits
    all_folds <- linear_folds
    all_orders <- linear_order
    noise_labels <- linear_noise
  } else {
    datasets <- quadratic_data
    splits <- quadratic_splits
    all_folds <- quadratic_folds
    all_orders <- quadratic_order
    noise_labels <- quadratic_noise
  }

  for (i in seq_along(datasets)) {
    dataset <- datasets[[i]]
    train_data <- dataset[splits[[i]]$train, ]
    test_data <- dataset[splits[[i]]$test, ]
    if (relationship == "linear") {
      true_y <- 2 * test_data$x + 3
    } else {
      true_y <- 2 * test_data$x^2 + 3 * test_data$x + 4
    }

    for (position in seq_along(all_orders[[i]])) {
      model <- all_orders[[i]][position]
      if (model == "OLS") {
        ols <- lm(y ~ x, data = train_data)
        predictions <- predict(ols, newdata = test_data)
        model_result <- data.frame(
          noise_free_rmse = sqrt(mean((true_y - predictions)^2)),
          test_rmse = sqrt(mean((test_data$y - predictions)^2))
        )
      } else if (model == "kNN") {
        model_result <- fit_knn(train_data, test_data, all_folds[[i]], true_y)
      } else {
        model_result <- fit_rf(train_data, test_data, all_folds[[i]], true_y)
      }

      result_number <- result_number + 1
      results_list[[result_number]] <- cbind(
        dataset = i,
        dataset_id = paste(relationship, i, sep = "_"),
        model = model,
        relationship = relationship,
        noise = noise_labels[i],
        evaluation_order = position,
        model_result
      )
    }
  }
}

results <- dplyr::bind_rows(results_list)
write.csv(results, "results/tables/model_results.csv", row.names = FALSE)
