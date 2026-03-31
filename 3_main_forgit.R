library(progress)
library(data.table)
library(dplyr)
library(tidymodels)
library(rlearner)
library(grf)
library(stochtree)
library(xgboost)

n_obs      <- nrow(dt_model)
start_iter <- 1
end_iter   <- 500        

set.seed(202603)
all_valid_results <- list()
all_train_results <- list()
importance_list   <- list()

pb <- progress_bar$new(format = "  Bootstrap [:bar] :percent  ETA: :eta",
                       total = end_iter - start_iter + 1, clear = FALSE)

for (i in start_iter:end_iter) {
  print(paste("Iteration", i))
  
  train_idx <- sample(1:n_obs, n_obs, replace = TRUE)
  valid_idx <- setdiff(1:n_obs, unique(train_idx))
  
  result <- run_one_bootstrap(X, W, Y, train_idx, valid_idx)
  
  result$valid_df$iter <- i
  result$train_df$iter <- i
  
  all_valid_results[[i - start_iter + 1]] <- result$valid_df
  all_train_results[[i - start_iter + 1]] <- result$train_df
  if (length(result$imps) > 0) {
    importance_list[[i]] <- list(iter = i, imps = result$imps)
  }
  
  rm(result)
  gc()
  pb$tick()
}

# Fit the best model on the full dataset
set.seed(202603)
model_uboost <- uboost(x = X, 
w = W,  
y = Y,
k_folds = 5, 
num_search_rounds = 100,
nthread = 1)