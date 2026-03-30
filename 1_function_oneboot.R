run_one_bootstrap <- function(X, W, Y, train_idx, valid_idx) {
  
  X_train <- X[train_idx, , drop = FALSE]
  W_train <- W[train_idx]
  Y_train <- Y[train_idx]
  
  X_valid <- X[valid_idx, , drop = FALSE]
  W_valid <- W[valid_idx]
  Y_valid <- Y[valid_idx]
  
  preds_valid <- list()
  preds_train <- list()
  imps        <- list()     # 变量重要性保存     
  
  # ==================== 1. Causal Forest ====================
  cat("Causal Forest", "\n")
  model_cf <- tryCatch({
    causal_forest(X = X_train,
     Y = Y_train, 
     W = W_train,
     num.trees = 2000, 
     min.node.size = 10,
     honesty = TRUE, 
     seed = 2026)
  }, error = function(e) { message("CF failed"); NULL })
  
 if (!is.null(model_cf)) {
    preds_valid$cf <- predict(model_cf, X_valid)$predictions
    preds_train$cf <- predict(model_cf, X_train)$predictions
    imps$cf <- variable_importance(model_cf)              
    names(imps$cf) <- colnames(X)
  }
  
  # ==================== 2. rboost ====================
  cat("rboost", "\n")
  model_rboost <- tryCatch({
    rboost(x = X_train, 
    w = W_train, 
    y = Y_train,
    k_folds = 5, 
    nthread = 1, 
    num_search_rounds = 10)
  }, error = function(e) { message("rboost failed"); NULL })
  
if (!is.null(model_rboost)) {
    preds_valid$rboost <- predict(model_rboost, newx = X_valid)
    preds_train$rboost <- predict(model_rboost, newx = X_train)
    imps$rboost <- xgboost::xgb.importance(               
      feature_names = colnames(X),
      model = model_rboost$tau_fit$xgb_fit
    )
  }
  
  # ==================== 3. uboost ====================
  cat("uboost", "\n")
  model_uboost <- tryCatch({
    uboost(x = X_train, 
    w = W_train, 
    y = Y_train,
    k_folds = 5, 
    num_search_rounds = 10, 
    nthread = 1)
  }, error = function(e) { message("uboost failed"); NULL })
  
if (!is.null(model_uboost)) {
    preds_valid$uboost <- predict(model_uboost, newx = X_valid)
    preds_train$uboost <- predict(model_uboost, newx = X_train)
    imps$uboost <- xgboost::xgb.importance(
      feature_names = colnames(X),
      model = model_uboost$tau_fit$xgb_fit
    )
  }
  
  # ==================== 4. Bayesian Causal Forest ====================
  cat("bayesion causal forest", "\n")
  model_bcf <- tryCatch({
    stochtree::bcf(X_train = X_train, 
    Z_train = W_train, 
    y_train = Y_train,
    num_gfr = 15, 
    num_burnin = 200, 
    num_mcmc = 500,
    prognostic_forest_params = list(min_samples_leaf = 10),
    treatment_effect_forest_params = list(min_samples_leaf = 10))
  }, error = function(e) { message("BCF failed"); NULL })
  
if (!is.null(model_bcf)) {
    preds_valid$bcf <- rowMeans(predict(model_bcf, X = X_valid, Z = W_valid, terms = "tau"))
    preds_train$bcf <- rowMeans(predict(model_bcf, X = X_train, Z = W_train, terms = "tau"))
    p <- ncol(X)
    tau_splits <- model_bcf$forests_tau$get_aggregate_split_counts(p)
    imps$bcf <- tau_splits / sum(tau_splits)               
    names(imps$bcf) <- colnames(X)
  }
  
  # ==================== 5. xlasso ====================
  cat("xlasso", "\n")
  model_xlasso <- tryCatch({
    xlasso(x = X_train, 
    w = W_train, 
    y = Y_train,
    alpha = 0.5,
    k_folds_p = 5, 
    k_folds_mu1 = 5, 
    k_folds_mu0 = 5,
    lambda_choice = "lambda.min")
  }, error = function(e) { message("xlasso failed"); NULL })
  
if (!is.null(model_xlasso)) {
    preds_valid$xlasso <- predict(model_xlasso, newx = X_valid)
    preds_train$xlasso <- predict(model_xlasso, newx = X_train)
    
    coefs <- coef(model_xlasso$t_1_fit, s = "lambda.min") # 以t1为例
    coef_df <- as.data.frame(as.matrix(coefs))
    colnames(coef_df) <- "Coefficient"
    res_df <- data.frame(Variable = c("Intercept", colnames(X)),
                         Coefficient = coef_df$Coefficient)
    important_vars <- res_df[res_df$Variable != "Intercept" & res_df$Coefficient != 0, ]
    important_vars$AbsCoef <- abs(important_vars$Coefficient)
    imps$xlasso <- important_vars                         # ← 原始 coef + AbsCoef
  }
  
  # ====================== 返回结果 ======================
  valid_df <- data.table(
    iter = NA_integer_,                    # 外层循环再填
    id = valid_idx,
    observed_Y = Y_valid,
    assigned_W = W_valid,
    pred_ite_cf     = preds_valid$cf,
    pred_ite_rboost = preds_valid$rboost,
    pred_ite_bcf    = preds_valid$bcf,
    pred_ite_uboost = preds_valid$uboost,
    pred_ite_xlasso = preds_valid$xlasso
  )

  train_df <- data.table(iter = NA_integer_, id = train_idx,
                         observed_Y = Y_train, assigned_W = W_train,
                         pred_ite_cf = preds_train$cf,
                         pred_ite_rboost = preds_train$rboost,
                         pred_ite_bcf = preds_train$bcf,
                         pred_ite_uboost = preds_train$uboost,
                         pred_ite_xlasso = preds_train$xlasso)

  list(valid_df = valid_df, train_df = train_df, imps = imps)
}
