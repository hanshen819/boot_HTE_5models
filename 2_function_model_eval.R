# ====================== compute_qini  ======================
# Ref: Individualized Treatment Effects of Oxygen Targets  in Mechanically Ventilated Critically Ill Adults
library(tools4uplift)
compute_qini <- function(df, model_col) {
  df <- as.data.table(df)
  pred_col <- grep(paste0("pred_ite_", model_col), names(df), value = TRUE)[1]

  df_trans <- df %>%
    mutate(survival = ifelse(observed_Y == 1, 0, 1)) %>%
    mutate(!!sym(paste0("opp_", model_col)) := -df[[pred_col]])
 
  perf <- PerformanceUplift(data = df_trans,
                            treat = "assigned_W",
                            outcome = "observed_Y",
                            prediction = paste0("opp_", model_col),
                            nb.group = 3, 
                            rank.precision = 2)

  data.table(
    QINI     = QiniArea(perf, adjusted = FALSE),
    Adj_QINI = QiniArea(perf, adjusted = TRUE)
  )
}

# ====================== c-statistic  ======================
# Ref: Individualized Treatment Effects of Oxygen Targets  in Mechanically Ventilated Critically Ill Adults
c_statistic <- function(pred_rr, y, w) {
  set.seed(202603)
  stopifnot(length(pred_rr) == length(y), length(y) == length(w))
  
  dt <- data.table(pred_rr = pred_rr, y = y, w = w)
  untreated <- dt[w == 0]
  treated   <- dt[w == 1]
  
  n_untreated <- nrow(untreated)
  n_treated   <- nrow(treated)
  
  if (n_treated < n_untreated) {
    untreated <- untreated[sample(.N, n_treated)]
  } else if (n_untreated < n_treated) {
    treated <- treated[sample(.N, n_treated)]
  }
  
  stopifnot(nrow(untreated) == nrow(treated))
  
  setorder(untreated, pred_rr)
  setorder(treated, pred_rr)
  
  pairs_dt <- data.table(
    obs_ben = treated$y - untreated$y,
    pred_ben = (treated$pred_rr + untreated$pred_rr) / 2
  )
  
  obs_benefit <- pairs_dt$obs_ben
  pred_benefit <- pairs_dt$pred_ben
  n_pairs <- nrow(pairs_dt)
  
  count <- 0
  total <- 0
  
  # 如果配对数少于2，无法进行组合比较
  if (n_pairs < 2) return(NA_real_)
  for (i in 1:(n_pairs - 1)) {
    for (j in (i + 1):n_pairs) {
      # 只有当观测收益不等时才进入比较
      if (obs_benefit[i] != obs_benefit[j]) {
        if ((obs_benefit[i] < obs_benefit[j] && pred_benefit[i] < pred_benefit[j]) ||
            (obs_benefit[i] > obs_benefit[j] && pred_benefit[i] > pred_benefit[j])) {
          count <- count + 1
        }
        total <- total + 1
      }
    }
  }
  # 返回最终的一致性比例
  if (total == 0) return(NA_real_)
  return(count / total)
}

# ====================== q-statistic  ======================
# Ref: Estimating individual treatment effects on COPD  exacerbations by causal machine learning on  randomised controlled trials
q_score <- function(t, y, w, 
                    num_quantiles = c(2, 4, 8, 16, 32, 64), 
                    weights = NULL) {
  set.seed(202603)
  
  N <- length(t)
  stopifnot(length(y) == N, length(w) == N, all(w %in% c(0, 1)))
  
  ord <- order(t)
  t_sorted <- t[ord]
  y_sorted <- y[ord]
  w_sorted <- w[ord]
  
  treated <- w_sorted == 1
  global_ate <- mean(y_sorted[treated]) - mean(y_sorted[!treated])
  
  if (is.null(weights)) weights <- rep(1, length(num_quantiles))
  stopifnot(length(weights) == length(num_quantiles))
  
  L_Q_total   <- 0
  L_ATE_total <- 0
  
  for (m in seq_along(num_quantiles)) {
    k <- num_quantiles[m]
    if (k >= N) break
    
    group <- cut(1:N, breaks = k, labels = FALSE, include.lowest = TRUE)
    
    lq_vec   <- numeric(0)
    late_vec <- numeric(0)
    
    for (j in 1:k) {
      idx <- which(group == j)
      if (length(idx) < 2) next
      
      sub_treated <- y_sorted[idx][w_sorted[idx] == 1]
      sub_control <- y_sorted[idx][w_sorted[idx] == 0]
      
      if (length(sub_treated) == 0 || length(sub_control) == 0) next
      
      obs_sub_ate <- mean(sub_treated) - mean(sub_control)      
      pred_ate    <- mean(t_sorted[idx])                        
      
      lq_vec   <- c(lq_vec,   obs_sub_ate - pred_ate)
      late_vec <- c(late_vec, global_ate  - pred_ate)
    }
    
    if (length(lq_vec) == 0) next
    
    lq_k   <- sqrt(mean(lq_vec^2))
    late_k <- sqrt(mean(late_vec^2))
    
    L_Q_total   <- L_Q_total   + weights[m] * lq_k
    L_ATE_total <- L_ATE_total + weights[m] * late_k
  }
  
  if (L_ATE_total == 0) {
    warning("L_ATE_total equals 0, cannot calculate Q-score")
    return(NA_real_)
  }
  
  Q <- 1 - L_Q_total / L_ATE_total
  return(Q)
}

compute_qscore_for_model <- function(df, model_col) {
  df <- as.data.table(df)
  pred_col <- grep(paste0("pred_ite_", model_col), names(df), value = TRUE)[1]
  
  t_val <- df[[pred_col]]
  y_val <- df$observed_Y
  w_val <- df$assigned_W
  
  qs <- q_score(t = -t_val, y = y_val, w = w_val)
  
  return(data.table(Q_Score = qs))
}

