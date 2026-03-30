# ====================== compute_qini 函数 ======================

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



# ====================== c-for-benefit 函数 ======================
library(data.table)

c_for_benefit <- function(ITE, Y, W) {
  set.seed(202603)
  # ITE 是预测的个体治疗效应（负向结局：越负越好）
  # 所以我们先取负号，让越高越好
  benefit_score <- -ITE
  dt <- data.table(benefit_score = benefit_score, Y = Y, W = W)
  
  # 分离 treated 和 untreated
  treated   <- dt[W == 1]
  untreated <- dt[W == 0]
  # 按 benefit_score 从低到高排序
  setorder(treated,   benefit_score)
  setorder(untreated, benefit_score)
  
  # 随机子采样使两组人数相等
  n <- min(nrow(treated), nrow(untreated))
  if (n == 0) return(NA_real_)
  
  treated   <- treated[sample(.N, n)]
  untreated <- untreated[sample(.N, n)]
  
  # 配对（按排序配对）
  pairs <- data.table(
    benefit_t = treated$benefit_score,
    benefit_u = untreated$benefit_score,
    Y_t = treated$Y,
    Y_u = untreated$Y
  )
  
  # observed benefit
  pairs[, obs_ben := fcase(
    Y_t == 0 & Y_u == 0, 0,
    Y_t == 0 & Y_u == 1, -1,   # treated 好于 untreated
    Y_t == 1 & Y_u == 0,  1,   # treated 差于 untreated
    Y_t == 1 & Y_u == 1,  0
  )]
  
  # predicted benefit（取平均）
  pairs[, pred_ben := (benefit_t + benefit_u) / 2]
  
  # 计算 concordance（所有 obs_ben 不同的配对中，预测排序一致的比例）
  diff_idx <- which(pairs$obs_ben != 0)
  if (length(diff_idx) == 0) return(NA_real_)
  
  concordant <- sum(
    (pairs$obs_ben[diff_idx] ==  1 & pairs$pred_ben[diff_idx] > 0) |
      (pairs$obs_ben[diff_idx] == -1 & pairs$pred_ben[diff_idx] < 0)
  )
  
  return(concordant / length(diff_idx))
}
