# 参照Throx计算每组下干预组和对照组的RD值
calc_stats <- function(df, n_boot = 1000) {
  # 定义计算分组发生率的函数（用于bootstrap）
  set.seed(202603)
  calc_rate <- function(data, group_val) {
    group_data <- data[data$Group == group_val, ]
    n <- nrow(group_data)
    if (n == 0) return(0) # 避免空分组报错
    evt <- sum(group_data$Y)
    return(evt / n)
  }
  # 基础统计量计算
  # 治疗组
  n_trt <- sum(df$Group == "1") 
  evt_trt <- sum(df$Y[df$Group == "1"])
  rate_trt <- evt_trt / n_trt
  # 对照组
  n_ctrl <- sum(df$Group != "1")
  evt_ctrl <- sum(df$Y[df$Group != "1"])
  rate_ctrl <- evt_ctrl / n_ctrl
  # 分位数法构建95%置信区间
  # 方法：通过bootstrap重抽样获取发生率的分布，再取2.5%和97.5%分位数
  # 1. 治疗组bootstrap
  set.seed(202603) #
  boot_trt <- replicate(n_boot, {
    boot_df <- df[sample(nrow(df), replace = TRUE), ]
    calc_rate(boot_df, "1")
  })
  
  # 2. 对照组bootstrap
  boot_ctrl <- replicate(n_boot, {
    boot_df <- df[sample(nrow(df), replace = TRUE), ]
    calc_rate(boot_df, "0") # 对照组Group标记为0（与原逻辑一致）
  })
  # 提取分位数置信区间（2.5%和97.5%）
  trt_lower <- quantile(boot_trt, 0.025)
  trt_upper <- quantile(boot_trt, 0.975)
  ctrl_lower <- quantile(boot_ctrl, 0.025)
  ctrl_upper <- quantile(boot_ctrl, 0.975)
  
  # 输出结果
  data.frame(
    Group_Type = c("1", "0"),
    Rate = c(rate_trt, rate_ctrl),
    Lower = c(trt_lower, ctrl_lower),
    Upper = c(trt_upper, ctrl_upper), 
    N = c(n_trt, n_ctrl)
  )
}


# 参照JAMA计算每个获益组的RD值
calc_boot_effect <- function(sub_df, target_col, n_boot = 1000, seed = 123) {
  # 设置随机数种子以保证结果可复现
  set.seed(202603)
  # --- 内部辅助函数：计算单次 RD 和 RR ---
  compute_stats <- function(df) {
    # 联合干预组 (Group == "1") vs 对照组
    n_trt <- sum(df$Group == "1")
    n_ctrl <- sum(df$Group != "1")
    # 防止某一分组在重抽样中完全消失导致的报错
    if (n_trt == 0 || n_ctrl == 0) return(c(NA, NA))
    
    ev_trt <- sum(df[[target_col]] == 1 & df$Group == "1")
    ev_ctrl <- sum(df[[target_col]] == 1 & df$Group != "1")
    
    rate_trt <- ev_trt / n_trt
    rate_ctrl <- ev_ctrl / n_ctrl
    
    # 计算 RD
    rd <- rate_trt - rate_ctrl  
    # 计算 RR (防止分母为0，如果对照组发生率为0，RR 设为 NA 或无穷大)
    if (rate_ctrl == 0) {
      rr <- NA 
    } else {
      rr <- rate_trt / rate_ctrl
    }
    return(c(rd, rr))
  }
  
  # 1. 计算原始数据的点估计 (Point Estimates)
  orig_stats <- compute_stats(sub_df)
  est_rd <- orig_stats[1]
  est_rr <- orig_stats[2]
  
  # 2. 计算 P值 (保持原有的卡方检验或Fisher精确检验)
  # 注意：Bootstrap 也可以算 P 值，但保留卡方检验通常更方便对比
  n_trt_orig <- sum(sub_df$Group == "1")
  n_ctrl_orig <- sum(sub_df$Group != "1")
  ev_trt_orig <- sum(sub_df[[target_col]] == 1 & sub_df$Group == "1")
  ev_ctrl_orig <- sum(sub_df[[target_col]] == 1 & sub_df$Group != "1")
  
  # 使用 tryCatch 防止样本量过小导致的报错
  p_val <- tryCatch({
    prop.test(c(ev_trt_orig, ev_ctrl_orig), c(n_trt_orig, n_ctrl_orig))$p.value
  }, error = function(e) NA)
  
  # 3. 执行 Bootstrap 重抽样
  # replicate 会运行 n_boot 次，每次重新采样并计算
  boot_results <- replicate(n_boot, {
    # 有放回抽样
    idx <- sample(nrow(sub_df), replace = TRUE)
    boot_df <- sub_df[idx, ]
    compute_stats(boot_df)
  })
  
  # boot_results 是一个 2 x n_boot 的矩阵
  # 第一行是 RD，第二行是 RR
  boot_rd <- boot_results[1, ]
  boot_rr <- boot_results[2, ]
  
  # 4. 基于分位数计算 95% CI (排除 NA 值)
  rd_ci <- quantile(boot_rd, probs = c(0.025, 0.975), na.rm = TRUE)
  rr_ci <- quantile(boot_rr, probs = c(0.025, 0.975), na.rm = TRUE)
  
  return(data.frame(
    N = nrow(sub_df),
    Rate_Trt = ev_trt_orig / n_trt_orig,
    Rate_Ctrl = ev_ctrl_orig / n_ctrl_orig,
    # 点估计
    RD = est_rd, 
    # Bootstrap CI
    RD_Low = rd_ci[1], 
    RD_High = rd_ci[2],
    # 点估计
    RR = est_rr, 
    # Bootstrap CI
    RR_Low = rr_ci[1], 
    RR_High = rr_ci[2],
    P_Value = p_val
  ))
}




