# =============================================================================
# Remotely Close Associations: Openness to Experience and Semantic Memory Structure
# 适配 SemNetCleaner 1.3.7 / R 4.5.3
# 包含交互式清洗、节点标准化、Jaccard相似度、igraph网络指标及统计检验
# =============================================================================

# -------------------------- 加载包 -----------------------------------------
library(SemNetCleaner)   # 数据清洗 (version 1.3.7)
library(psych)           # Cronbach's alpha
library(NetworkToolbox)  # 提供 jaccard, TMFG
library(igraph)          # 网络指标计算及导出
library(SemNeT)
clean_bom <- function(df) {
  colnames(df)[1] <- gsub("^\ufeff", "", colnames(df)[1])
  return(df)
}


# =============================================================================
# 第一步：交互式数据清洗（使用 textcleaner）
# =============================================================================
cat("请选择 FINAL fluency.csv 文件\n")
raw <- clean_bom(read.csv(file.choose(), as.is = TRUE))

cat("正在启动 textcleaner 交互式校正，请按提示操作...\n")
cleaned <- textcleaner(raw, miss = 99, partBY = "row")

# 提取二进制矩阵（行为被试，列为清洗后的唯一响应词）
con <- as.matrix(cleaned$responses$binary)
cat("清洗后二进制矩阵维度：", dim(con), "\n")

# =============================================================================
# 第二步：节点标准化（Node Cleaning）
# =============================================================================

cat("\n================ 节点标准化 =================\n")

old_names <- colnames(con)

new_names <- tolower(old_names)

# Christensen (2018) 数据中的已知拼写变体
new_names[new_names == "tiger"]   <- "tiger"
new_names[new_names == "Tiger"]   <- "tiger"

new_names[new_names == "racoon"]  <- "raccoon"

new_names[new_names == "cheeta"]  <- "cheetah"

# 合并重复列
dup_names <- unique(new_names[duplicated(new_names)])

if(length(dup_names) > 0){
  
  cat("发现重复节点：\n")
  print(dup_names)
  
  temp_df <- as.data.frame(con)
  
  merged_df <- sapply(
    unique(new_names),
    function(nm){
      
      cols <- which(new_names == nm)
      
      if(length(cols) == 1){
        
        temp_df[, cols]
        
      } else {
        
        as.numeric(rowSums(temp_df[, cols, drop=FALSE]) > 0)
        
      }
    }
  )
  
  con <- as.matrix(merged_df)
  
  colnames(con) <- unique(new_names)
  
} else {
  
  colnames(con) <- new_names
  
}

cat("节点标准化完成\n")
cat("最终节点数:", ncol(con), "\n")

# =============================================================================
# 第三步：手动修正特定响应（拆分 catefrog，补全缺失列与个案）
# =============================================================================

# 3.1 拆分 "catefrog" 为 "cat" 和 "frog"（如果存在）
if ("catefrog" %in% colnames(con)) {
  idx_catefrog <- which(colnames(con) == "catefrog")
  rows_with_catefrog <- which(con[, idx_catefrog] != 0)
  # 确保 "cat" 和 "frog" 列存在
  for (word in c("cat", "frog")) {
    if (!word %in% colnames(con)) {
      con <- cbind(con, 0)
      colnames(con)[ncol(con)] <- word
    }
  }
  con[rows_with_catefrog, c("cat", "frog")] <- 1
  con <- con[, -idx_catefrog]  # 删除原列
  cat("已拆分 'catefrog' 并删除原列\n")
}

# 3.2 补全缺失的动物词列（确保所有可能出现的词都在矩阵中）
animal_words <- c("cat", "dog", "fish", "elephant", "tiger", "zebra", 
                  "monkey", "giraffe", "lion", "dolphin", "chicken", "squirrel",
                  "mouse", "moose", "horse", "bear", "deer", "pig", "cow")
missing_cols <- setdiff(animal_words, colnames(con))
if (length(missing_cols) > 0) {
  new_cols <- matrix(0, nrow = nrow(con), ncol = length(missing_cols))
  colnames(new_cols) <- missing_cols
  con <- cbind(con, new_cols)
  cat("已添加缺失列：", paste(missing_cols, collapse = ", "), "\n")
}

# 3.3 补全全零个案（case 386 和 499）
zero_rows <- which(rowSums(con) == 0)
if (length(zero_rows) > 0) {
  cat("发现全零行：", zero_rows, "，将依据原始数据补全。\n")
  # 个案386
  if (386 %in% zero_rows) {
    con[386, c("cat", "dog", "fish", "elephant", "tiger", "zebra", 
               "monkey", "giraffe", "lion", "dolphin", "chicken", "squirrel")] <- 1
    cat("Case 386 已补全\n")
  }
  # 个案499
  if (499 %in% zero_rows) {
    con[499, c("dog", "cat", "mouse", "moose", "horse", "lion", 
               "tiger", "bear", "deer", "pig", "cow")] <- 1
    cat("Case 499 已补全\n")
  }
  # 再次检查
  zero_rows_after <- which(rowSums(con) == 0)
  if (length(zero_rows_after) > 0) {
    warning("补全后仍存在全零行：", zero_rows_after, "，请手动检查。")
  } else {
    cat("✅ 所有个案均有效。\n")
  }
} else {
  cat("✅ 无缺失个案。\n")
}

finalFull <- con
# 可选保存
# save(finalFull, file = "FINAL_Cleaned_Matrix.Rdata")

# =============================================================================
# 第四步：加载潜变量（Openness）并计算 Cronbach's alpha
# =============================================================================
cat("请选择 FINAL open.csv 文件\n")

latent <- read.csv(
  file.choose(),
  stringsAsFactors = FALSE
)
# =============================================================================
# 【新增】数据对齐校验与合并（按被试 ID）
# =============================================================================
# 检查两个数据框是否具有共同的 ID 列（不区分大小写）
find_id_col <- function(df) {
  ids <- grep("^id$|^subject$|^participant$", colnames(df), ignore.case = TRUE, value = TRUE)
  if (length(ids) == 1) return(ids[1]) else return(NULL)
}
id_fluency <- find_id_col(as.data.frame(raw))   # 注意 raw 是原始 fluency 数据，包含 ID
id_latent  <- find_id_col(latent)

if (!is.null(id_fluency) && !is.null(id_latent)) {
  # 两个文件都有 ID 列，按 ID 合并
  cat("按 ID 列 '", id_fluency, "' 和 '", id_latent, "' 合并数据\n", sep = "")
  # 从 finalFull 重建 ID 信息：finalFull 的行顺序与 raw 一致
  # 因此将 ID 列附加到 finalFull
  finalFull_with_id <- cbind(raw[[id_fluency]], finalFull)
  colnames(finalFull_with_id)[1] <- id_fluency
  # 合并 latent 和 finalFull
  comb <- merge(latent, finalFull_with_id, by.x = id_latent, by.y = id_fluency, all = FALSE)
  # 确保顺序不影响后续分组，按 latent 排序
  comb <- comb[order(comb$no_int), ]   # 假设 latent 中有 no_int
  # 提取 latent 列和响应矩阵
  latent_col <- comb$no_int
  finalFull_aligned <- comb[, colnames(finalFull_with_id)[colnames(finalFull_with_id) != id_fluency]]
  # 注意：如果 latent 还有其他列，我们只保留 no_int，其他忽略
  # 重新构建 comb 矩阵用于后续分析
  comb <- as.data.frame(cbind(latent_col, finalFull_aligned))
  colnames(comb)[1] <- "latent"
} else {
  # 没有找到 ID 列，退回到按行顺序拼接（但增加校验）
  cat("未找到共同的 ID 列，将按行顺序拼接。\n")
  if (nrow(finalFull) != nrow(latent)) {
    stop("❌ 行数不匹配：fluency 矩阵有 ", nrow(finalFull), " 行，latent 有 ", nrow(latent), " 行。必须一致。")
  }
  comb <- as.data.frame(cbind(latent$no_int, finalFull))
  colnames(comb)[1] <- "latent"
  warning("⚠️ 未使用 ID 列进行合并，请确保两个文件的行顺序完全对应。")
}
# ====== 在这里插入断言 ======
stopifnot(nrow(comb) == 516)
cat("✅ 合并后样本量校验通过，nrow(comb) =", nrow(comb), "\n")


# =============================================================================
# 第五步：Cronbach's alpha 及分组
# =============================================================================
# 注意：latent 原始数据可能已被合并，但我们仍用原始 latent 计算 alpha（因为未改动）
bfasInt  <- latent[, 10:19]   # bfasi01~10 → Intellect
bfasOpen <- latent[, 20:29]   # bfaso01~10 → Openness
neo      <- latent[, 30:41] 
cat("BFAS Openness alpha =", alpha(bfasOpen)$total$std.alpha, "\n")
cat("BFAS Intellect alpha =", alpha(bfasInt)$total$std.alpha, "\n")
cat("NEO Openness alpha   =", alpha(neo)$total$std.alpha, "\n")

# 现在 comb 已经包含 latent 和响应矩阵，按 latent 排序
# 分组（低分组前258，高分组后258）
low  <- comb[1:258, ]
high <- comb[259:516, ]

deLow  <- low[, -1]                         # 低分组响应矩阵
deHigh <- high[, -1]                        # 高分组响应矩阵

# =============================================================================
# 第六步：唯一响应分析（McNemar 检验）
# =============================================================================
onlyH <- deHigh[, which(colSums(deHigh) >= 1)]
onlyL <- deLow[,  which(colSums(deLow)  >= 1)]

uniH <- colnames(onlyH)
uniL <- colnames(onlyL)
cat("高分组独有响应数：", length(setdiff(uniH, uniL)), "\n")
cat("低分组独有响应数：", length(setdiff(uniL, uniH)), "\n")

uniT <- unique(c(uniH, uniL))
cat("总唯一响应数：", length(uniT), "\n")

oneH <- match(uniT, uniH)
oneL <- match(uniT, uniL)
chitest <- matrix(0, nrow = length(uniT), ncol = 2)
chitest[, 1] <- ifelse(!is.na(oneH), 1, 0)
chitest[, 2] <- ifelse(!is.na(oneL), 1, 0)
cat("高分组唯一响应占比：", colSums(chitest)[1] / length(uniT), "\n")
cat("低分组唯一响应占比：", colSums(chitest)[2] / length(uniT), "\n")

mcnemar_result <- mcnemar.test(chitest[, 1], chitest[, 2])
print(mcnemar_result)
effect_size <- sqrt(mcnemar_result$statistic / length(uniT))
cat("McNemar 效应量 (phi)：", effect_size, "\n")

# =============================================================================
# 第七步：响应总数与 Openness 的相关性及 t 检验
# =============================================================================
sumAll <- rowSums(comb[, -1])
cor_test <- cor.test(sumAll, comb$latent)
print(cor_test)

sumLow  <- rowSums(deLow)
sumHigh <- rowSums(deHigh)
# 使用等方差 t 检验（pooled），以便后续 Cohen's d 严格成立
t_test_total <- t.test(sumHigh, sumLow, var.equal = TRUE)
print(t_test_total)

# =============================================================================
# 第八步：网络构建
# =============================================================================
jaccard_similarity <- function(mat){
  
  mat <- as.matrix(mat)
  
  p <- ncol(mat)
  
  sim <- matrix(
    NA,
    nrow = p,
    ncol = p
  )
  
  colnames(sim) <- colnames(mat)
  rownames(sim) <- colnames(mat)
  
  for(i in seq_len(p)){
    
    for(j in i:p){
      
      inter <- sum(mat[, i] == 1 & mat[, j] == 1)
      
      union <- sum(mat[, i] == 1 | mat[, j] == 1)
      
      value <- ifelse(union == 0, NA, inter / union)
      
      sim[i, j] <- value
      sim[j, i] <- value
    }
  }
  
  return(sim)
}
# 8.1 保留出现不少于2次的响应
finLow  <- finalize(deLow,  minCase = 2)
finHigh <- finalize(deHigh, minCase = 2)

# 8.2 equate 使两组拥有相同的节点集合

eq_res <- equate(finLow, finHigh)

lowMat  <- as.matrix(eq_res$finLow)
highMat <- as.matrix(eq_res$finHigh)

cat("lowMat维度:", dim(lowMat), "\n")
cat("highMat维度:", dim(highMat), "\n")

cat("共同节点数:", ncol(lowMat), "\n")

# 8.3 计算 Jaccard 相似度（节点×节点）
simLow  <- jaccard_similarity(lowMat)
simHigh <- jaccard_similarity(highMat)

# 8.4 处理 NaN 并保留共同节点
valid_low  <- !is.nan(rowSums(simLow))
valid_high <- !is.nan(rowSums(simHigh))

common_idx <- which(valid_low & valid_high)

simLow  <- simLow[common_idx, common_idx]
simHigh <- simHigh[common_idx, common_idx]

common_nodes <- colnames(lowMat)[common_idx]
cat("最终节点数：", length(common_nodes), "\n")


# 8.5 使用 TMFG 构建网络
# =============================================================================

netLow  <- NetworkToolbox::TMFG(simLow)$A
netHigh <- NetworkToolbox::TMFG(simHigh)$A

cat("netLow维度:", dim(netLow), "\n")
cat("netHigh维度:", dim(netHigh), "\n")

cat("netLow非零元素数:",
    sum(netLow > 0),
    "\n")

cat("netHigh非零元素数:",
    sum(netHigh > 0),
    "\n")
# =============================================================================
# 第九步：网络指标计算与统计检验

# 9.1 二值化
# =============================================================================

bin_low <- netLow
bin_low[bin_low != 0] <- 1

bin_high <- netHigh
bin_high[bin_high != 0] <- 1
# =============================================================================
# 9.2 转换为igraph对象
# =============================================================================

g_low <- igraph::graph_from_adjacency_matrix(
  bin_low,
  mode = "undirected",
  diag = FALSE
)

g_high <- igraph::graph_from_adjacency_matrix(
  bin_high,
  mode = "undirected",
  diag = FALSE
)
safe_modularity <- function(g){
  
  if(igraph::vcount(g) < 2)
    return(0)
  
  if(igraph::ecount(g) == 0)
    return(0)
  
  out <- tryCatch(
    {
      igraph::modularity(
        igraph::cluster_louvain(g)
      )
    },
    error = function(e){
      NA_real_
    }
  )
  
  return(out)
}

calc_metrics <- function(g){
  
  list(
    
    ASPL =
      igraph::mean_distance(
        g,
        directed = FALSE,
        unconnected = TRUE
      ),
    
    CC =
      igraph::transitivity(
        g,
        type = "global"
      ),
    
    Q =
      safe_modularity(g)
    
  )
}

metrics_low  <- calc_metrics(g_low)
metrics_high <- calc_metrics(g_high)
cat("\n=== 真实网络指标 ===\n")
cat("低分组：ASPL =", metrics_low$ASPL, ", CC =", metrics_low$CC, ", Q =", metrics_low$Q, "\n")
cat("高分组：ASPL =", metrics_high$ASPL, ", CC =", metrics_high$CC, ", Q =", metrics_high$Q, "\n")

# 9.3 随机网络 Z-test（1000次模拟）
set.seed(123)
n_nodes  <- length(common_nodes)
n_edges_low  <- sum(bin_low) / 2
n_edges_high <- sum(bin_high) / 2

simulate_random <- function(n, e, iter = 1000) {
  res <- data.frame(ASPL = numeric(iter), CC = numeric(iter), Q = numeric(iter))
  for (i in 1:iter) {
    g_rand <- sample_gnm(n, as.integer(e), directed = FALSE)
    res$ASPL[i] <- mean_distance(g_rand, unconnected = TRUE)
    res$CC[i]   <- transitivity(g_rand, type = "global")
    res$Q[i]    <- modularity(cluster_louvain(g_rand))
  }
  return(res)
}

rand_low  <- simulate_random(n_nodes, n_edges_low)
rand_high <- simulate_random(n_nodes, n_edges_high)

z_test <- function(real, dist) {
  z <- (real - mean(dist)) / sd(dist)
  p <- 2 * pnorm(-abs(z))
  return(c(Z = z, p = p))
}

cat("\n=== 随机网络 Z-test ===\n")
cat("低分组：\n")
cat("  ASPL Z =", z_test(metrics_low$ASPL, rand_low$ASPL)["Z"], 
    ", p =", z_test(metrics_low$ASPL, rand_low$ASPL)["p"], "\n")
cat("  CC   Z =", z_test(metrics_low$CC, rand_low$CC)["Z"], 
    ", p =", z_test(metrics_low$CC, rand_low$CC)["p"], "\n")
cat("  Q    Z =", z_test(metrics_low$Q, rand_low$Q)["Z"], 
    ", p =", z_test(metrics_low$Q, rand_low$Q)["p"], "\n")

cat("高分组：\n")
cat("  ASPL Z =", z_test(metrics_high$ASPL, rand_high$ASPL)["Z"], 
    ", p =", z_test(metrics_high$ASPL, rand_high$ASPL)["p"], "\n")
cat("  CC   Z =", z_test(metrics_high$CC, rand_high$CC)["Z"], 
    ", p =", z_test(metrics_high$CC, rand_high$CC)["p"], "\n")
cat("  Q    Z =", z_test(metrics_high$Q, rand_high$Q)["Z"], 
    ", p =", z_test(metrics_high$Q, rand_high$Q)["p"], "\n")

# 9.4 部分节点 Bootstrap（t-test，五个比例）
cat("\n=== 部分节点 Bootstrap（t-test，1000次） ===\n")
manual_partboot <- function(mat_low, mat_high, keep_ratio, iter = 1000) {
  nodes <- rownames(mat_low)
  keep_n <- round(length(nodes) * keep_ratio)
  res_low  <- data.frame(ASPL = numeric(iter), CC = numeric(iter), Q = numeric(iter))
  res_high <- data.frame(ASPL = numeric(iter), CC = numeric(iter), Q = numeric(iter))
  
  for (i in 1:iter) {
    samp <- sample(nodes, keep_n)
    sub_low  <- mat_low[samp, samp]; sub_low[sub_low != 0] <- 1
    sub_high <- mat_high[samp, samp]; sub_high[sub_high != 0] <- 1
    g_l <- igraph::graph_from_adjacency_matrix(
      sub_low,
      mode = "undirected"
    )
    
    g_h <- igraph::graph_from_adjacency_matrix(
      sub_high,
      mode = "undirected"
    )
    
    res_low$ASPL[i] <- igraph::mean_distance(
      g_l,
      directed = FALSE,
      unconnected = TRUE
    )
    
    res_low$CC[i] <- igraph::transitivity(
      g_l,
      type = "global"
    )
    
    res_low$Q[i] <- safe_modularity(g_l)
    
    res_high$ASPL[i] <- igraph::mean_distance(
      g_h,
      directed = FALSE,
      unconnected = TRUE
    )
    
    res_high$CC[i] <- igraph::transitivity(
      g_h,
      type = "global"
    )
    
    res_high$Q[i] <- safe_modularity(g_h)
   
  }
  return(list(low = res_low, high = res_high))
}

# 主循环
ratios <- c(0.90, 0.80, 0.70, 0.60, 0.50)
for (r in ratios) {
  cat("\n保留比例：", r*100, "%\n")
  boots <- manual_partboot(netLow, netHigh, keep_ratio = r, iter = 1000)
  
  # 使用 var.equal=TRUE 的 t.test（因为两组样本量相同，且来自同一抽样方案）
  t_aspl <- t.test(boots$high$ASPL, boots$low$ASPL, var.equal = TRUE)
  t_cc   <- t.test(boots$high$CC, boots$low$CC, var.equal = TRUE)
  t_q    <- t.test(boots$high$Q, boots$low$Q, var.equal = TRUE)
  
  # 计算 Cohen's d 使用 pooled SD 公式（精确）
  pooled_sd <- function(x, y) {
    n1 <- length(x); n2 <- length(y)
    s1 <- sd(x); s2 <- sd(y)
    sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n2-2))
  }
  d_aspl <- (mean(boots$high$ASPL) - mean(boots$low$ASPL)) / pooled_sd(boots$high$ASPL, boots$low$ASPL)
  d_cc   <- (mean(boots$high$CC) - mean(boots$low$CC)) / pooled_sd(boots$high$CC, boots$low$CC)
  d_q    <- (mean(boots$high$Q) - mean(boots$low$Q)) / pooled_sd(boots$high$Q, boots$low$Q)
  
  cat("  ASPL: t =", round(t_aspl$statistic, 3), ", df =", round(t_aspl$parameter, 1),
      ", p =", format(t_aspl$p.value, digits = 3),
      ", Cohen's d =", round(d_aspl, 2), "\n")
  cat("  CC:   t =", round(t_cc$statistic, 3), ", df =", round(t_cc$parameter, 1),
      ", p =", format(t_cc$p.value, digits = 3),
      ", Cohen's d =", round(d_cc, 2), "\n")
  cat("  Q:    t =", round(t_q$statistic, 3), ", df =", round(t_q$parameter, 1),
      ", p =", format(t_q$p.value, digits = 3),
      ", Cohen's d =", round(d_q, 2), "\n")
}

# =============================================================================
# 第十步：导出 Cytoscape 文件（边列表 + GraphML）
# =============================================================================
g_high_cyto <- graph_from_adjacency_matrix(netHigh, mode = "undirected", weighted = TRUE)
g_low_cyto  <- graph_from_adjacency_matrix(netLow,  mode = "undirected", weighted = TRUE)

high_edges <- as_data_frame(g_high_cyto, what = "edges")
low_edges  <- as_data_frame(g_low_cyto,  what = "edges")

write.csv(high_edges, "high_open_edges.csv", row.names = FALSE)
write.csv(low_edges,  "low_open_edges.csv",  row.names = FALSE)

write.graph(g_high_cyto, "high_open.graphml", format = "graphml")
write.graph(g_low_cyto,  "low_open.graphml",  format = "graphml")

cat("\n✅ 所有分析完成！Cytoscape 文件已导出。\n")

