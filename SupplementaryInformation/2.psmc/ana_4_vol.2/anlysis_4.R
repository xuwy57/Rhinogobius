library(ape)
library(geiger)
library(dplyr)

# 1. 读取树
phy <- read.tree("../tree.tre")

# 2. 读取物种名单 
meta <- read.table("../get_clade_inf.txt", header=TRUE, sep="")
rownames(meta) <- meta$species

# 3. 检查匹配情况
check <- name.check(phy, meta)

if(length(check$tree_not_data) > 0) {
  cat("以下物种在树里有，但数据里没有 (会被自动丢弃，包括外群):\n")
  print(head(check$tree_not_data))
}

if(length(check$data_not_tree) > 0) {
  cat("\n⚠️ 警告：以下物种在数据里有，但在树里找不到 (会被丢弃):\n")
  print(check$data_not_tree)
  cat("请检查名称拼写是否完全一致！(例如 Rhinogobius_formosanus)\n")
} else {
  cat("\n完美！所有数据中的物种都在树里找到了。\n")
}


#!/usr/bin/env Rscript
  
# ==============================================================================
# 0. 环境设置与包加载
# ==============================================================================
packages <- c("ape", "nlme", "dplyr", "readr", "geiger")
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

suppressPackageStartupMessages({
  library(ape)    # 系统发育树处理
  library(nlme)   # 广义最小二乘法 (GLS)
  library(geiger) # 数据与树的匹配
  library(dplyr)
  library(readr)
})

# ==============================================================================
# 1. 配置文件路径 
# ==============================================================================
setwd("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/plot/psmc/ana_4_vol.2") 

# 输入文件
load_file    <- "high_gene_species_counts.txt"   # 累积有害突变文件
sample_file  <- "../get_clade_inf.txt"              # 物种信息
psmc_folder  <- "../psmc_result/"                 # PSMC 结果目录
tree_file    <- "../tree.tre" 

# 参数
nMu <- 2.6440E-9
# 计算 Ne 的时间窗口 (Years Ago)
# 选取近期窗口，因为这与当前的突变负荷最相关
time_window <- c(20000, 500000) 

# ==============================================================================
# 2. 准备数据：Mutation Load (Y)
# ==============================================================================
cat("正在计算 Mutation Load...\n")
df_load <- read_tsv(load_file, show_col_types = FALSE)
if (!("gene_id" %in% colnames(df_load))) colnames(df_load)[1] <- "gene_id"

# 计算每个物种的总负荷 (Sum)
load_sums <- colSums(df_load %>% select(-gene_id), na.rm = TRUE)
data_load <- data.frame(species = names(load_sums), Load = as.numeric(load_sums))

# ==============================================================================
# 3. 准备数据：Harmonic Mean Ne (X1)
# ==============================================================================
cat("正在计算 Harmonic Mean Ne (窗口: ", paste(time_window, collapse="-"), ")... \n")

# 读取物种列表
meta_df <- read.table(sample_file, header = TRUE, sep = "", stringsAsFactors = FALSE)
meta_df <- meta_df %>% filter(!is.na(lifestyle))

ne_results <- list()

for (i in 1:nrow(meta_df)) {
  sPop <- meta_df$species[i]
  
  path1 <- paste0(psmc_folder, sPop, ".bam.final.txt")
  path2 <- paste0(psmc_folder, sPop, ".final.txt")
  sInFile <- NULL
  if (file.exists(path1)) {
    sInFile <- path1
  } else if (file.exists(path2)) {
    sInFile <- path2
  }

  if (!is.null(sInFile)) {
    tryCatch({
      dat <- read.table(sInFile, header = TRUE, sep = "\t")
      dat$Years <- dat$left_time_boundary / nMu
      dat$Ne    <- (1 / dat$lambda) / (2 * nMu)
      dat_win <- dat %>% filter(Years >= time_window[1] & Years <= time_window[2])
      if (nrow(dat_win) > 0) {
        h_mean_ne <- 1 / mean(1 / dat_win$Ne)
        ne_results[[sPop]] <- h_mean_ne
      } else {
        ne_results[[sPop]] <- NA
        cat("警告: 物种", sPop, "在指定时间窗口内无有效数据点。\n")
      }
    }, error = function(e) {
      cat("错误: 读取或处理文件失败:", sInFile, "\n")
      ne_results[[sPop]] <- NA
    })
    
  } else {
    cat("警告: 找不到文件 ->", sPop, "\n")
    ne_results[[sPop]] <- NA
  }
}
ne_results <- ne_results[!is.na(ne_results)]

if (length(ne_results) == 0) {
  stop("错误：没有计算出任何有效的 Ne 值，请检查 PSMC 结果路径和时间窗口设置！")
}

data_ne <- data.frame(species = names(ne_results), Ne = unlist(ne_results))
# ==============================================================================
# 4. 准备数据：Habitat (X2) 与合并
# ==============================================================================
data_habitat <- meta_df %>% 
  select(species, lifestyle, clade) %>%
  mutate(
    Habitat = ifelse(grepl("amphidromous", lifestyle, ignore.case=T), "Amphidromous", "Landlocked")
  )

full_data <- data_load %>%
  inner_join(data_ne, by = "species") %>%
  inner_join(data_habitat, by = "species") %>%
  filter(!is.na(Ne)) # 去掉没有 Ne 数据的

rownames(full_data) <- full_data$species

cat("数据准备完成，共", nrow(full_data), "个物种。\n")

# ==============================================================================
# 5. 系统发育树处理
# ==============================================================================
cat("正在读取并修剪系统发育树...\n")
if (!file.exists(tree_file)) stop("找不到树文件，请检查脚本中的 tree_file 路径！")

phy_tree <- read.tree(tree_file)

cleaned <- treedata(phy_tree, full_data, sort = TRUE, warnings = FALSE)

phy_tree_clean <- cleaned$phy
data_clean <- as.data.frame(cleaned$data)

data_clean$Load <- as.numeric(data_clean$Load)
data_clean$Ne   <- as.numeric(data_clean$Ne)
data_clean$Habitat <- factor(data_clean$Habitat, levels = c("Amphidromous", "Landlocked"))

cat("最终纳入分析的物种数:", nrow(data_clean), "\n")

# ==============================================================================
# 6. PGLS 模型构建
# ==============================================================================
cat("\n=== 运行 PGLS 模型 ===\n")
cat("Formula: Load ~ log10(Ne) + Habitat\n")
cat("Correlation: Pagel's Lambda (自动估计系统发育信号)\n\n")

# 核心模型
# correlation = corPagel(1, ...) 表示假设 lambda=1 开始拟合，允许优化
model_pgls <- gls(Load ~ log10(Ne) + Habitat, 
                  correlation = corPagel(1, phy = phy_tree_clean, fixed = FALSE),
                  data = data_clean)

# 输出摘要
sum_res <- summary(model_pgls)
print(sum_res)

# ==============================================================================
# 7. 结果解析与导出
# ==============================================================================
coefs <- sum_res$tTable
p_val_ne <- coefs["log10(Ne)", "p-value"]
p_val_hab <- coefs["HabitatLandlocked", "p-value"]
est_hab <- coefs["HabitatLandlocked", "Value"]

cat("\n--------------------------------------------------------------\n")
cat("分析结论:\n")
cat("1. 种群大小 (logNe) 的 P值:", sprintf("%.4f", p_val_ne), "\n")
if(p_val_ne < 0.05) {
  cat("   -> 显著! 说明漂变确实影响了突变积累。\n")
} else {
  cat("   -> 不显著。可能样本量不足或效应被掩盖。\n")
}

cat("2. 生境类型 (Habitat) 的 P值:", sprintf("%.4f", p_val_hab), "\n")
cat("   -> 系数 (Coefficient):", sprintf("%.4f", est_hab), "\n")

if (p_val_hab < 0.05 && est_hab > 0) {
  cat("\n[重要发现] !!!\n")
  cat("即使控制了 Ne 和亲缘关系，Landlocked 生境依然与更高的 Mutation Load 显著正相关。\n")
  cat("这支持了‘生境转变(放松选择/适应代价)独立于种群变小，导致了额外的遗传负荷’的假说。\n")
} else if (p_val_hab >= 0.05) {
  cat("\n[结论]\n")
  cat("生境类型没有显著的独立效应。这意味着 Landlocked 物种的高遗传负荷\n")
  cat("主要完全可以通过它们较小的种群大小(漂变)来解释，不需要引入放松选择假说。\n")
}
cat("--------------------------------------------------------------\n")

# 保存简单的统计结果
write.csv(coefs, "PGLS_Result_Coefficients.csv")
cat("详细统计表已保存至: PGLS_Result_Coefficients.csv\n")


# ==============================================================================
# 8. 可视化 PGLS 结果 
# ==============================================================================
library(ggplot2)

# 1. 提取 PGLS 模型的系数，用于画回归线
# PGLS 回归方程: Y = Intercept + Slope * X + Habitat_Effect
# Amphidromous (基准): Y = Intercept + Slope * X
# Landlocked:          Y = (Intercept + Habitat_Effect) + Slope * X

coeffs <- coef(model_pgls)
intercept_base <- coeffs["(Intercept)"]
slope_ne <- coeffs["log10(Ne)"]
habitat_effect <- coeffs["HabitatLandlocked"]

# 2. 准备绘图数据 (确保 logNe 已计算)
plot_data <- data_clean
plot_data$logNe <- log10(plot_data$Ne)

# 3. 计算预测线的数据框
# 需要画两条平行的线，斜率都是 slope_ne
# 截距分别是 intercept_base 和 intercept_base + habitat_effect

lines_data <- data.frame(
  Habitat = c("Amphidromous", "Landlocked"),
  Intercept = c(intercept_base, intercept_base + habitat_effect),
  Slope = c(slope_ne, slope_ne)
)

# 4. 绘图
# 颜色设置
my_colors <- c("Amphidromous" = "#377EB8", "Landlocked" = "#E41A1C")

p <- ggplot(plot_data, aes(x = logNe, y = Load, color = Habitat)) +
  # A. 散点
  geom_point(size = 3, alpha = 0.7) +
  
  # B. PGLS 拟合线 (注意：这里不能用 geom_smooth(method=lm)，因为那是 OLS 线)
  # 我们使用 geom_abline 手动画出 PGLS 算出的线
  geom_abline(data = lines_data, 
              aes(intercept = Intercept, slope = Slope, color = Habitat, linetype = Habitat),
              size = 1) +
  
  # C. 展示 "生境效应" 的辅助箭头 
  # 找一个 X 位置 (X轴的中位数) 画垂直双向箭头
  annotate("segment", 
           x = median(plot_data$logNe), xend = median(plot_data$logNe),
           y = (intercept_base + slope_ne * median(plot_data$logNe)),
           yend = (intercept_base + habitat_effect + slope_ne * median(plot_data$logNe)),
           arrow = arrow(ends = "both", length = unit(0.2, "cm")),
           color = "black", alpha = 0.6) +
  annotate("text", 
           x = median(plot_data$logNe) + 0.05, 
           y = (intercept_base + habitat_effect/2 + slope_ne * median(plot_data$logNe)),
           label = paste0("Habitat Effect\n+", round(habitat_effect), " mutations"),
           hjust = 0, color = "black", fontface = "italic") +
  
  # D. 美化
  scale_color_manual(values = my_colors) +
  scale_linetype_manual(values = c("solid", "solid")) + # 两条都是实线
  theme_classic() +
  labs(
    title = "PGLS Regression: Habitat Effect on Genetic Load",
    subtitle = "Landlocked species carry higher load even at the same Ne",
    x = expression(bold(log[10](Harmonic~Mean~N[e]))),
    y = "Accumulated Mutation Load (High Impact)",
    caption = paste0("PGLS Model: Load ~ log(Ne) + Habitat\n",
                     "Ne p-value < 0.0001 | Habitat p-value < 0.0001\n",
                     "Lambda = ", round(summary(model_pgls)$modelStruct$corStruct[1], 3))
  ) +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10, color = "black")
  )

# 5. 保存
ggsave("PGLS_Regression_Plot.pdf", p, width = 8, height = 6)

print(p)
