setwd("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/snpEff/snpeff_0413")

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ggpubr) # 建议安装这个包，用于方便地添加统计标注，如果没有可以手动计算

file_high <- "high_gene_species_counts.txt"
file_mod  <- "moderate_gene_species_counts.txt"
file_low  <- "low_gene_species_counts.txt"
species_file <- "get_clade_inf.txt"

output_prefix <- "Mutation_Load_With_Stats"

if (!file.exists(file_high) || !file.exists(file_mod) || !file.exists(file_low) || !file.exists(species_file)) {
  stop("错误: 找不到输入文件，请检查目录！")
}

# 1. 数据处理函数
process_load <- function(file_path, impact_label) {
  df <- suppressMessages(read_tsv(file_path, show_col_types = FALSE))
  if (!("gene_id" %in% colnames(df))) colnames(df)[1] <- "gene_id"
  species_data <- df %>% select(-gene_id)
  
  species_vals <- colSums(species_data, na.rm = TRUE)
  
  data.frame(
    species = names(species_vals),
    load_val = as.numeric(species_vals),
    impact_type = impact_label,
    stringsAsFactors = FALSE
  )
}

# 2. 读取并合并数据
df_high <- process_load(file_high, "High")
df_mod  <- process_load(file_mod,  "Moderate")
df_low  <- process_load(file_low,  "Low")

all_data <- bind_rows(df_high, df_mod, df_low)
all_data$impact_type <- factor(all_data$impact_type, levels = c("Low", "Moderate", "High"))

# 3. 物种信息匹配
sp_info <- read.table(species_file, header = TRUE, sep = "", stringsAsFactors = FALSE)
sp_map <- sp_info %>% select(species, clade, lifestyle)

merged_df <- left_join(all_data, sp_map, by = "species")

# 4. 筛选并分组
plot_data <- merged_df %>%
  mutate(
    clade_lower = tolower(clade),
    life_lower = tolower(lifestyle),
    group_label = case_when(
      grepl("amphidromous", life_lower) ~ "Amphidromous",
      clade_lower %in% c("clade2", "clade5", "clade7") ~ "Landlocked",
      TRUE ~ "Others"
    )
  ) %>%
  filter(group_label != "Others")

plot_data$group_label <- factor(plot_data$group_label, levels = c("Amphidromous", "Landlocked"))

# ==========================================
# 5. 统计检验 (Wilcoxon Test)
# ==========================================
stats_results <- plot_data %>%
  group_by(impact_type) %>%
  summarise(
    p_val = wilcox.test(load_val ~ group_label)$p.value,
    .groups = 'drop'
  ) %>%
  mutate(
    p_label = case_when(
      p_val < 0.001 ~ "***",
      p_val < 0.01  ~ "**",
      p_val < 0.05  ~ "*",
      TRUE          ~ "ns"
    )
  )

print("统计检验结果 (Amphidromous vs Landlocked):")
print(stats_results)

# 6. 计算绘图摘要
summary_df <- plot_data %>%
  group_by(impact_type, group_label) %>%
  summarise(
    mean_value = mean(load_val),
    sd_value = sd(load_val),
    n = n(),
    se_value = sd_value / sqrt(n),
    max_top = mean_value + se_value, # 用于确定显著性标签的高度
    .groups = 'drop'
  )

# 合并 P 值标签到摘要中以便绘图
label_df <- summary_df %>%
  group_by(impact_type) %>%
  summarise(y_pos = max(max_top) * 1.05) %>% # 标签显示在柱子上方 5% 的位置
  left_join(stats_results, by = "impact_type")

# 7. 绘图
soft_colors <- c("Amphidromous" = "#8DA0CB", "Landlocked" = "#FC8D62")
bar_width <- 0.6
dodge_width <- 0.7

p <- ggplot(summary_df, aes(x = impact_type, y = mean_value, fill = group_label)) +
  # 柱状图
  geom_bar(stat = "identity",
           position = position_dodge(width = dodge_width),
           width = bar_width) +
  # 误差棒
  geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value),
                position = position_dodge(width = dodge_width),
                width = 0.2, size = 0.5, color = "black") +
  
  # 添加显著性标记 (星号)
  geom_text(data = label_df, 
            aes(x = impact_type, y = y_pos, label = p_label), 
            inherit.aes = FALSE, 
            size = 5, fontface = "bold", vjust = 0) +
  
  scale_fill_manual(values = soft_colors, name = "Ecological Type") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + # 顶部留出 15% 空间放标签
  theme_classic() +
  labs(
    title = "Genetic Load Comparison by Ecological Type",
    subtitle = "Stats: Wilcoxon Rank Sum Test (* p<0.05, ** p<0.01, *** p<0.001)",
    y = "Average Total Mutations",
    x = "Mutation Impact Category"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 9, color = "gray20", hjust = 0.5),
    axis.text.x = element_text(size = 11, face = "bold", color = "black"),
    legend.position = "top"
  )

# 保存
out_pdf <- paste0(output_prefix, ".pdf")
ggsave(out_pdf, p, width = 7, height = 8)
cat(paste0("分析完成，结果已保存至: ", out_pdf, "\n"))
