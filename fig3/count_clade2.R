library(ggplot2)
library(dplyr)
library(tidyr)

sRefSp <- "RhinogobiusFormosanus";
nFDRCutoff <- 0.05
nRawPCutoff <- 1
setwd("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/hyphy_1016/GO")
datRelax_clade2 <- read.table("../relax/2clade/sum_Clade2.txt", header=F, sep="\t", fill = T, quote = "")
datRelax_clade5 <- read.table("../relax/5clade/sum_Clade5.txt", header=F, sep="\t", fill = T, quote = "")
datRelax_clade7 <- read.table("../relax/7clade/sum_Clade7.txt", header=F, sep="\t", fill = T, quote = "")
datRelax_clade2$fdr <- p.adjust(datRelax_clade2$V5, method = "fdr")
datRelax_clade5$fdr <- p.adjust(datRelax_clade5$V5, method = "fdr")
datRelax_clade7$fdr <- p.adjust(datRelax_clade7$V5, method = "fdr")
datGroupID2TransID <- read.table("../exportAln_genespace/groupid2transid.txt", sep="\t", header=T)
colnames(datGroupID2TransID) <- c("OrthoID", 'Gene');
datRelax_clade2 <- merge(datRelax_clade2, datGroupID2TransID, by.x = "V1", by.y="OrthoID", all.x=T, all.y=F)
datRelax_clade5 <- merge(datRelax_clade5, datGroupID2TransID, by.x = "V1", by.y="OrthoID", all.x=T, all.y=F)
datRelax_clade7 <- merge(datRelax_clade7, datGroupID2TransID, by.x = "V1", by.y="OrthoID", all.x=T, all.y=F)

# 添加 group 列，根据 k 值判断是 relaxed 或 intensified
datRelax_clade7 <- datRelax_clade7 %>%
  mutate(group = ifelse(V9 < 1, "Relaxed", "Intensified"))

# 设置 bin 宽度
bin_width <- 0.05
breaks <- seq(0, 1, by = bin_width)

# 添加 bin，统计每组每个 bin 的数量
datRelax_clade7_bins <- datRelax_clade7 %>%
  mutate(bin = cut(fdr, breaks = breaks, include.lowest = TRUE, right = FALSE)) %>%
  group_by(group, bin) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(
    bin_start = as.numeric(gsub("\\[|\\)|\\[|\\]", "", sub(",.*", "", as.character(bin)))),
    bin_mid = bin_start + bin_width / 2,
    count = ifelse(group == "Intensified", -count, count)
  )
datRelax_clade7$bin <- cut(datRelax_clade7$fdr, breaks = breaks, include.lowest = TRUE, right = FALSE)
# 统计 Intensified / Relaxed 各 bin 的数量
datRelax_clade7_wide <- datRelax_clade7 %>%
  group_by(bin, group) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = count, values_fill = 0) %>%
  mutate(
    log_ratio = log10((Relaxed) / (Intensified)),  # +1 防止除以0
    bin_start = as.numeric(gsub("\\[|\\)|\\[|\\]", "", sub(",.*", "", as.character(bin)))),
    bin_mid = bin_start + bin_width / 2
  )

scaling_factor <- max(abs(datRelax_clade7_bins$count)) / max(abs(datRelax_clade7_wide$log_ratio)) * 0.8
ggplot() +
  # 主体条形图
  geom_bar(data = datRelax_clade7_bins,
           aes(x = bin_mid, y = count, fill = group),
           stat = "identity", position = "identity", width = bin_width, alpha = 0.7) +
  
  # 添加 log(Intensified / Relaxed) 曲线（右轴）
  geom_line(data = datRelax_clade7_wide,
            aes(x = bin_mid, y = log_ratio * scaling_factor),
            color = "black", size = 0.8) +
  
  scale_fill_manual(values = c("Intensified" = rgb(183, 034, 048, max = 255), "Relaxed" = rgb(049, 124, 183, max = 255))) +
  
  # 重新设计双 y 轴
  scale_y_continuous(
    name = "Count",
    sec.axis = sec_axis(
      trans = ~ . / scaling_factor,
      name = "Log10(relaxed : intensified)"
    ),
    breaks = pretty(c(-max(abs(datRelax_clade7_bins$count)), max(abs(datRelax_clade7_bins$count))), n = 5),
    labels = function(x) abs(x)
  ) +
  
  labs(x = "p value (FDR)", fill = NULL) +
  theme_classic(base_size = 10) +  # 基础字体稍小，适合文章图
  theme(
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.4),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 11, face = "bold"),
    axis.title.y.right = element_text(color = "black", size = 11, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    panel.grid = element_blank(),        # 无网格线
    panel.background = element_blank(),  # 白色背景
    plot.background = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    plot.margin = margin(t = 15, r = 5, b = 5, l = 5)  # top, right, bottom, left
  )+ 
  ggtitle("\nclade7") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  )

