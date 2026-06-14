library(dplyr)

clades <- list(
  Clade2 = c("Z5_R-niger_Zhejianglongquanshizhulongzhen", "Z2_R-aporus_Zhejiangpananxianshuangxizhen",
             "A8_R-cf-wuyiensis_Anhuiqimenlikouzhen", "F6_R-cf-wuyiensis_Fujianwuyishan", "J2_R-sp_Jiangxiwuyuan",
             "Z13_R-sp_Zhejiangkaihuaxiansuzhuangzhen", "Z1_R-sp_Zhejianglongquanshizhulongzhen",
             "J1_R-sp_Jiangxipingxiangshianyuanqu", "A1_R-sp_Anhuidongzhixiannixizhen",
             "A3_R-sp_Anhuiqimendatanxiang", "A6_R-lentiginis_Anhuishitai",
             "Z16_R-cf-lentiginis_Zhejiangchangshanxian", "Z14_R-lentiginis_Zhejiangyuyaoliangnongzhen",
             "Z19_R-aporus_Zhejiangwuyixianchengqu", "Z3_R-wuyiensis_Zhejianglongquanshijinxizhen",
             "Z20_R-aporus_Zhejiangwuyixianxilianxiang"),
  Clade5 = c("S1_R-szechuanensis_Sichuan", "HB1_R-houheensis_Hubeixuandongxian",
             "A7_R-sp_Anhuishitai", "Z12_R-davidi_Zhejiangchangshanxian",
             "Z7_R-davidi_Zhejianganjixiantianhuangpingzhen", "Z11_R-sp_Zhejiangchangshanxian",
             "Z21_R-sp_Zhejiangchangshanxiansianzhen"),
  Clade7 = c("HN2_R-wanchuangensis_Hainanqionghaishihuishanzhen", "HN5_R-linshuiensis_Hainanbaotingshilingzhen",
             "GX2_R-sp_Guangxifangchenggangmaluzhen",
             "HN1_R-nandujiangensis_Hainanhaikoushijiuzhouzhen",
             "GD1_R-duospilus_Guangdonglongmen", "GX3_R-sp_Guangxifangchenggangmaluzhen",
             "F7a_R-sp_Fujianyunxiao", "F7b_R-sp_Fujianyunxiao", "F2_R-xianshuiensis_Fujianputianshi",
             "Z4_R-cf-changtingensis_Zhejianglongquanshizhulongz", "F4_R-reticulatus_Fujianfuzhou",
             "F3_R-reticulatus_Fujianpuchengshifulingzhen", "H1_R-sp_Hunanlianyuanshimaotangzhen",
             "H2_R-sp_Hunanloudilouxingqu", "Z9_R-sp_Zhejiangkaihuaxiansuzhuangzhen",
             "Z8_R-immaculatus_Zhejiangtongluxianlucixiang", "Z10_R-sp_Zhejianglongquanshizhulongzhen")
)
Amphidromous <- c('L1_R-cf-brunneus_Liaoningkuandian',
                  'HN3_R-leavelli_Hainanzhanzhoudachengzhen', 'HN4_R-leavelli_Hainanbaotingshilingzhen',
                  'Z15_R-leavelli_Zhejiangchangshanxian',
                  'Z17_R-cliffordpopei_Zhejiangchangshanxian', 'Z18_R-cliffordpopeiR-leavelli_Zhejiangchangshanxia', 'Z6_R-giurinus_Zhejianglongquanshizhulongzhen')



setwd("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/cafe5/family_test0401")
input_file <- "orthogroups.genecount.txt"
input_file <- "Rformosanus.orthogroups.genecount.txt"
data <- read.table(input_file, header = TRUE, sep = "\t", check.names = FALSE)

# 删除非数值列（如 "Desc"）
data <- data[,-1]
data <- data[,-c(2:9)]
data <- data[rowSums(data[,-1]) != 0, ]
rownames(data) <- data[[1]]
data <- data[, -1]


# 显著性检验函数
run_significance_test <- function(data, clade_name, clade_members, comparison_group, output_file, output_file2) {
  # 初始化存储结果的表格
  result <- data.frame(Orthogroup = rownames(data), P_value = NA, Clade_Avg = NA, Amphidromous_Avg = NA, Expansion_Status = NA)
  
  # 循环每个 Orthogroup 进行检验
  for (i in 1:nrow(data)) {
    clade_counts <- as.numeric(data[i, clade_members, drop = FALSE])  # Clade 的数据
    comparison_counts <- as.numeric(data[i, comparison_group, drop = FALSE])  # 对比组的数据
    
    # 计算平均值
    clade_avg <- mean(clade_counts)
    amphidromous_avg <- mean(comparison_counts)
    
    # 检查是否存在非零值，确保能够进行检验
    if (length(clade_counts) > 0 && length(comparison_counts) > 0) {
      # 使用 Wilcoxon 检验（非参数检验）
      test <- wilcox.test(clade_counts, comparison_counts, exact = FALSE)
      result$P_value[i] <- test$p.value
    }
    # 存储平均值和扩张/收缩信息
    result$Clade_Avg[i] <- clade_avg
    result$Amphidromous_Avg[i] <- amphidromous_avg
    result$Expansion_Status[i] <- ifelse(clade_avg > amphidromous_avg, "Expanded", "Contracted")
    
  }
  
  # 添加 FDR 校正列
  result$FDR <- p.adjust(result$P_value, method = "fdr")
  
  # 输出显著性结果
  significant_results <- result %>%
    filter(!is.na(FDR) & FDR < 0.05) %>%
    arrange(FDR)
  
  # 保存结果到文件
  write.csv(significant_results, output_file, row.names = FALSE)
  write.csv(result, output_file2, row.names = FALSE, quote = FALSE)
  
  # 返回结果
  return(significant_results)
}

# 运行显著性检验
clade2_results <- run_significance_test(data, "Clade2", clades$Clade2, Amphidromous, "clade2_vs_amphidromous_fdr0.05_results.csv", "clade2_vs_amphidromous_all_results.csv")
clade5_results <- run_significance_test(data, "Clade5", clades$Clade5, Amphidromous, "clade5_vs_amphidromous_fdr0.05_results.csv", "clade5_vs_amphidromous_all_results.csv")
clade7_results <- run_significance_test(data, "Clade7", clades$Clade7, Amphidromous, "clade7_vs_amphidromous_fdr0.05_results.csv", "clade7_vs_amphidromous_all_results.csv")

# 打印显著性结果
print("Clade2 vs Amphidromous significant results:")
print(clade2_results)

print("Clade5 vs Amphidromous significant results:")
print(clade5_results)

print("Clade7 vs Amphidromous significant results:")
print(clade7_results)
r
