setwd("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/snpEff")
#洄游型7，clade2 16，clade5 7，clade7 17, 取0.6
#frame shift
frameshift_data <- read.table("frameshift_gene_stats.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
clade2_frameshift <- frameshift_data[frameshift_data$clade2_count > 6 & frameshift_data$amphidromous_count == 0, ]
clade5_frameshift <- frameshift_data[frameshift_data$clade5_count > 2 & frameshift_data$amphidromous_count == 0, ]
clade7_frameshift <- frameshift_data[frameshift_data$clade7_count > 6 & frameshift_data$amphidromous_count == 0, ]
arrTargetGenes_clade2 <- clade2_frameshift$gene_id
arrTargetGenes_clade5 <- clade5_frameshift$gene_id
arrTargetGenes_clade7 <- clade7_frameshift$gene_id
#high impact
clade2 <- read.table("highimpact_clade2_vs_amphidromous.txt", 
                     header = TRUE, sep = "\t", stringsAsFactors = FALSE)
clade5 <- read.table("highimpact_clade5_vs_amphidromous.txt", 
                     header = TRUE, sep = "\t", stringsAsFactors = FALSE)
clade7 <- read.table("highimpact_clade7_vs_amphidromous.txt", 
                     header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sig_clade2 <- subset(clade2, p_value < 0.05, select = gene_id)
sig_clade5 <- subset(clade5, p_value < 0.05, select = gene_id)
sig_clade7 <- subset(clade7, p_value < 0.05, select = gene_id)
all_clade2 <- clade2$gene_id
all_clade5 <- clade5$gene_id
all_clade7 <- clade7$gene_id
all_genes <- intersect(all_clade2, all_clade5)
all_genes <- intersect(all_genes, all_clade7)
arrTargetGenes_clade2 <- sig_clade2$gene_id
arrTargetGenes_clade5 <- sig_clade5$gene_id
arrTargetGenes_clade7 <- sig_clade7$gene_id
arrTargetGenes_clade25 <- intersect(arrTargetGenes_clade2, arrTargetGenes_clade5)
arrTargetGenes_clade_all <- intersect(arrTargetGenes_clade25, arrTargetGenes_clade7)
clade2$change <- clade2$clade_mean - clade2$amphidromous_mean
clade2 <- subset(clade2, p_value < 0.05 & change > 0)
clade5$change <- clade5$clade_mean - clade5$amphidromous_mean
clade5 <- subset(clade5, p_value < 0.05 & change > 0)
clade7$change <- clade7$clade_mean - clade7$amphidromous_mean
clade7 <- subset(clade7, p_value < 0.05 & change > 0)
clade2 <- fnAddGeneSymbols(clade2)
clade5 <- fnAddGeneSymbols(clade5)
clade7 <- fnAddGeneSymbols(clade7)
clade2_renamed <- clade2 %>% rename_with(~ ifelse(.x == "gene_id", .x, paste0(.x, "_clade2")))
clade5_renamed <- clade5 %>% rename_with(~ ifelse(.x == "gene_id", .x, paste0(.x, "_clade5")))
clade7_renamed <- clade7 %>% rename_with(~ ifelse(.x == "gene_id", .x, paste0(.x, "_clade7")))
merged_data <- clade2_renamed %>%
  inner_join(clade5_renamed, by = "gene_id") %>%
  inner_join(clade7_renamed, by = "gene_id")
sorted_data <- merged_data %>%
  mutate(mean_change = rowMeans(select(., contains("change")), na.rm = TRUE)) %>%
  arrange(desc(mean_change))
df_top30 <- sorted_data %>%
  filter(!is.na(genesymbols_clade7) & genesymbols_clade7 != "") %>%
  head(30)
# 提取 geneid 与 symbol 对应表
gene_map <- df_top30[, c("gene_id", "genesymbols_clade7")]


fig2 <- read.table("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/plot/plot2/fig2.csv", header = TRUE, sep = ",", check.names = FALSE)
snpeff_data <- read.table("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/snpEff/high_gene_species_counts.txt", header = TRUE, sep = "\t", check.names = FALSE)
snpeff_data <- t(snpeff_data)
colnames(snpeff_data) <- snpeff_data[1,]
snpeff_data <- snpeff_data[-1,]
snpeff_data <- as.data.frame(snpeff_data)
snpeff_data$"Node IDs" <- rownames(snpeff_data)
library(dplyr)
snpeff_data <- snpeff_data %>% select("Node IDs", everything())
# 从 df2 中选择第一列 + 这30个基因列
# 确保 geneid 匹配的是列名
genes_to_keep <- gene_map$gene_id
cols_to_keep <- c(colnames(snpeff_data)[1], intersect(genes_to_keep, colnames(snpeff_data)))
df_selected <- snpeff_data[, cols_to_keep, drop = FALSE]
# 替换列名为对应 symbol（保留第一列原名）
new_colnames <- colnames(df_selected)
symbol_map <- setNames(gene_map$genesymbols_clade7, gene_map$gene_id)
for (i in seq_along(new_colnames)) {
  if (new_colnames[i] %in% names(symbol_map)) {
    new_colnames[i] <- symbol_map[[new_colnames[i]]]
  }
}
colnames(df_selected) <- new_colnames
print(df_selected)

merged_table <- merge(fig2, df_selected, by = "Node IDs", all = TRUE)
fnAddGeneSymbols <- function(datOut, geneID_col = "gene_id", symbol_col = "genesymbols") {
  # 添加调试信息
  cat("函数开始执行...\n")
  cat("输入数据框维度:", dim(datOut), "\n")
  cat("基因ID列名:", geneID_col, "\n")
  
  # 检查输入数据
  if (nrow(datOut) == 0) {
    cat("输入数据框为空\n")
    return(datOut)
  }
  
  # 检查必需的列是否存在
  if (!geneID_col %in% colnames(datOut)) {
    stop(paste("错误: 数据框中不存在列", geneID_col, "。可用列名:", paste(colnames(datOut), collapse = ", ")))
  }
  
  # 检查datSymbol是否存在
  if (!exists("datSymbol")) {
    stop("错误: 未找到基因符号映射表 'datSymbol'")
  }
  
  cat("datSymbol 维度:", dim(datSymbol), "\n")
  cat("datSymbol 列名:", paste(colnames(datSymbol), collapse = ", "), "\n")
  
  # 检查datSymbol的列结构
  if (!"Gene" %in% colnames(datSymbol)) {
    stop("错误: datSymbol 中缺少 'Gene' 列")
  }
  
  if (ncol(datSymbol) < 1) {
    stop("错误: datSymbol 至少需要有一列包含基因符号")
  }
  
  # 初始化基因符号向量
  arrSymbolStr <- character(nrow(datOut))  # 预分配长度，提高性能
  
  # 遍历每一行
  for(nRow in 1:nrow(datOut)) {
    # 获取当前行的基因ID，确保为字符型
    geneID_value <- as.character(datOut[[geneID_col]][nRow])
    
    # 检查geneID值
    if (is.na(geneID_value) || geneID_value == "") {
      arrSymbolStr[nRow] <- ""
      next
    }
    
    # 用'/'分割成多个基因ID
    arrPGIDs <- unlist(strsplit(geneID_value, '/'))
    
    # 从datSymbol中查找对应的基因符号
    arrGeneSymbols <- datSymbol[datSymbol$Gene %in% arrPGIDs, 1]
    
    # 如果找到对应的基因符号，用'/'连接
    if (length(arrGeneSymbols) > 0) {
      arrSymbolStr[nRow] <- paste(as.character(arrGeneSymbols), collapse = '/')
    } else {
      arrSymbolStr[nRow] <- ""
    }
  }
  
  # 添加基因符号列到数据框
  datOut[[symbol_col]] <- arrSymbolStr
  
  cat("函数执行完成，添加了列:", symbol_col, "\n")
  cat("结果数据框维度:", dim(datOut), "\n")
  
  # 返回修改后的数据框
  return(datOut)
}
write.table(merged_table, "/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/plot/fig4/final_fig4data.csv", sep = ",", col.names = T, row.names = F, quote = F)







df <- read.table("SV.polarized.out.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(df)[1:8] <- c("CHROM", "POS", "REF", "ALT", "GENE", "FUNCTION", "IMPACT", "AA_CHANGE")
df_high <- df[df$IMPACT == "HIGH", ]
# 检查哪些列全是NA
na_columns <- colnames(df_high)[colSums(!is.na(df_high)) == 0]
df <- read.table("high_gene_species_counts.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
num_df <- df[, -1]
# 找出全为0的列
zero_cols <- names(num_df)[colSums(num_df) == 0]
zero_cols