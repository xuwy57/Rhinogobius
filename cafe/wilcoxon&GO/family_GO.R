setwd("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/cafe5/family_test0401");
library(../family_test0125/org.Rformosanus.clade2.eg.db)
library(org.Rformosanus.clade5.eg.db)
library(org.Rformosanus.clade7.eg.db)
library(clusterProfiler)
library(stringr)
library(readr)

clade2_significant_results <- read_csv("clade2_vs_amphidromous_fdr0.2_results.csv")
clade5_significant_results <- read_csv("clade5_vs_amphidromous_fdr0.2_results.csv")
clade7_significant_results <- read_csv("clade7_vs_amphidromous_fdr0.2_results.csv")
sRefSp <- "RhinogobiusFormosanus"
datFam2Gene_clade2 <- read.table("N0_modified_clade2.tsv", header=T, sep="\t", stringsAsFactors = F);
datFam2Gene_clade5 <- read.table("N0_modified_clade5.tsv", header=T, sep="\t", stringsAsFactors = F);
datFam2Gene_clade7 <- read.table("N0_modified_clade7.tsv", header=T, sep="\t", stringsAsFactors = F);

fnFam2Gene <- function(arr, dat) {
  arrGenes <- dat[dat$OG %in% arr, sRefSp];
  return(unique(str_trim(unlist(strsplit(x = arrGenes, split = ',')))));
}

Expansion_Status <- 'Expanded'
Expansion_Status <- 'Contracted'
arrFamClade2 <- clade2_significant_results$Orthogroup[clade2_significant_results$Expansion_Status == Expansion_Status]
arrFamClade5 <- clade5_significant_results$Orthogroup[clade5_significant_results$Expansion_Status == Expansion_Status]
arrFamClade7 <- clade7_significant_results$Orthogroup[clade7_significant_results$Expansion_Status == Expansion_Status]
arrFamClade2Genes <- fnFam2Gene(arrFamClade2, datFam2Gene_clade2);
arrFamClade5Genes <- fnFam2Gene(arrFamClade5, datFam2Gene_clade5);
arrFamClade7Genes <- fnFam2Gene(arrFamClade7, datFam2Gene_clade7);
arrControlGenes_clade2 <- fnFam2Gene(datFam2Gene_clade2$OG, datFam2Gene_clade2)
arrControlGenes_clade5 <- fnFam2Gene(datFam2Gene_clade5$OG, datFam2Gene_clade5)
arrControlGenes_clade7 <- fnFam2Gene(datFam2Gene_clade7$OG, datFam2Gene_clade7)

###enrichment test:
go <- enrichGO(gene = arrFamClade7Genes,
               keyType="GID",
               OrgDb = org.Rformosanus.clade7.eg.db,
               ont = "ALL",
               universe = arrControlGenes_clade7 ,
               pAdjustMethod = "BH",
               pvalueCutoff = 1,
               qvalueCutoff =1)
write.table(go, paste0("Clade7", Expansion_Status, ".interprotGOwithEnsembl.fdr0.2.csv"), sep = ",", col.names = T, row.names = F, quote = F)


Contracted_interprotGOwithEnsembl_fdr0_05_overlapped <- merge(Clade2Contracted_interprotGOwithEnsembl_fdr0_05, Clade7Contracted_interprotGOwithEnsembl_fdr0_05, by = "Description")
Expanded_interprotGOwithEnsembl_fdr0_05_overlapped <- merge(Clade2Expanded_interprotGOwithEnsembl_fdr0_05, Clade7Expanded_interprotGOwithEnsembl_fdr0_05, by = "Description")
Contracted_interprotGOwithEnsembl_fdr0_2_overlapped <- merge(Clade2Contracted_interprotGOwithEnsembl_fdr0_2, Clade7Contracted_interprotGOwithEnsembl_fdr0_2, by = "Description")
Expanded_interprotGOwithEnsembl_fdr0_2_overlapped <- merge(Clade2Expanded_interprotGOwithEnsembl_fdr0_2, Clade7Expanded_interprotGOwithEnsembl_fdr0_2, by = "Description")
write.table(Contracted_interprotGOwithEnsembl_fdr0_05_overlapped, paste0("Contracted_interprotGOwithEnsembl_fdr0_05_overlapped.csv"), sep = ",", col.names = T, row.names = F, quote = F)
write.table(Contracted_interprotGOwithEnsembl_fdr0_2_overlapped, paste0("Contracted_interprotGOwithEnsembl_fdr0_2_overlapped.csv"), sep = ",", col.names = T, row.names = F, quote = F)
write.table(Expanded_interprotGOwithEnsembl_fdr0_05_overlapped, paste0("Expanded_interprotGOwithEnsembl_fdr0_05_overlapped.csv"), sep = ",", col.names = T, row.names = F, quote = F)
write.table(Expanded_interprotGOwithEnsembl_fdr0_2_overlapped, paste0("Expanded_interprotGOwithEnsembl_fdr0_2_overlapped.csv"), sep = ",", col.names = T, row.names = F, quote = F)


# 加载必要的包
library(tidyverse)

# 读取 CSV 结果
clade2_significant_results <- read_csv("clade2_vs_amphidromous_fdr0.2_results.csv")
clade7_significant_results <- read_csv("clade7_vs_amphidromous_fdr0.2_results.csv")

# 读取 OG 到基因的映射文件（每行: OG, 参考物种基因ID列表用逗号分隔）
datFam2Gene_clade2 <- read.table("N0_modified_clade2.tsv", header=T, sep="\t", stringsAsFactors = F)
datFam2Gene_clade7 <- read.table("N0_modified_clade7.tsv", header=T, sep="\t", stringsAsFactors = F)

# 参考物种名称
sRefSp <- "RhinogobiusFormosanus"

# 函数：从 OG 列表获取该物种的所有基因ID（去重）
fnFam2Gene <- function(arr, dat) {
  # arr: OG 名称向量
  # dat: OG->基因映射数据框，其中有一列是参考物种的基因ID（逗号分隔）
  arrGenes <- dat[dat$OG %in% arr, sRefSp]
  # 拼接所有基因ID，然后按逗号分割，去除空白，去重
  all_genes <- unlist(strsplit(paste(arrGenes, collapse = ","), ","))
  all_genes <- str_trim(all_genes)
  all_genes <- all_genes[all_genes != ""]
  return(unique(all_genes))
}

fnAddGeneSymbols_flex <- function(input, 
                                  input_gene_col = "GeneID",   # 输入数据框中的基因ID列名
                                  ref_gene_col = "Gene",       # datSymbol 中的基因ID列名
                                  ref_symbol_col = "Symbol",   # datSymbol 中的符号列名
                                  output_csv = NULL) {
  if (!exists("datSymbol")) stop("datSymbol not found")
  if (!all(c(ref_gene_col, ref_symbol_col) %in% colnames(datSymbol))) {
    stop("datSymbol missing required columns")
  }
  
  gene2symbol <- setNames(datSymbol[[ref_symbol_col]], datSymbol[[ref_gene_col]])
  
  map_gene <- function(gid) {
    sym <- gene2symbol[as.character(gid)]
    if (is.na(sym)) return(as.character(gid)) else return(sym)
  }
  
  if (is.character(input)) {
    symbols <- sapply(input, map_gene, USE.NAMES = FALSE)
    res <- data.frame(GeneID = input, GeneSymbol = symbols, stringsAsFactors = FALSE)
    if (!is.null(output_csv)) write.csv(res, output_csv, row.names = FALSE)
    return(res)
  }
  
  if (is.data.frame(input)) {
    if (!input_gene_col %in% colnames(input)) {
      stop("Input data frame lacks column: ", input_gene_col)
    }
    input$genesymbols <- sapply(input[[input_gene_col]], map_gene)
    return(input)
  }
  
  stop("input must be character vector or data frame")
}
# ================= 选择扩张或收缩 =================
Expansion_Status <- "Expanded"   # 可选 "Contracted"
Expansion_Status <- "Contracted"

# 获取显著 OG
arrFamClade2 <- clade2_significant_results$Orthogroup[clade2_significant_results$Expansion_Status == Expansion_Status]
arrFamClade7 <- clade7_significant_results$Orthogroup[clade7_significant_results$Expansion_Status == Expansion_Status]

# 取交集
overlap_ogs <- intersect(arrFamClade2, arrFamClade7)

# 从映射表获取重叠OG对应的所有基因ID（每个OG内可能有多个基因，用逗号分隔）
# 方法：先按OG提取基因字符串，再拆分成单个基因，保持OG与基因的对应关系
gene_list <- list()
for (og in overlap_ogs) {
  # 从 clade2 的映射表获取该OG的基因字符串（两个表内容应一致，任选一个）
  gene_str <- datFam2Gene_clade2[datFam2Gene_clade2$OG == og, sRefSp]
  if (length(gene_str) == 0) next
  # 分割基因ID
  genes <- unlist(strsplit(gene_str, ","))
  genes <- str_trim(genes)
  genes <- genes[genes != ""]
  if (length(genes) > 0) {
    gene_list[[og]] <- data.frame(OG = og, GeneID = genes, stringsAsFactors = FALSE)
  }
}
# 合并所有基因，每个基因一行
df_genes <- bind_rows(gene_list)

colnames(datSymbol)[colnames(datSymbol) == "symbol"] <- "Symbol"

# 添加基因符号（使用改进后的函数，假设 datSymbol 已经存在）
if (exists("datSymbol") && nrow(df_genes) > 0) {
  df_genes_with_symbol <- fnAddGeneSymbols_flex(df_genes, 
                                                input_gene_col = "GeneID",
                                                ref_gene_col = "Gene",
                                                ref_symbol_col = "Symbol")
  # 重命名列：原函数返回的数据框有基因ID列名为 geneID（小写？不，它用的是 input 的列名）
  # 注意 fnAddGeneSymbols 当输入是数据框时，会添加 genesymbols 列，原列保持不变
  # 所以最终列名为：OG, GeneID, genesymbols
} else {
  warning("datSymbol 未找到，无法添加基因符号")
  df_genes_with_symbol <- df_genes
  df_genes_with_symbol$genesymbols <- NA
}

# 按所有列去重（即完全相同的行只保留一条）
df_genes_unique <- distinct(df_genes_with_symbol)
df_genes_unique_deg <- df_genes_unique %>% filter(GeneID %in% DEGs)

# 输出结果文件
output_file <- paste0("overlap_", tolower(Expansion_Status), "_genes_clade2_7DEG.csv")
write.csv(df_genes_unique_deg, output_file, row.names = FALSE)

# 查看结果
print(head(df_genes_with_symbol))
