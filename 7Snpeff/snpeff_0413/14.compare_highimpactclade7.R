#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

# Helper functions
find_col_ci <- function(df, name) {
  cols <- colnames(df)
  idx <- which(tolower(cols) == tolower(name))
  if (length(idx) == 0) return(NA_character_)
  cols[idx[1]]
}

safe_wilcox <- function(x, y, alternative = "greater") {
  if (length(x) == 0 || length(y) == 0) return(NA_real_)
  if (length(unique(c(x, y))) <= 1) return(1.0)
 
  test <- tryCatch(
    wilcox.test(x, y, alternative = alternative),
    error = function(e) NULL
  )
  if (is.null(test)) return(NA_real_)
  test$p.value
}

safe_mean <- function(x) {
  if (length(x) == 0) return(NA_real_)
  mean(x, na.rm = TRUE)
}

passes_na_filter <- function(gene_row, clade_species, amphidromous_species, min_clade) {
  clade_vals <- as.numeric(gene_row[clade_species])
  amphi_vals <- as.numeric(gene_row[amphidromous_species])
 
  clade_non_na <- sum(!is.na(clade_vals))
  amphi_non_na <- sum(!is.na(amphi_vals))
 
  clade_non_na >= min_clade && amphi_non_na >= 1
}

# ====================== 主程序 ======================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("用法: Rscript compare_groups.R <count_file> <species_file> <output_prefix>")
}

count_file <- args[1]
species_file <- args[2]
output_prefix <- args[3]

# 配置参数
target_clades <- c("clade2", "clade5", "clade7")
min_species_df <- list(clade2 = 8, clade5 = 3, clade7 = 8)
fdr_threshold <- 0.2
pvalue_threshold <- 0.05

message("\n=== 读取数据 ===")
df <- read_tsv(count_file, col_types = cols())

if (!("gene_id" %in% colnames(df))) {
  colnames(df)[1] <- "gene_id"
  message("已将第一列重命名为 'gene_id'")
}

species_info <- read.table(species_file, header = TRUE, sep = "",
                          stringsAsFactors = FALSE, comment.char = "")

# 找到关键列
species_col <- find_col_ci(species_info, "species")
clade_col   <- find_col_ci(species_info, "clade")
lifestyle_col <- find_col_ci(species_info, "lifestyle")

amphi_mask <- grepl("amphidromous", species_info[[lifestyle_col]], ignore.case = TRUE)
amphidromous_species_all <- unique(species_info[[species_col]][amphi_mask])
all_count_species <- colnames(df)[-1]

# ====================== 主循环 ======================
for (clade in target_clades) {
  message("\n", paste(rep("=", 60), collapse = ""))
  message("处理: ", toupper(clade))
 
  clade_mask <- tolower(as.character(species_info[[clade_col]])) == tolower(clade)
  clade_species_all <- unique(species_info[[species_col]][clade_mask])
 
  clade_species <- intersect(clade_species_all, all_count_species)
  amphidromous_species <- intersect(amphidromous_species_all, all_count_species)
 
  min_required <- min_species_df[[clade]]
 
  if (length(clade_species) < min_required) {
    message("⚠️ 跳过: ", clade, " 物种数不足 (", length(clade_species), " < ", min_required, ")")
    next
  }
 
  if (length(amphidromous_species) == 0) {
    message("⚠️ 跳过: 无 amphidromous 物种")
    next
  }
 
  # 过滤基因
  valid_mask <- apply(df, 1, function(row) {
    passes_na_filter(row, clade_species, amphidromous_species, min_required)
  })
 
  df_valid <- df[valid_mask, , drop = FALSE]
  message("通过过滤的基因: ", nrow(df_valid), " / ", nrow(df))
 
  if (nrow(df_valid) == 0) next
 
  # 执行 Wilcoxon 检验
  results <- vector("list", nrow(df_valid))
 
  for (i in seq_len(nrow(df_valid))) {
    gene_id <- df_valid$gene_id[i]
   
    clade_vals <- as.numeric(df_valid[i, clade_species, drop = TRUE])
    amphi_vals <- as.numeric(df_valid[i, amphidromous_species, drop = TRUE])
   
    clade_vals <- clade_vals[!is.na(clade_vals)]
    amphi_vals <- amphi_vals[!is.na(amphi_vals)]
   
    clade_mean <- safe_mean(clade_vals)
    amphi_mean <- safe_mean(amphi_vals)
   
    p_value <- safe_wilcox(clade_vals, amphi_vals, alternative = "greater")
   
    results[[i]] <- data.frame(
      gene_id = gene_id,
      clade_mean = clade_mean,
      amphidromous_mean = amphi_mean,
      p_value = p_value,
      stringsAsFactors = FALSE
    )
  }
 
  res_df <- bind_rows(results)
  res_df$fdr <- p.adjust(res_df$p_value, method = "fdr")
  res_df$significant_p   <- !is.na(res_df$p_value) & res_df$p_value < pvalue_threshold
  res_df$significant_fdr <- !is.na(res_df$fdr) & res_df$fdr < fdr_threshold
 
  res_df <- res_df %>%
    arrange(p_value) %>%
    select(gene_id, clade_mean, amphidromous_mean, p_value, fdr,
           significant_p, significant_fdr)
 
  out_file <- paste0(output_prefix, "_", clade, "_vs_amphidromous.txt")
 
  # ==================== 准备摘要信息 ====================
  summary_lines <- c(
    "# ================================================",
    paste0("# Comparison: ", clade, " vs Amphidromous"),
    paste0("# Date: ", Sys.time()),
    paste0("# Minimum clade species required: ", min_required),
    paste0("# FDR threshold: ", fdr_threshold),
    paste0("# P-value threshold: ", pvalue_threshold),
    "#",
    "# Sample counts:",
    paste0("# ", clade, " species in count table: ", length(clade_species)),
    paste0("# Amphidromous species in count table: ", length(amphidromous_species)),
    "#",
    "# Gene statistics:",
    paste0("# Total genes (before filtering): ", nrow(df)),
    paste0("# Genes passed NA filter: ", nrow(df_valid)),
    paste0("# Genes tested (Wilcoxon): ", nrow(res_df)),
    paste0("# Significant (p < 0.05): ", sum(res_df$significant_p, na.rm = TRUE)),
    paste0("# Significant (FDR < 0.2): ", sum(res_df$significant_fdr, na.rm = TRUE)),
    "# ================================================",
    "#",
    "# Column descriptions:",
    "# gene_id             : Gene identifier",
    paste0("# clade_mean          : Mean expression in ", clade, " species"),
    "# amphidromous_mean   : Mean expression in amphidromous species",
    "# p_value             : Wilcoxon rank-sum test p-value (clade > amphidromous)",
    "# fdr                 : FDR adjusted p-value",
    "# significant_p       : Significant at p < 0.05",
    "# significant_fdr     : Significant at FDR < 0.2",
    "# ================================================",
    ""
  )
 
  # ==================== 写入结果 ====================
  writeLines(summary_lines, out_file)
 
  write_delim(res_df,
              file = out_file,
              delim = "\t",
              col_names = TRUE,
              quote = "none",
              append = TRUE,
              na = "NA")
 
  message("✓ 结果已保存: ", out_file)
  message("  总测试基因数: ", nrow(res_df))
  message("  p < 0.05 显著: ", sum(res_df$significant_p, na.rm = TRUE))
  message("  FDR < 0.2 显著: ", sum(res_df$significant_fdr, na.rm = TRUE))
}

message("\n所有分析完成！")