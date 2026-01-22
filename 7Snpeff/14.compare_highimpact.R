#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# Helper: find a column name in a data.frame ignoring case
find_col_ci <- function(df, name) {
  cols <- colnames(df)
  idx <- which(tolower(cols) == tolower(name))
  if (length(idx) == 0) return(NA_character_)
  cols[idx[1]]
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("用法: Rscript compare_groups.R <count_file> <species_file> <output_prefix>")
}
count_file <- args[1]
species_file <- args[2]
output_prefix <- args[3]

# 目标 clade 列表
target_clades <- c("clade2", "clade5", "clade7")
# ✅ 定义每个 clade 的最小物种数要求
min_species_df <- data.frame(
  clade2 = 8,
  clade5 = 3,
  clade7 = 8
)

# 读取计数表（基因 x 物种），要求第一列是 gene_id 或脚本会把第一列当作 gene_id
df <- read_tsv(count_file, col_types = cols())
if (!("gene_id" %in% colnames(df))) {
  # 如果头部不是 gene_id，就把第一列名改为 gene_id
  colnames(df)[1] <- "gene_id"
}
message("计数表: ", nrow(df), " 个基因，", ncol(df)-1, " 个物种列 (第一列 gene_id)")

# 读取物种信息（空格分隔/可变列名）
species_info <- read.table(species_file, header = TRUE, sep = "", stringsAsFactors = FALSE)
message("读取物种信息: ", nrow(species_info), " 条记录; 列名: ", paste(colnames(species_info), collapse = ", "))

# 找到必需列名（大小写不敏感）
species_col <- find_col_ci(species_info, "species")
clade_col   <- find_col_ci(species_info, "clade")
lifestyle_col <- find_col_ci(species_info, "lifestyle")

if (is.na(species_col)) stop("物种信息文件必须包含 'species' 列 (不区分大小写)")
if (is.na(clade_col))   stop("物种信息文件必须包含 'clade' 列 (不区分大小写)")
if (is.na(lifestyle_col)) stop("物种信息文件必须包含 'lifestyle' 列 (不区分大小写)")

# 获取 amphidromous 物种（lifestyle 包含 amphidromous，忽略大小写）
amphi_mask <- grepl("amphidromous", species_info[[lifestyle_col]], ignore.case = TRUE)
amphidromous_species_all <- unique(species_info[[species_col]][amphi_mask])
message("找到 amphidromous 物种 (总): ", length(amphidromous_species_all))

# 物种列在计数表中的实际列名（交集）
all_count_species <- colnames(df)[-1]  # 去掉 gene_id
if (length(all_count_species) == 0) stop("计数表中没有物种列 (只有 gene_id?)")

# 对每个 target_clade 做比较
for (clade in target_clades) {
  # 按 clade 取物种
  clade_mask <- tolower(as.character(species_info[[clade_col]])) == tolower(clade)
  clade_species_all <- unique(species_info[[species_col]][clade_mask])
  message("=== 处理: ", clade, " ===")
  message("物种表中标记为 ", clade, " 的物种数 (原始): ", length(clade_species_all))
  message("amphidromous (原始): ", length(amphidromous_species_all))
  
  clade_species_in_df <- intersect(clade_species_all, all_count_species)
  amphidromous_in_df   <- intersect(amphidromous_species_all, all_count_species)
  
  message("clade 在计数表中的物种数 (使用): ", length(clade_species_in_df))
  message("amphidromous 在计数表中的物种数 (使用): ", length(amphidromous_in_df))
  
  min_required <- as.numeric(min_species_df[[clade]])
  if (length(clade_species_in_df) < min_required) {
    message("⚠️ 跳过 ", clade, "：物种数(", length(clade_species_in_df),
            ") 少于设定阈值(", min_required, ")")
    next
  }
  if (length(amphidromous_in_df) == 0) {
    message("⚠️ 无 amphidromous 物种匹配，跳过 ", clade)
    next
  }

  # ✅ 非NA样本数过滤
  valid_mask <- apply(df, 1, function(x) {
    clade_vals <- suppressWarnings(as.numeric(x[clade_species_in_df]))
    amphi_vals <- suppressWarnings(as.numeric(x[amphidromous_in_df]))
    sum(!is.na(clade_vals)) >= min_required && sum(!is.na(amphi_vals)) >= 1
  })
  df_valid <- df[valid_mask, , drop = FALSE]
  message("通过非NA过滤的基因数: ", nrow(df_valid), " / ", nrow(df))

  if (nrow(df_valid) == 0) {
    message("⚠️ 没有基因符合非NA样本数要求，跳过 ", clade)
    next
  }

  # Wilcoxon 检验
  results <- vector("list", nrow(df_valid))
  for (i in seq_len(nrow(df_valid))) {
    gene_id <- df_valid$gene_id[i]
    clade_vals <- as.numeric(unlist(df_valid[i, clade_species_in_df, drop = TRUE]))
    amphi_vals <- as.numeric(unlist(df_valid[i, amphidromous_in_df, drop = TRUE]))
    clade_vals <- clade_vals[!is.na(clade_vals)]
    amphi_vals <- amphi_vals[!is.na(amphi_vals)]

    clade_mean <- if (length(clade_vals) == 0) NA_real_ else mean(clade_vals)
    amphi_mean <- if (length(amphi_vals) == 0) NA_real_ else mean(amphi_vals)
    
    pval <- 1.0
    if (length(clade_vals) > 0 && length(amphi_vals) > 0) {
      if (length(unique(c(clade_vals, amphi_vals))) <= 1) {
        pval <- 1.0
      } else {
        test <- tryCatch(wilcox.test(clade_vals, amphi_vals, alternative = "greater"), error = function(e) NULL)
        if (!is.null(test)) pval <- test$p.value
      }
    }
    
    results[[i]] <- data.frame(
      gene_id = gene_id,
      clade_mean = clade_mean,
      amphidromous_mean = amphi_mean,
      p_value = pval,
      stringsAsFactors = FALSE
    )
  }

  res_df <- bind_rows(results)
  res_df$fdr <- p.adjust(res_df$p_value, method = "fdr")
  res_df$significant <- ifelse(is.na(res_df$fdr), FALSE, res_df$fdr < 0.05)

  # ✅ 确保列顺序，gene_id 作为第一列
  res_df <- res_df[, c("gene_id", "clade_mean", "amphidromous_mean", "p_value", "fdr", "significant")]
  
  # ✅ 写出 summary + 表格
  out_file <- paste0(output_prefix, "_", clade, "_vs_amphidromous.txt")
  summary_text <- paste0(
    "# Summary for ", clade, "\n",
    "# Minimum species required (per gene): ", min_required, "\n",
    "# Clade species in count table: ", length(clade_species_in_df), "\n",
    "# Amphidromous species in count table: ", length(amphidromous_in_df), "\n",
    "# Total genes (before NA filter): ", nrow(df), "\n",
    "# Genes passed non-NA filter: ", nrow(df_valid), "\n",
    "# Genes tested (Wilcoxon): ", nrow(res_df), "\n",
    "# Significant (FDR<0.05): ", sum(res_df$significant, na.rm = TRUE), "\n",
    "#----------------------------------------------\n"
  )
  
  # 先写入注释信息
  writeLines(summary_text, out_file)
  
  # 再写入数据（包含表头）
  write_tsv(res_df, out_file, append = TRUE, col_names = TRUE)
  
  message("已写出: ", out_file, " (", nrow(res_df), " 行 + summary)")
  message("输出表头: ", paste(colnames(res_df), collapse = ", "))
}