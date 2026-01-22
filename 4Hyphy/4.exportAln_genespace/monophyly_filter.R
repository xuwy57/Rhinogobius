library(ape)
library(phytools)

# 单系群定义
clades <- list(
  Clade1 = c("HN3_R-leavelli_Hainanzhanzhoudachengzhen", "Z15_R-leavelli_Zhejiangchangshanxian", 
             "Z18_R-cliffordpopeiR-leavelli_Zhejiangchangshanxia", "Z17_R-cliffordpopei_Zhejiangchangshanxian", 
             "HN4_R-leavelli_Hainanbaotingshilingzhen"),
  Clade2 = c("Z5_R-niger_Zhejianglongquanshizhulongzhen", "Z2_R-aporus_Zhejiangpananxianshuangxizhen", 
             "A8_R-cf-wuyiensis_Anhuiqimenlikouzhen", "F6_R-cf-wuyiensis_Fujianwuyishan", "J2_R-sp_Jiangxiwuyuan", 
             "Z13_R-sp_Zhejiangkaihuaxiansuzhuangzhen", "Z1_R-sp_Zhejianglongquanshizhulongzhen", 
             "J1_R-sp_Jiangxipingxiangshianyuanqu", "A1_R-sp_Anhuidongzhixiannixizhen", 
             "A3_R-sp_Anhuiqimendatanxiang", "A6_R-lentiginis_Anhuishitai", 
             "Z16_R-cf-lentiginis_Zhejiangchangshanxian", "Z14_R-lentiginis_Zhejiangyuyaoliangnongzhen", 
             "Z19_R-aporus_Zhejiangwuyixianchengqu", "Z3_R-wuyiensis_Zhejianglongquanshijinxizhen", 
             "Z20_R-aporus_Zhejiangwuyixianxilianxiang"),
  Clade3 = c("*Clade1", "*Clade2"),
  Clade4 = c("*Clade3", "Rhinogobius_formosanus", "L1_R-cf-brunneus_Liaoningkuandian"),
  Clade5 = c("S1_R-szechuanensis_Sichuan", "HB1_R-houheensis_Hubeixuandongxian", 
             "A7_R-sp_Anhuishitai", "Z12_R-davidi_Zhejiangchangshanxian", 
             "Z7_R-davidi_Zhejianganjixiantianhuangpingzhen", "Z11_R-sp_Zhejiangchangshanxian", 
             "Z21_R-sp_Zhejiangchangshanxiansianzhen"),
  Clade6 = c("*Clade4", "*Clade5", "F1_R-cf-rubromaculatus_Fujianputianshi"),
  Clade7 = c("HN2_R-wanchuangensis_Hainanqionghaishihuishanzhen", "HN5_R-linshuiensis_Hainanbaotingshilingzhen", 
             "HN6_R-sp_Hainanhaikou", "GX2_R-sp_Guangxifangchenggangmaluzhen", 
             "GX1_R-sp_Guangxinanning", "HN1_R-nandujiangensis_Hainanhaikoushijiuzhouzhen", 
             "GD1_R-duospilus_Guangdonglongmen", "GX3_R-sp_Guangxifangchenggangmaluzhen", 
             "F7a_R-sp_Fujianyunxiao", "F7b_R-sp_Fujianyunxiao", "F2_R-xianshuiensis_Fujianputianshi", 
             "Z4_R-cf-changtingensis_Zhejianglongquanshizhulongz", "F4_R-reticulatus_Fujianfuzhou", 
             "F3_R-reticulatus_Fujianpuchengshifulingzhen", "H1_R-sp_Hunanlianyuanshimaotangzhen", 
             "H2_R-sp_Hunanloudilouxingqu", "Z9_R-sp_Zhejiangkaihuaxiansuzhuangzhen", 
             "Z8_R-immaculatus_Zhejiangtongluxianlucixiang", "Z10_R-sp_Zhejianglongquanshizhulongzhen"),
  Clade8 = c("*Clade6", "*Clade7"),
  Clade9 = c("*Clade8", "Z6_R-giurinus_Zhejianglongquanshizhulongzhen"),
  root = c("*Clade9", "Tridentiger_bifasciatus", "Tridentiger_radiatus", 
           "Tridentiger_barbatus", "Boleophthalmus_pectinirostris", 
           "Odontamblyopus_rebecca", "Glossogobius_giuris", "Gobiopsis_macrostoma")
)

# 展开单系群定义
expand_clade <- function(clade_name, clades) {
  members <- clades[[clade_name]]
  expanded <- unlist(lapply(members, function(x) {
    if (startsWith(x, "*")) {
      expand_clade(substring(x, 2), clades)
    } else {
      x
    }
  }))
  return(expanded)
}

#getMRCA对于无根树找最近共同祖先有歧义，需要先定根
# 定义外群列表，按优先级顺序
outgroup_list <- list(
  c("Boleophthalmus_pectinirostris", "Odontamblyopus_rebecca"),
  c("Glossogobius_giuris", "Gobiopsis_macrostoma"),
  c("Tridentiger_bifasciatus", "Tridentiger_radiatus", "Tridentiger_barbatus"),
  c("Z6_R-giurinus_Zhejianglongquanshizhulongzhen")
)

# 检查并扩展外群集合
expand_outgroup <- function(tree, outgroup_tips) {
  if (!is.monophyletic(tree, outgroup_tips)) {
    cat("The provided outgroup set is not monophyletic. Expanding to the smallest monophyletic clade.\n")
    mrca_node <- getMRCA(tree, outgroup_tips)  # 获取最近共同祖先节点
    expanded_outgroup <- extract.clade(tree, mrca_node)$tip.label  # 提取后代物种
    cat("Expanded outgroup set:\n", paste(expanded_outgroup, collapse = ", "), "\n")
    return(expanded_outgroup)
  } else {
    cat("The provided outgroup set is monophyletic. No expansion needed.\n")
    return(outgroup_tips)
  }
}

# 尝试依次用外群定根
set_root <- function(tree, outgroup_list, auto_expand = TRUE) {
  for (outgroup in outgroup_list) {
    # 检查外群中的物种是否存在于树中
    valid_outgroup <- outgroup[outgroup %in% tree$tip.label]
    
    if (length(valid_outgroup) > 0) {
      cat("Using outgroup:", paste(valid_outgroup, collapse = ", "), "\n")
      
      # 如果启用自动扩展，检查并扩展外群集合
      if (auto_expand) {
        valid_outgroup <- expand_outgroup(tree, valid_outgroup)
      }
      
      # 尝试用外群定根
      return(root(tree, outgroup = valid_outgroup, resolve.root = TRUE))
    }
  }
  
  return(NULL)  # 无法定根时返回 NULL
}



# 获取最近共同祖先及最小单系群的实际物种
get_actual_clade <- function(tree, tips) {
  # 检查物种是否在树中
  valid_tips <- tips[tips %in% tree$tip.label]
  
  # 如果没有有效的物种，发出警告并返回空
  if (length(valid_tips) <= 1) {
    warning("None of the provided tips are found in the tree.")
    return(valid_tips)
  }
  mrca_node <- tryCatch({
    getMRCA(tree, valid_tips)
  }, error = function(e) {
    warning("Failed to find MRCA for tips: ", paste(valid_tips, collapse = ", "))
    return(NA)
  })
  # 如果 MRCA 无效，则返回空
  if (is.na(mrca_node)) return(NA)
  
  actual_clade <- extract.clade(tree, mrca_node)$tip.label  # 提取最小单系群的物种
  return(actual_clade)
}



#遍历所有单系群，找到每个单系群实际的最小单系群。
#比较定义单系群与实际单系群，记录多余物种及其数量。

#按多余物种的数量从少到多排序，优先修复多余物种最少的单系群。

#在每次移除多余物种后，重新计算所有单系群的多余物种和数量。
#重复以上过程，直到所有单系群符合定义。


# 计算多余物种的列表
compute_extra_tips <- function(tree, clades) {
  extra_tips_list <- list()

  for (clade_name in names(clades)) {
    tips <- expand_clade(clade_name, clades)  # 定义的单系群物种
    actual_clade <- get_actual_clade(tree, tips)  # 实际单系群物种
    
    # 如果定义的单系群在树中不存在，跳过
    if (length(actual_clade) == 0) {
      cat("Skipping clade:", clade_name, "- None of its tips are found in the tree.\n")
      next
    }
    
    extra_tips <- setdiff(actual_clade, tips)  # 找出多余物种
    extra_tips_list[[clade_name]] <- extra_tips
    
  }
  return(extra_tips_list)
}


# 主循环：修复单系群
fix_clades <- function(tree, clades) {
  while (TRUE) {
    # 计算每个单系群多余的物种
    extra_tips_list <- compute_extra_tips(tree, clades)
    
    # 统计每个单系群多余物种的数量
    extra_tips_count <- sapply(extra_tips_list, length)
    
    # 找到多余物种最少的单系群
    valid_counts <- extra_tips_count[extra_tips_count > 0]  # 排除 0 的单系群
    if (length(valid_counts) == 0) {
      cat("Done\n")
      break
    }  # 如果没有多余物种直接退出循环
    
    target_clade <- names(valid_counts)[which.min(valid_counts)]
    extra_tips <- extra_tips_list[[target_clade]]
    
    # 移除该单系群的一个多余物种
    cat("Removing tip:", extra_tips[1], "from clade:", target_clade, "\n")
    tree <- drop.tip(tree, extra_tips[1])
  }
  return(tree)
}


# 检查特定单系群是否有至少 3 个物种
check_clade_members <- function(tree, clades_to_check) {
  for (clade_name in clades_to_check) {
    tips <- clades[[clade_name]]
    valid_tips <- tips[tips %in% tree$tip.label]
    if (length(valid_tips) >= 3) {
      return(TRUE)  # 至少一个单系群满足条件
    }
  }
  return(FALSE)  # 没有单系群满足条件
}


# 指定输入目录
input_dir <- "genetrees_improved"

# 获取所有文件
file_list <- list.files(path = input_dir, pattern = "^ret_part.*\\.txt$", full.names = TRUE)

final_result <- data.frame(ID = character(), RemovedTips = character(), FileName = character(), stringsAsFactors = FALSE)

# 循环处理每个文件
for (file in file_list) {
  cat("Processing file:", file, "\n")
  # 读取文件并检查列数
  data <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
  
  # 如果列数不足 6 则跳过该文件
  if (ncol(data) < 6) {
    cat("Skipping file:", file, "- Insufficient columns.\n")
    next
  }
  colnames(data) <- c("ID", "Value2", "Status", "Value4", "Value5", "Tree")
  
  filtered_data <- data[data$Status == "rejected_monophyly", ]
  
  for (i in 1:nrow(filtered_data)) {
    row <- filtered_data[i, ]
    tree <- read.tree(text = row$Tree)
    
    # 定根
    tree <- set_root(tree, outgroup_list)
    if (is.null(tree)) {
      final_result <- rbind(final_result, data.frame(
        ID = row$ID,
        RemovedTips = "N/A",
        FileName = file,
        Status = "NO_OUTGROUP"
      ))
      next
    }
    
    # 修复单系群
    original_tips <- tree$tip.label
    tree <- fix_clades(tree, clades)
    removed_tips <- setdiff(original_tips, tree$tip.label)
    
    # 检查 Clade2、Clade5 或 Clade7 是否有至少 3 个物种
    status <- ifelse(check_clade_members(tree, c("Clade2", "Clade5", "Clade7")), "PASS", "FAIL")
    
    # 添加结果
    final_result <- rbind(final_result, data.frame(
      ID = row$ID, 
      RemovedTips = paste(removed_tips, collapse = ","),
      FileName = file,
      Status = status
    ))
  }
}

# 保存结果
write.table(final_result, file = "final_output_table.txt", sep = "\t", row.names = FALSE, quote = FALSE)
cat("All files processed. Results saved to 'final_output_table.txt'.\n")