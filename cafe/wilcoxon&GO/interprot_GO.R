library(clusterProfiler)#进行GO富集和KEGG富集
library(dplyr) #进行数据转换
library(ggplot2)#绘图
library(stringr)
library(enrichplot)
library(tidyverse)
library(GO.db)
library(AnnotationForge) #加载
options(stringsAsFactors = F)

# with danio_rerio oryzias_latipes
sSp <- "RhinogobiusFormosanus";
sTERM2GO <- paste0(sSp, "_term2go.tsv");
datAllGO <- read.table(sTERM2GO, sep="\t", header=T, stringsAsFactors = T);
colnames(datAllGO) <- c("GO", "GID")
pather <- read.csv("./interprotscan/Rform.prot_no_stop.fasta.tsv", sep='\t', header=T, stringsAsFactors=FALSE)

setwd("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/cafe5/family_test0125")
go_long <- pather %>%
  separate_rows(GO, sep = "\\|") %>%
  rename(GeneID = Protein_id, GO = GO)
gene2go <- go_long %>%
  distinct(GeneID, GO) %>%
  rename(GID = GeneID, GO = GO)
gene2go <- gene2go %>%
  filter(GO != "-")
gene2go$GO <- gsub("\\s*\\(.*\\)$", "", gene2go$GO)
gene2go <- rbind(gene2go, datAllGO)
duplicated_rows <- duplicated(gene2go)
if (any(duplicated_rows)) {
  cat("Found duplicated rows in gene2go:\n")
  print(gene2go[duplicated_rows, ])
}
gene2go <- gene2go[!duplicated(gene2go), ]
colnames(gene2go) <- c("GID", "GO")
gene2go$EVIDENCE <- "PANTHER"


#根据clade基因拷贝数进行修改
# 读取 N0.tsv 文件
n0_data <- read.table("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/genespace/rundir/orthofinder/Results_Nov05/Phylogenetic_Hierarchical_Orthogroups/N0.tsv", header = TRUE, sep = "\t")
# 提取 OG 和 RhinogobiusFormosanus 列
og_rhinogobius_data <- n0_data[, c("OG", "RhinogobiusFormosanus")]
# 读取 CSV 文件，假设文件包含 Orthogroup 和 Clade_Avg 列
csv_data <- read.csv("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/cafe5/family_test0122/clade7_vs_amphidromous_all_results.csv")
set.seed(123)  # 全局设置随机种子
# 预处理数据 ----
og_rhinogobius <- n0_data %>%
  group_by(OG) %>%
  summarise(RhinogobiusFormosanus = paste(RhinogobiusFormosanus, collapse = ", ")) %>%
  ungroup() %>%
  mutate(RhinogobiusFormosanus = sapply(strsplit(RhinogobiusFormosanus, ",\\s*"), function(x) {
    # 去除空字符串
    x <- x[x != ""]
    paste(x, collapse = ", ")
  }))
orthogroup_clade <- csv_data %>% select(Orthogroup, Clade_Avg)
# 生成modified_data ----
modified_list <- lapply(1:nrow(og_rhinogobius), function(i) {
  og <- og_rhinogobius$OG[i]
  genes <- strsplit(og_rhinogobius$RhinogobiusFormosanus[i], ",")[[1]] %>% trimws()
  # 如果 genes 为空，跳过此行
  if (length(genes) == 0) return(NULL)
  clade_avg <- orthogroup_clade %>% 
    filter(Orthogroup == og) %>% 
    pull(Clade_Avg)
  # 如果 Clade_Avg 为空，跳过此行
  if (length(clade_avg) == 0) return(NULL)
  desired_count <- round(clade_avg)
  # 基因数量调整
  if (length(genes) < desired_count) {
    need_add <- desired_count - length(genes)
    pseudo_gids <- paste0("pseudo_", og, "_", seq_len(need_add))
    genes <- c(genes, pseudo_gids)
  } else if (length(genes) > desired_count) {
    genes <- sample(genes, desired_count)
  }
  # 返回数据框
  # 确保genes不为空，避免生成空的data.frame
  if (length(genes) > 0) {
    return(data.frame(OG = og, GID = genes, stringsAsFactors = FALSE))  # 如果有基因，创建data.frame
  } else {
    return(NULL)  # 如果没有基因，则返回NULL
  }
}) %>% bind_rows()
# 5. 合并GO信息 ----
final_data <- modified_list %>%
  left_join(gene2go, by = "GID") %>%
  rename(GO = GO, EVIDENCE = EVIDENCE)  # 调整列名匹配实际数据
# 补充缺失的GO和EVIDENCE信息 ----
final_data <- final_data %>%
  group_by(OG) %>%
  mutate(
    # 如果GO或EVIDENCE缺失，则从同OG的其他基因复制GO和EVIDENCE
    GO = ifelse(is.na(GO), paste(GO[!is.na(GO)], collapse = ", "), GO),
    EVIDENCE = ifelse(is.na(EVIDENCE), paste(EVIDENCE[!is.na(EVIDENCE)], collapse = ", "), EVIDENCE)
  ) %>%
  ungroup()
final_data_gene2go <- final_data %>%
  select(-OG) %>%  # 删除 OG 列
  separate_rows(GO, sep = ",\\s*") %>%
  filter(GO != "") %>%
  mutate(EVIDENCE = "PANTHER") %>%  # 将 EVIDENCE 全部设置为 "PANTHER"
  distinct()
# 更新N0数据 ----
n0_modified <- n0_data %>%
  select(-RhinogobiusFormosanus) %>%
  left_join(
    final_data %>%
      group_by(OG) %>%
      summarise(RhinogobiusFormosanus = paste(GID, collapse = ", ")),
    by = "OG"
  )
# 保存结果 ----
write.table(n0_modified, "N0_modified_clade7.tsv", sep="\t", row.names=FALSE, quote=FALSE)









# 创建数据库
AnnotationForge::makeOrgPackage(
  gene_data = final_data_gene2go,
  go = final_data_gene2go,
  version = "3.0",
  maintainer = 'Wanying Xu <xuwy57@mail2.sysu.edu.cn>',
  author = "Wanying Xu",
  outputDir = ".",  # 输出目录
  tax_id = "508009",
  genus = "Rhinogobius",
  species = "formosanus.clade7",
  goTable = "go"
)
install.packages('org.Rformosanus.clade7.eg.db',repos = NULL, type="source") #安装包
library(org.Rformosanus.clade2.eg.db) #加载包
remove.packages("org.Rformosanus.v2.eg.db" , lib = file.path("/public3/group_crf/home/b20xuwy57/R/x86_64-pc-linux-gnu-library"))





orthogroup <- clade2_vs_amphidromous_fdr0_05_results$Orthogroup
orthogroup <- c(orthogroup, "OG0020698")
fig2data <- subset(data, data$Orthogroup %in% orthogroup)
fig2data <- t(fig2data)
write.table(fig2data, "/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/plot/plot2/fig2data.csv", sep=",", row.names=TRUE, quote=FALSE)
