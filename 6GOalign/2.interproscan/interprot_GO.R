library(clusterProfiler)#进行GO富集和KEGG富集
library(dplyr) #进行数据转换
library(ggplot2)#绘图
library(stringr)
library(enrichplot)
library(tidyverse)
library(GO.db)
library(AnnotationForge) #加载
options(stringsAsFactors = F)

setwd("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/GO")

# only interprotscan
pather <- read.csv("./interprotscan/Rform.prot_no_stop.fasta.tsv", sep='\t', header=T, stringsAsFactors=FALSE)

go_long <- pather %>%
  separate_rows(GO, sep = "\\|") %>%
  rename(GeneID = Protein_id, GO = GO)
gene2go <- go_long %>%
  distinct(GeneID, GO) %>%
  rename(GID = GeneID, GOID = GO)
gene2go <- gene2go %>%
  filter(GOID != "-")
gene2go$GOID <- gsub("\\s*\\(.*\\)$", "", gene2go$GOID)
duplicated_rows <- duplicated(gene2go)
if (any(duplicated_rows)) {
  cat("Found duplicated rows in gene2go:\n")
  print(gene2go[duplicated_rows, ])
}
gene2go <- gene2go[!duplicated(gene2go), ]
colnames(gene2go) <- c("GID", "GO")
gene2go$EVIDENCE <- "PANTHER"
# 创建数据库
AnnotationForge::makeOrgPackage(
  gene_data = gene2go,
  go = gene2go,
  version = "1.0",
  maintainer = 'Wanying Xu <xuwy57@mail2.sysu.edu.cn>',
  author = "Wanying Xu",
  outputDir = ".",  # 输出目录
  tax_id = "508009",
  genus = "Rhinogobius",
  species = "formosanus",
  goTable = "go"
)
install.packages('org.Rformosanus.eg.db',repos = NULL, type="source") #安装包
library(org.Rformosanus.eg.db) #加载包
remove.packages("org.Rformosanus.eg.db" , lib = file.path("/public3/group_crf/home/b20xuwy57/R/x86_64-pc-linux-gnu-library"))




# with danio_rerio oryzias_latipes
sSp <- "RhinogobiusFormosanus";
sTERM2GO <- paste0(sSp, "_term2go.tsv");
datAllGO <- read.table(sTERM2GO, sep="\t", header=T, stringsAsFactors = T);
colnames(datAllGO) <- c("GO", "GID")
pather <- read.csv("./interprotscan/Rform.prot_no_stop.fasta.tsv", sep='\t', header=T, stringsAsFactors=FALSE)

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
# 创建数据库
AnnotationForge::makeOrgPackage(
  gene_data = gene2go,
  go = gene2go,
  version = "2.0",
  maintainer = 'Wanying Xu <xuwy57@mail2.sysu.edu.cn>',
  author = "Wanying Xu",
  outputDir = ".",  # 输出目录
  tax_id = "508009",
  genus = "Rhinogobius",
  species = "formosanus.v2",
  goTable = "go"
)
install.packages('org.Rformosanus.v2.eg.db',repos = NULL, type="source") #安装包
library(org.Rformosanus.v2.eg.db) #加载包
remove.packages("org.Rformosanus.v2.eg.db" , lib = file.path("/public3/group_crf/home/b20xuwy57/R/x86_64-pc-linux-gnu-library"))
