library(ggplot2)
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(Glimma)
#library(ComplexHeatmap)
library(circlize)
library(dendextend)
#vignette("DESeq2")
library(paletteer)
library(devtools)
library(ggbiplot)
library(factoextra)
library(Rling); library(fields); library(rgl); library(MASS)

#################
setwd("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/deseq2")
count.data <- read.table("counts.txt.renamed.txt", header = T, row.names = 1L)
keep_feature <- rowSums(count.data) > 17
table(keep_feature)  
count.data <- count.data[keep_feature, ] 
expressed <- row.names(count.data)

#################
index<-rowSums(count.data)
count.data_index <- count.data[order(rowSums(count.data)),]  #order对count.data按行总和进行排序

coldata<-data.frame(name=colnames(count.data_index),
                    salinity=c(rep("Fresh",3),rep("Fresh",3),rep("Fresh",3),rep("Salty",2),rep("Salty",3),rep("Salty",3)),
                      group=c(rep("F1",3),rep("F2",3),rep("F3",3),rep("S1",2),rep("S2",3),rep("S3",3)))
coldata$salt <- factor(coldata$salinity, 
                            levels = c("Fresh", "Salty"))
coldata$group <- factor(coldata$group, 
                           levels = c("F1", "F2", "F3", "S1", "S2", "S3"))
dds<-DESeqDataSetFromMatrix(countData = count.data_index,colData = coldata,
                              design = ~ salinity)
dds<-DESeq(dds)

counts <- as.data.frame(counts(dds, normalized=TRUE))
DEGs <- read.table("../deseq2/DEGs.csv", header=T, sep=",", fill = T, quote = "")
DEGs <- DEGs[,1]
counts$sum <- rowSums(counts)
counts_sig <- counts[as.character(DEGs),]
counts_sig <- na.omit(counts_sig)
sum_sig <- counts_sig$sum
sum_sig <- sort(sum_sig, decreasing=T) 
index <- round(length(sum_sig) * 95 / 100)
min_5 <- sum_sig[index]
counts_sim <- subset(counts, min_5 <= sum)
counts_sim_id <- rownames(counts_sim)
write.csv(counts_sim,file = "sim_genes.csv", quote = F)

sim_up <- res[which(res$log2FoldChange > 0), ]
sim_up <- sim_up[rownames(sim_up) %in% counts_sim_id, ]
write.csv(sim_up,file = "sim_up.csv", quote = F)
sim_down <- res[which(res$log2FoldChange < 0), ]
sim_down <- sim_down[rownames(sim_down) %in% counts_sim_id, ]
write.csv(sim_down,file = "sim_down.csv", quote = F)

## Fresh vs Salty
res <- results(dds, contrast=c("salinity", "Salty", "Fresh"))
DEG = res[which(abs(res$log2FoldChange) > 1 & res$padj < 0.05), ]
DEGs_up = res[which(res$log2FoldChange > 1 & res$padj < 0.05), ]
DEGs_down = res[which(res$log2FoldChange < -1 & res$padj < 0.05), ]
write.csv(DEG,file = "DEGs.csv", quote = F)
write.csv(DEGs_down,file = "DEGs_down.csv", quote = F)
write.csv(DEGs_up,file = "DEGs_up.csv", quote = F)
res_sig <- subset(res, padj < 0.05)

write.csv(expressed,file = "expressed.csv", quote = F)


#### PCA
dat.pca <- PCA(t(dat) , graph = F) 
pca <- fviz_pca_ind(dat.pca,
                    title = "Principal Component Analysis",
                    legend.title = "Groups",
                    geom.ind = c("point", "text"), 
                    pointsize = 1.5,
                    labelsize = 4,
                    col.ind = group_list, # 分组上色
                    axes.linetype=NA,  # remove axeslines
                    mean.point=F#去除分组中心点
) + 
  coord_fixed(ratio = 1)+  #坐标轴的纵横比
  xlab(paste0("PC1 (",round(percentVar[1,'variance.percent'],1),"%)")) +
  ylab(paste0("PC2 (",round(percentVar[2,'variance.percent'],1),"%)"))

# 将所有样本转换为 rlog
dds_rlog <- rlog(dds, blind = FALSE)
plotPCA(dds_rlog, intgroup = "salinity", ntop = 500) +
  theme_bw() +
  geom_point(size = 5) + 
  scale_y_continuous(limits = c(-50, 50)) +
  scale_x_continuous(limits = c(-50, 50)) +
  ggtitle(label = "Principal Component Analysis (PCA)", 
          subtitle = "Top 500 most variable genes") 
text(pcaData[,1],pcaData[,2]+0.2,row.names(pcaData),cex=0.5)




# 收集40个显著基因，制作矩阵
mat <- assay(rld_mat[sigLRT_genes5$gene])[1:40, ]
# 选择您要用来注释列的列变量。
annotation_col = data.frame(
  Group = factor(colData(rld_mat)$salinity),
  row.names = colData(rld_mat)$name
)
# 指定要用来注释列的颜色。
ann_colors = list(
  Group = c(Control = "#C16D58", Low = "#F2EBE5", Moderate = "#C8D6E7", High="#375093"))
# 使用 pheatmap 功能制作热图。
pheatmap(mat = mat,
         color = colorRampPalette(c("#375093", "#F2EBE5", "#831A21"))(256),
         scale = "row",
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         fontsize = 6.5,
         cellwidth = 10,
         show_colnames = F)

#### PCA plot
# 按列变量绘制 PCA
plotPCA(rld_mat, intgroup = "salinity") +
  theme_bw() +
  geom_point(size = 1) +
  scale_y_continuous(limits = c(-100, 100)) +
  ggtitle(label = "Principal Component Analysis (PCA)",
          subtitle = "Top 500 most variable genes")




normalized <- counts(dds, normalized=TRUE)
normalized_t <- t(normalized)
data.pca <- prcomp(normalized_t)
# 查看PCA的分析结果
summary(data.pca)
# 绘制主成分的碎石图
screeplot(data.pca, npcs = 10, type = "lines")
#查看主成分的结果
pca.scores <- as.data.frame(data.pca$x)
head(pca.scores)

color_scheme <-  paletteer_c("grDevices::Cividis",n=4)
group=c(c(rep("Fresh", 9), rep("Salty", 8)))
pdf("PCA_allgenes.pdf", width=12, height=8)
fviz_pca_ind(data.pca, col.ind=group, mean.point=F, addEllipses = T,
             legend.title="salinity", ellipse.type="confidence",
             ellipse.level=0.9, palette = color_scheme) +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))#加个边框
dev.off()

pdf("PCA_allsamples2.pdf", width=12, height=8)
ggbiplot(data.pca,obs.scale = 1,var.scale = 1,
         groups = group,ellipse = T,var.axes = F,
         #labels = row.names(pca.scores.rm)
) +
  scale_color_brewer(palette = "Set1") +
  theme(legend.direction = 'horizontal',legend.position = 'top') +
  geom_vline(xintercept = c(0), linetype = 6,color = "black") +
  geom_hline(yintercept = c(0),linetype = 6,color = "black") +
  theme_bw() +
  theme(panel.grid = element_line(colour = NA))
dev.off()

normalized_a.rm <- counts(dds_salty.rm, normalized=TRUE)
normalized_a_t.rm <- t(normalized_a.rm)
data.pca.rm <- prcomp(normalized_a_t.rm)
# 查看PCA的分析结果
summary(data.pca.rm)
# 绘制主成分的碎石图
screeplot(data.pca.rm, npcs = 10, type = "lines")
#查看主成分的结果
pca.scores.rm <- as.data.frame(data.pca.rm$x)
head(pca.scores.rm)

library(factoextra)
color_scheme <-  paletteer_c("grDevices::Cividis",n=4)
color_scheme <-  paletteer_c("awtools::Cividis",n=4)
group.rm=c(c(rep("Control", 3), rep("Low", 10), rep("Moderate", 11), rep("High", 5)))
pdf("PCA_rm2samples1.pdf", width=12, height=8)
fviz_pca_ind(data.pca.rm, col.ind=group.rm, mean.point=F, addEllipses = T,
             legend.title="salty", ellipse.type="confidence",
             ellipse.level=0.9, palette = color_scheme) +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))#加个边框
dev.off()
pdf("PCA_rm2samples2.pdf", width=12, height=8)
ggbiplot(data.pca.rm,obs.scale = 1,var.scale = 1,
         groups = group.rm,ellipse = T,var.axes = F,
         #labels = row.names(pca.scores.rm)
) +
  scale_color_brewer(palette = "Set1") +
  theme(legend.direction = 'horizontal',legend.position = 'top') +
  geom_vline(xintercept = c(0), linetype = 6,color = "black") +
  geom_hline(yintercept = c(0),linetype = 6,color = "black") +
  theme_bw() +
  theme(panel.grid = element_line(colour = NA))
dev.off()



# 加载包
library(DESeq2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# 1. 读取数据 -----------------------------------------------------------------
# 假设:
# - dds: DESeqDataSet对象（已进行DESeq分析）
# - deg_results: 差异表达结果（data.frame格式）
# - gene_list: 目标基因列表（字符向量）
# 2. 提取目标基因表达矩阵 -----------------------------------------------------
# 获取标准化表达矩阵（推荐使用vst转换）
vsd <- vst(dds, blind = FALSE)
expr_matrix <- assay(vsd)
gene_list <- arrTargetGenes_clade257_DEG
# 筛选目标基因（处理基因名匹配问题）
target_genes <- gene_list[gene_list %in% rownames(expr_matrix)]
if (length(target_genes) == 0) stop("无匹配基因！检查基因命名格式")
# 提取目标基因表达矩阵
target_expr <- expr_matrix[target_genes, ]
# 3. 数据标准化（行Z-score） --------------------------------------------------
scaled_expr <- t(scale(t(target_expr)))

# 4. 准备样本注释 ------------------------------------------------------------
# 1. 创建只包含condition的样本注释
condition_vector <- as.character(colData(dds)$salinity)
col_anno <- HeatmapAnnotation(
  Condition = condition_vector,  # 直接使用condition向量
  col = list(
    Condition = c("Salty" = rgb(173, 31, 31, max = 255), "Fresh" =  rgb(15, 62, 138, max = 255))  # 自定义颜色
  ),
  annotation_legend_param = list(
    Condition = list(title = "Condition")  # 设置图例标题
  ),
  show_annotation_name = TRUE  # 显示注释名称
)

# 5. 创建基因注释（添加差异表达信息） -----------------------------------------
# 添加log2FC和p-value信息
gene_anno <- res[target_genes, c("log2FoldChange", "padj")]
gene_anno$abs_log2FC <- abs(gene_anno$log2FoldChange)
heatmap_colors <- colorRamp2(
  breaks = c(-2, -1, 0, 1, 2),  # 根据数据范围调整
  colors = c("#FFFFFF", "#FFFFCC", "#FFFF00", "#FFA500", "#FF0000")
)
row_anno <- rowAnnotation(
  `|log2FC|` = anno_simple(
    gene_anno$abs_log2FC,
    col = heatmap_colors, # 蓝-白-红渐变
    width = unit(0.5, "cm")
  ),
  `-log10(padj)` = anno_barplot(
    -log10(gene_anno$padj),
    bar_width = 0.6,
    gp = gpar(fill = "grey30", col = NA),
    width = unit(1.5, "cm")
  ),
  show_annotation_name = TRUE,
  annotation_name_rot = 0
)
heatmap_colors <- colorRamp2(
  breaks = c(-2, -1, 0, 1, 2),  # 根据数据范围调整
  colors = c("#FFFFFF", "#FFFFCC", "#FFFF00", "#FFA500", "#FF0000")
)
# 6. 绘制热图 -----------------------------------------------------------------
ht <- Heatmap(
  scaled_expr,
  name = "Z-score", # 图例标题
  col = heatmap_colors, # 热图颜色
  
  # 行列设置
  row_names_gp = gpar(fontsize = 9), # 基因名字体
  column_names_gp = gpar(fontsize = 8), # 样本名字体
  show_column_names = TRUE, # 样本多时隐藏名称
  
  # 聚类设置
  cluster_rows = TRUE, # 基因聚类
  clustering_distance_rows = "pearson", # 使用相关系数聚类
  cluster_columns = TRUE, # 样本聚类
  show_row_dend = FALSE, # 隐藏基因聚类树
  
  # 添加注释
  top_annotation = col_anno, # 列注释（样本分组）
  right_annotation = row_anno, # 行注释（基因DE信息）
  
  # 热图分割
  row_split = 2, # 按基因表达模式分成2组
  row_title = NULL, # 隐藏行标题
  
  # 热图尺寸
  width = unit(15, "cm"), 
  height = unit(15, "cm")
)

# 7. 保存结果 -----------------------------------------------------------------
pdf("DEG_Heatmap.pdf", width = 10, height = 12)
draw(ht, heatmap_legend_side = "right")
dev.off()
