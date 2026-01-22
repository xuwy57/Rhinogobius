.libPaths("/public3/group_crf/home/b20xuwy57/R/x86_64-pc-linux-gnu-library/")
setwd("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/deseq2")
library(DESeq2)  
library(tidyverse)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(DEGreport)

count.data <- read.table("counts.txt.renamed.txt", header = T, row.names = 1L)
keep_feature <- rowSums(count.data) > 17 
table(keep_feature)  
count.data <- count.data[keep_feature, ] 

coldata_a<-data.frame(name=colnames(count.data),
                      vs=c(rep("Fresh",9),rep("Salinity",8)),
                      salinity=c(rep("0.Fresh",9),rep("1.5ppt",2),rep("2.15ppt",3),rep("3.20ppt",3)))

coldata_a$vs <- as.factor(coldata_a$vs)
dds_salinity <- DESeqDataSetFromMatrix(countData=count.data, colData=coldata_a, design=~ salinity)
dds_salinity <- DESeq(dds_salinity, test="LRT", reduced = ~ 1)
res_LRT <- results(dds_salinity)
rld_mat <- rlog(dds_salinity, blind = FALSE)

# Create a tibble for LRT results
res_LRT_tb <- res_LRT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
# Subset to return genes with padj < 0.05
sigLRT_genes <- res_LRT_tb %>% 
  dplyr::filter(padj < 0.05)
# Get number of significant genes
nrow(sigLRT_genes)
write.table(sigLRT_genes,file = "sigLRT_genes.csv", quote = F, sep=",", col.names = T, row.names = F)

# Obtain rlog values for those significant genes
rld_mat <- rld_mat@assays@data@listData[[1]]
cluster_rlog <- rld_mat[sigLRT_genes$gene, ]
rownames(coldata_a) <- colnames(count.data)
coldata_a$salinity <- factor(coldata_a$salinity)
pdf("cluster_sigLRT_all.pdf", width=12, height=12)
clusters <- degPatterns(cluster_rlog, metadata = coldata_a, time = "salinity", col = NULL)
dev.off()

clusters$df$genes <- gsub('[.]', '-', clusters$df$genes)
group <- clusters$df %>%
  dplyr::filter(cluster == 4)
group$genes <- gsub('[.]', '-', group$genes)

all_group <- clusters$df
write.csv(all_group, "siglrt_all_output.csv", row.names = FALSE, quote = F)


counts <- as.data.frame(counts(dds_salinity, normalized=TRUE))
counts$sum <- rowSums(counts)
counts_sig <- subset(counts, rownames(counts) %in% siglrt_all_output$genes)
sum_sig <- counts_sig$sum
sum_sig <- sort(sum_sig, decreasing=T)
index <- round(length(sum_sig) * 95 / 100)
min_5 <- sum_sig[index]
counts_sim <- subset(counts, min_5 <= sum)
counts_sim_id <- rownames(counts_sim)
write.csv(counts_sim,file = "sim_genes_lrt.csv", quote = F)













