setwd("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/hyphy_1016/GO")
library(org.Dr.eg.db)
library(org.Hs.eg.db)
library(org.Rformosanus.eg.db)
library(org.Rformosanus.v2.eg.db)
library(clusterProfiler)
library(ggplot2)

nFDRCutoff <- 0.05
nRawPCutoff <- 1

datRelax_clade2 <- read.table("../relax/2clade/sum_Clade2.txt", header=F, sep="\t", fill = T, quote = "")
datRelax_clade5 <- read.table("../relax/5clade/sum_Clade5.txt", header=F, sep="\t", fill = T, quote = "")
datRelax_clade7 <- read.table("../relax/7clade/sum_Clade7.txt", header=F, sep="\t", fill = T, quote = "")
datRelax_clade2$fdr <- p.adjust(datRelax_clade2$V5, method = "fdr")
datRelax_clade5$fdr <- p.adjust(datRelax_clade5$V5, method = "fdr")
datRelax_clade7$fdr <- p.adjust(datRelax_clade7$V5, method = "fdr")
datGroupID2TransID <- read.table("../exportAln_genespace/groupid2transid.txt", sep="\t", header=T)
colnames(datGroupID2TransID) <- c("OrthoID", 'Gene');
datRelax_clade2 <- merge(datRelax_clade2, datGroupID2TransID, by.x = "V1", by.y="OrthoID", all.x=T, all.y=F)
datRelax_clade5 <- merge(datRelax_clade5, datGroupID2TransID, by.x = "V1", by.y="OrthoID", all.x=T, all.y=F)
datRelax_clade7 <- merge(datRelax_clade7, datGroupID2TransID, by.x = "V1", by.y="OrthoID", all.x=T, all.y=F)

sRefSp <- "RhinogobiusFormosanus";
sCategory <- 'BP' ; #biological process

datPGID2Ref <- read.table("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/GO/genespace/synorthos.txt", header = T, stringsAsFactors = F)

datRNAID2ZebraFishHuman <- datPGID2Ref[, c(sRefSp, paste0(sRefSp,'.1'), 'Human', 'Human.1', 'Zebrafish', 'Zebrafish.1') ];
datRNAID2ZebraFishHuman[is.na(datRNAID2ZebraFishHuman[, sRefSp]) , sRefSp ] <- datRNAID2ZebraFishHuman[is.na(datRNAID2ZebraFishHuman[, sRefSp]) , paste0(sRefSp,'.1') ]
datRNAID2ZebraFishHuman[is.na(datRNAID2ZebraFishHuman[, 'Human']) , 'Human' ] <- datRNAID2ZebraFishHuman[is.na(datRNAID2ZebraFishHuman[, "Human"]) , 'Human.1' ]
datRNAID2ZebraFishHuman[is.na(datRNAID2ZebraFishHuman[, 'Zebrafish']) , 'Zebrafish' ] <- datRNAID2ZebraFishHuman[is.na(datRNAID2ZebraFishHuman[, 'Zebrafish']) , 'Zebrafish.1' ]
datRNAID2ZebraFishHuman <- datRNAID2ZebraFishHuman[, c(sRefSp, 'Human', 'Zebrafish') ]
datRNAID2ZebraFishHuman <- datRNAID2ZebraFishHuman[!is.na(datRNAID2ZebraFishHuman[,sRefSp]) , ];
datRNAID2ZebraFishHuman <- datRNAID2ZebraFishHuman[!(is.na(datRNAID2ZebraFishHuman[, 'Human']) & is.na(datRNAID2ZebraFishHuman[, 'Zebrafish']) ), ];
datRNAID2ZebraFishHuman$Human <- sub(";.*", "", datRNAID2ZebraFishHuman$Human)
datRNAID2ZebraFishHuman$Zebrafish <- sub(";.*", "", datRNAID2ZebraFishHuman$Zebrafish)
colnames(datRNAID2ZebraFishHuman)[1] <- "Gene";
datOnlyCompHighImpact <- datRNAID2ZebraFishHuman; #merge(datConsurf, datRNAID2ZebraFishHuman, by.x = "Gene", by.y=1 , all.x = T, all.y =F)
datZebrafishTrans2GeneID <- as.data.frame( org.Dr.egENSEMBLTRANS)
colnames(datZebrafishTrans2GeneID)[1] <- "ZebrafishGeneID";
datOnlyCompHighImpact <- merge(datOnlyCompHighImpact, datZebrafishTrans2GeneID, by.x="Zebrafish", by.y = "trans_id", all.x =T, all.y =F);
datHumanTrans2GeneID <- as.data.frame( org.Hs.egENSEMBLTRANS)
colnames(datHumanTrans2GeneID)[1] <- "HumanGeneID";
datOnlyCompHighImpact <- merge(datOnlyCompHighImpact, datHumanTrans2GeneID, by.x="Human", by.y = "trans_id", all.x =T, all.y =F);
datHumanGOMap <- as.data.frame(org.Hs.egGO)
datZebrafishGOMap <- as.data.frame(org.Dr.egGO)
datMap <- datOnlyCompHighImpact[, c('Gene', "ZebrafishGeneID")]
datMap <- datMap[complete.cases(datMap) , ];
datMap <- merge(datMap, datZebrafishGOMap[datZebrafishGOMap$Ontology == sCategory,1:2], by.y="gene_id", by.x="ZebrafishGeneID", all.x=T, all.y=F )
datMap <- datMap[complete.cases(datMap) , c(3,2)];
datMapHs <- datOnlyCompHighImpact[, c('Gene', "HumanGeneID")]
datMapHs <- datMapHs[complete.cases(datMapHs) , ];
datMapHs <- merge(datMapHs, datHumanGOMap[datHumanGOMap$Ontology ==sCategory ,1:2], by.y="gene_id", by.x="HumanGeneID", all.x=T, all.y=F )
datMapHs <- datMapHs[complete.cases(datMapHs) , c(3,2)];
datMap <- rbind(datMap, datMapHs);
datMap <- datMap[!duplicated(datMap),];
colnames(datMap) <- c("term", "gene");

#get gene symbol map
datHumanSymbolMap <- as.data.frame(org.Hs.egSYMBOL)
datZebrafishSymbolMap <- as.data.frame(org.Dr.egSYMBOL)
datSymbol <- datOnlyCompHighImpact[, c('Gene', "ZebrafishGeneID")]
datSymbol <- datSymbol[complete.cases(datSymbol) , ];
datSymbol <- merge(datSymbol, datZebrafishSymbolMap[,1:2], by.y="gene_id", by.x="ZebrafishGeneID", all.x=T, all.y=F )
datSymbol <- datSymbol[complete.cases(datSymbol) , c(3,2)];
datSymbolHs <- datOnlyCompHighImpact[, c('Gene', "HumanGeneID")]
datSymbolHs <- datSymbolHs[complete.cases(datSymbolHs) , ];
datSymbolHs <- merge(datSymbolHs, datHumanSymbolMap[,1:2], by.y="gene_id", by.x="HumanGeneID", all.x=T, all.y=F )
datSymbolHs <- datSymbolHs[complete.cases(datSymbolHs) , c(3,2)];
datSymbol <- rbind(datSymbol , datSymbolHs[ !(datSymbolHs$Gene %in% datSymbol$Gene),] );
datSymbol <- datSymbol[!duplicated(datSymbol),];
fnAddGeneSymbols <- function(datOut) {
  arrSymbolStr <- c();
  if (nrow(datOut) ==0) {
    return(datOut)
  }
  for(nRow in 1:nrow(datOut)) {
    arrPGIDs <- unlist(strsplit(datOut$geneID[nRow] ,'/'));
    arrGeneSymbols <- datSymbol[datSymbol$Gene %in% arrPGIDs, 1];
    sStr <- "";
    if (length(arrGeneSymbols)>0) {
      sStr <- paste(arrGeneSymbols, collapse = '/');
    }
    arrSymbolStr <- c(arrSymbolStr, sStr);
  }
  datOut$genesymbols <- arrSymbolStr;
  return(datOut)
}



bRelaxed <- T; #True  - relaxed genes, False - intensified genes
bRelaxed <- F; #True  - relaxed genes, False - intensified genes
#for(sHap in arrSamples) {
arrIsRelaxed <- datRelax_clade2$V9 < 1;
sRelaxed <- "relaxed";
if (  !bRelaxed ) {
  arrIsRelaxed <- datRelax_clade2$V9 > 1;
  sRelaxed <- "intensified";
}
arrTargetGenes_clade2 <- unique(datRelax_clade2[ arrIsRelaxed , 'Gene']);
arrTargetGenes_clade2 <- unique(datRelax_clade2[ datRelax_clade2$V5<= nRawPCutoff & datRelax_clade2$fdr <= nFDRCutoff & arrIsRelaxed , 'Gene']);
arrTargetGenes_clade2 <- as.character(arrTargetGenes_clade2)
arrAllGenes_clade2 <- unique(datRelax_clade2[, 'Gene'])

bRelaxed <- T; #True  - relaxed genes, False - intensified genes
#bRelaxed <- F; #True  - relaxed genes, False - intensified genes
#for(sHap in arrSamples) {
arrIsRelaxed <- datRelax_clade5$V9 < 1;
sRelaxed <- "relaxed";
if (  !bRelaxed ) {
  arrIsRelaxed <- datRelax_clade5$V9 > 1;
  sRelaxed <- "intensified";
}
arrTargetGenes_clade5 <- unique(datRelax_clade5[ arrIsRelaxed , 'Gene']);
arrTargetGenes_clade5 <- unique(datRelax_clade5[ datRelax_clade5$V5<= nRawPCutoff & datRelax_clade5$fdr <= nFDRCutoff & arrIsRelaxed , 'Gene']);
arrTargetGenes_clade5 <- as.character(arrTargetGenes_clade5)
arrAllGenes_clade5 <- unique(datRelax_clade5[, 'Gene'])

bRelaxed <- T; #True  - relaxed genes, False - intensified genes
#bRelaxed <- F; #True  - relaxed genes, False - intensified genes
#for(sHap in arrSamples) {
arrIsRelaxed <- datRelax_clade7$V9 < 1;
sRelaxed <- "relaxed";
if (  !bRelaxed ) {
  arrIsRelaxed <- datRelax_clade7$V9 > 1;
  sRelaxed <- "intensified";
}
arrTargetGenes_clade7 <- unique(datRelax_clade7[ arrIsRelaxed , 'Gene']);
arrTargetGenes_clade7 <- unique(datRelax_clade7[ datRelax_clade7$V5<= nRawPCutoff & datRelax_clade7$fdr <= nFDRCutoff & arrIsRelaxed , 'Gene']);
arrTargetGenes_clade7 <- as.character(arrTargetGenes_clade7)
arrAllGenes_clade7 <- unique(datRelax_clade7[, 'Gene'])


#snpeff分析部分
arrTargetGenes_clade27 <- intersect(arrTargetGenes_clade2, arrTargetGenes_clade7)
arrTargetGenes_clade25 <- intersect(arrTargetGenes_clade2, arrTargetGenes_clade5)
arrTargetGenes_clade57 <- intersect(arrTargetGenes_clade5, arrTargetGenes_clade7)
arrTargetGenes_clade257 <- intersect(arrTargetGenes_clade2, arrTargetGenes_clade57)
arrAllGenes_clade27 <- intersect(arrAllGenes_clade2, arrAllGenes_clade7)
arrAllGenes_clade25 <- intersect(arrAllGenes_clade2, arrAllGenes_clade5)
arrAllGenes_clade57 <- intersect(arrAllGenes_clade5, arrAllGenes_clade7)
arrAllGenes_clade257 <- intersect(arrAllGenes_clade2, arrAllGenes_clade57)

all <- union(arrAllGenes_clade2, arrAllGenes_clade5, arrAllGenes_clade7)

expressed <- read.table("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/deseq2/sim_genes.csv", header=T, sep=",", fill = T, quote = "")
expressed <- expressed[,1]
expressed_up <- read.table("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/deseq2/sim_up.csv", header=T, sep=",", fill = T, quote = "")
expressed_up <- expressed_up[,1]
expressed_down <- read.table("./fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/deseq2/sim_down.csv", header=T, sep=",", fill = T, quote = "")
expressed_down <- expressed_down[,1]
DEGs <- read.table("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/deseq2/DEGs.csv", header=T, sep=",", fill = T, quote = "")
DEGs <- DEGs[,1]
DEGs_up <- read.table("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/deseq2/DEGs_up.csv", header=T, sep=",", fill = T, quote = "")
DEGs_up <- DEGs_up[,1]
DEGs_down <- read.table("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/deseq2/DEGs_down.csv", header=T, sep=",", fill = T, quote = "")
DEGs_down <- DEGs_down[,1]

all <- intersect(all, expressed)

arrTargetGenes_clade2_DEG <- intersect(arrTargetGenes_clade2, DEGs)
arrTargetGenes_clade2_DEGup <- intersect(arrTargetGenes_clade2, DEGs_up)
arrTargetGenes_clade2_DEGdown <- intersect(arrTargetGenes_clade2, DEGs_down)
arrTargetGenes_clade2_expressedup <- intersect(arrTargetGenes_clade2, expressed_up)
arrTargetGenes_clade2_expresseddown <- intersect(arrTargetGenes_clade2, expressed_down)
arrAllGenes_clade2_expressed <- intersect(arrAllGenes_clade2, expressed)
arrAllGenes_clade2_expressedup <- intersect(arrAllGenes_clade2, expressed_up)
arrAllGenes_clade2_expresseddown <- intersect(arrAllGenes_clade2, expressed_down)
arrTargetGenes_clade5_DEG <- intersect(arrTargetGenes_clade5, DEGs)
arrTargetGenes_clade5_DEGup <- intersect(arrTargetGenes_clade5, DEGs_up)
arrTargetGenes_clade5_DEGdown <- intersect(arrTargetGenes_clade5, DEGs_down)
arrTargetGenes_clade5_expressedup <- intersect(arrTargetGenes_clade5, expressed_up)
arrTargetGenes_clade5_expresseddown <- intersect(arrTargetGenes_clade5, expressed_down)
arrAllGenes_clade5_expressed <- intersect(arrAllGenes_clade5, expressed)
arrAllGenes_clade5_expressedup <- intersect(arrAllGenes_clade5, expressed_up)
arrAllGenes_clade5_expresseddown <- intersect(arrAllGenes_clade5, expressed_down)
arrTargetGenes_clade7_DEG <- intersect(arrTargetGenes_clade7, DEGs)
arrTargetGenes_clade7_DEGup <- intersect(arrTargetGenes_clade7, DEGs_up)
arrTargetGenes_clade7_DEGdown <- intersect(arrTargetGenes_clade7, DEGs_down)
arrTargetGenes_clade7_expressedup <- intersect(arrTargetGenes_clade7, expressed_up)
arrTargetGenes_clade7_expresseddown <- intersect(arrTargetGenes_clade7, expressed_down)
arrAllGenes_clade7_expressed <- intersect(arrAllGenes_clade7, expressed)
arrAllGenes_clade7_expressedup <- intersect(arrAllGenes_clade7, expressed_up)
arrAllGenes_clade7_expresseddown <- intersect(arrAllGenes_clade7, expressed_down)
arrTargetGenes_clade27_DEG <- intersect(arrTargetGenes_clade27, DEGs)
arrTargetGenes_clade27_DEGup <- intersect(arrTargetGenes_clade27, DEGs_up)
arrTargetGenes_clade27_DEGdown <- intersect(arrTargetGenes_clade27, DEGs_down)
arrAllGenes_clade27_expressed <- intersect(arrAllGenes_clade27, expressed)
arrAllGenes_clade27_expressedup <- intersect(arrAllGenes_clade27, expressed_up)
arrAllGenes_clade27_expresseddown <- intersect(arrAllGenes_clade27, expressed_down)
arrTargetGenes_clade57_DEG <- intersect(arrTargetGenes_clade57, DEGs)
arrTargetGenes_clade57_DEGup <- intersect(arrTargetGenes_clade57, DEGs_up)
arrTargetGenes_clade57_DEGdown <- intersect(arrTargetGenes_clade57, DEGs_down)
arrAllGenes_clade57_expressed <- intersect(arrAllGenes_clade57, expressed)
arrAllGenes_clade57_expressedup <- intersect(arrAllGenes_clade57, expressed_up)
arrAllGenes_clade57_expresseddown <- intersect(arrAllGenes_clade57, expressed_down)
arrTargetGenes_clade25_DEG <- intersect(arrTargetGenes_clade25, DEGs)
arrTargetGenes_clade25_DEGup <- intersect(arrTargetGenes_clade25, DEGs_up)
arrTargetGenes_clade25_DEGdown <- intersect(arrTargetGenes_clade25, DEGs_down)
arrAllGenes_clade25_expressed <- intersect(arrAllGenes_clade25, expressed)
arrAllGenes_clade25_expressedup <- intersect(arrAllGenes_clade25, expressed_up)
arrAllGenes_clade25_expresseddown <- intersect(arrAllGenes_clade25, expressed_down)
arrTargetGenes_clade257_DEG <- intersect(arrTargetGenes_clade257, DEGs)
arrTargetGenes_clade257_DEGup <- intersect(arrTargetGenes_clade257, DEGs_up)
arrTargetGenes_clade257_DEGdown <- intersect(arrTargetGenes_clade257, DEGs_down)
arrAllGenes_clade257_expressed <- intersect(arrAllGenes_clade257, expressed)
arrAllGenes_clade257_expressedup <- intersect(arrAllGenes_clade257, expressed_up)
arrAllGenes_clade257_expresseddown <- intersect(arrAllGenes_clade257, expressed_down)


arrTargetGenes_clade2_expressed <- intersect(arrTargetGenes_clade2, expressed)
DEGs_Clade2 <- intersect(arrAllGenes_clade2, DEGs)
arrTargetGenes_clade5_expressed <- intersect(arrTargetGenes_clade5, expressed)
DEGs_Clade5 <- intersect(arrAllGenes_clade5, DEGs)
arrTargetGenes_clade7_expressed <- intersect(arrTargetGenes_clade7, expressed)
DEGs_Clade7 <- intersect(arrAllGenes_clade7, DEGs)

DEGsup_Clade2 <- intersect(arrAllGenes_clade2, DEGs_up) 
DEGsup_Clade5 <- intersect(arrAllGenes_clade5, DEGs_up) 
DEGsup_Clade7 <- intersect(arrAllGenes_clade7, DEGs_up) 
DEGsdown_Clade2 <- intersect(arrAllGenes_clade2, DEGs_down) 
DEGsdown_Clade5 <- intersect(arrAllGenes_clade5, DEGs_down) 
DEGsdown_Clade7 <- intersect(arrAllGenes_clade7, DEGs_down) 
all <- intersect(all_genes, expressed)
setwd("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/snpEff/GO")
genes <- all_genes
nName <- "clade7highimpact"
nName <- "clade2highimpact_DEG"
filter_genes <- genes
arrTargetGenes <- arrTargetGenes_clade7

oEnrichRet <- enricher(
  arrTargetGenes,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  filter_genes,
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 1,
  TERM2GENE = datMap ,
  TERM2NAME = go2term(unique(datMap$term) )
)
sOut <- paste0(nName,"_",sRelaxed,"_fdr", nFDRCutoff, "_p", nRawPCutoff, "_GO_",sCategory,".csv");
sOut <- paste0(nName,"_", "p0.05",sCategory,".csv");
datOut <- oEnrichRet@result[ oEnrichRet@result$p.adjust <=1, ]
View(fnAddGeneSymbols(datOut))
write.table(fnAddGeneSymbols(datOut), file=sOut, quote = F, sep=",", col.names = T, row.names = F)


length(intersect(arrAllGenes_clade257, arrTargetGenes_clade2))

## phyper test
phyper(599, 1490, 12651 - 1490, 1294, lower.tail = FALSE) #25
phyper(1511, 2296, 16899 - 2296, 3823, lower.tail = FALSE) #27
phyper(797, 1323, 12876 - 1323, 2640, lower.tail = FALSE) #57
phyper(16, 172, 12651 - 172, 120, lower.tail = FALSE) #25
phyper(30, 236, 16899 - 236, 299, lower.tail = FALSE) #27
phyper(19, 121, 12876 - 121, 212, lower.tail = FALSE) #57

phyper(15, 42, 5754 - 42, 35, lower.tail = FALSE) #25
phyper(46, 65, 8230 - 65, 122, lower.tail = FALSE) #27
phyper(21, 35, 5859 - 35, 80, lower.tail = FALSE) #57
phyper(1, 5, 5754 - 5, 5, lower.tail = FALSE) #25
phyper(1, 6, 8230 - 6, 13, lower.tail = FALSE) #27
phyper(1, 5, 5859 - 5, 10, lower.tail = FALSE) #57

phyper(69, 1111, 8243 - 1111, 581, lower.tail = FALSE) #2
phyper(36, 623, 5861 - 623, 431, lower.tail = FALSE) #5
phyper(144, 2079, 8890 - 2079, 626, lower.tail = FALSE) #7


##snpeff phyper test 
#with DEG
phyper(61, 736, 16213 - 736, 1182, lower.tail = FALSE) #2
phyper(52, 559, 16213 - 559, 1182, lower.tail = FALSE) #5
phyper(104, 1126, 16213 - 1126, 1182, lower.tail = FALSE) #7

# 多重检验校正
adjusted_p <- p.adjust(c(p_AB, p_AC, p_BC), method = "BH")
names(adjusted_p) <- c("A-B", "A-C", "B-C")
print(adjusted_p)


# 清除旧版本残留
remove.packages("SuperExactTest")
# 重新安装
devtools::install_github("mw201608/SuperExactTest")
library(SuperExactTest)
library(gridExtra)
library(ggplot2)  # 用于后期修饰
set.seed(1234)
# 创建存储图形的列表
plot_list <- list()
#generate random strings
n=100000
input = list(clade2 = arrTargetGenes_clade2,
             clade5 = arrTargetGenes_clade5,
             clade7 = arrTargetGenes_clade7)
input = list(clade2 = arrTargetGenes_clade2_DEG,
             clade5 = arrTargetGenes_clade5_DEG,
             clade7 = arrTargetGenes_clade7_DEG)
Result=supertest(input,n=n)
res_adj <- summary(Result, method="BH")
png(tempfile(), width=1200, height=1000, res=200)  # 使用临时文件
par(mar=c(3, 12, 4, 12))  # 调整边距适应组合布局
order2=c("001","010","100","011","101","110","111")
print(order2)
plot(Result, Layout="landscape", sort.by=order2, keep=FALSE,
#     bar.split=c(50,400), 
     show.elements=FALSE, elements.cex=0.7,
     elements.list=subset(summary(Result)$Table,Observed.Overlap <= 20),
     show.expected.overlap=TRUE,expected.overlap.style="hatchedBox",color.on = NULL, 
     color.expected.overlap='red',
     elements.list = subset(summary(Result)$Table, Observed.Overlap <= 20),
     title = 'b. Relaxed & DEG') # 调整图例位置
plot_list[[2]] <- grid::grid.grab()  # 捕获当前图形
dev.off()
# 组合图形并添加全局元素
final_plot <- grid.arrange(
  grobs = plot_list,
  ncol = 2,  # 2x2布局
  top = textGrob("Comparative Gene Set Intersection Analysis", 
                 gp = gpar(fontsize=18, fontface="bold")),
#  bottom = textGrob("Figure 1. Four independent analyses showing...", 
#                    gp = gpar(fontsize=10)),
  padding = unit(1.5, "cm")  # 整体边距
)
# 输出高分辨率图片
ggsave("fig2.tiff", final_plot, 
       width=16, height=12, units="in", dpi=300, 
       compression = "lzw")
ggsave("fig2.pdf", device=cairo_pdf)





#氨基酸突变检验
library(dplyr)
#洄游型7，clade2 16，clade5 7，clade7 17
clade2 <- read.table("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/plot/plot_protein/ALL/out_filter_new/c2.stats.txt", header=T, sep="\t", fill = T, quote = "")
clade5 <- read.table("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/plot/plot_protein/ALL/out_filter_new/c5.stats.txt", header=T, sep="\t", fill = T, quote = "")
clade7 <- read.table("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/plot/plot_protein/ALL/out_filter_new/c7.stats.txt", header=T, sep="\t", fill = T, quote = "")
colnames(clade2) <- c("genes", "loc", "counts_a", "counts_c2", "allcounts_a", "allcounts_c2")  
colnames(clade5) <- c("genes", "loc", "counts_a", "counts_c5", "allcounts_a", "allcounts_c5")  
colnames(clade7) <- c("genes", "loc", "counts_a", "counts_c7", "allcounts_a", "allcounts_c7")  
clade2$ID <- paste(clade2[[1]], clade2[[2]], sep = "_")  # 使用$和paste创建新列
clade5$ID <- paste(clade5[[1]], clade5[[2]], sep = "_")  # 使用$和paste创建新列
clade7$ID <- paste(clade7[[1]], clade7[[2]], sep = "_")  # 使用$和paste创建新列
# 9       13      6       14      2577735
filtered_clade2 <- clade2 %>%
  filter(counts_a >= 9, counts_c2 >= 13) %>%
  pull(ID)
filtered_clade5 <- clade5 %>%
  filter(counts_a >= 9, counts_c5 >= 6) %>%
  pull(ID)
filtered_clade7 <- clade7 %>%
  filter(counts_a >= 9, counts_c7 >= 14) %>%
  pull(ID)
filtered_clade25 <- intersect(filtered_clade2, filtered_clade5)
filtered_clade25 <- str_remove(filtered_clade25, "_[^_]*$")
filtered_clade25_gene <- unique(filtered_clade25)
filtered_clade257 <- intersect(filtered_clade2, filtered_clade57)
filtered_clade257 <- str_remove(filtered_clade257, "_[^_]*$")
filtered_clade257_gene <- unique(filtered_clade257)

filtered_clade2_gene <- clade2 %>%
  filter(counts_a >= 9, counts_c2 >= 13) %>%
  pull(genes)
filtered_clade2_gene <- unique(filtered_clade2_gene)
filtered_clade5_gene <- clade5 %>%
  filter(counts_a >= 9, counts_c5 >= 6) %>%
  pull(genes)
filtered_clade5_gene <- unique(filtered_clade5_gene)
filtered_clade7_gene <- clade7 %>%
  filter(counts_a >= 9, counts_c7 >= 14) %>%
  pull(genes)
filtered_clade7_gene <- unique(filtered_clade7_gene)

arrTargetGenes_filteredAAclade2<- intersect(arrTargetGenes_clade2, filtered_clade2_gene)
arrTargetGenes_filteredAAclade5 <- intersect(arrTargetGenes_clade5, filtered_clade5_gene)
arrTargetGenes_filteredAAclade7 <- intersect(arrTargetGenes_clade7, filtered_clade7_gene)
test_overlap_filteredAA <-  intersect(filtered_clade2_gene, arrAllGenes_filteredAAclade2)
test_overlap_arrTargetGenes <-  intersect(arrTargetGenes_clade2, arrAllGenes_filteredAAclade2)
phyper(449, 1254, 7779 - 1254, 2390, lower.tail = FALSE) #2
phyper(231, 683, 5430 - 683, 1416, lower.tail = FALSE) #5
phyper(2000, 3346, 13258 - 3346, 7702, lower.tail = FALSE) #7
phyper(30, 110, 7779 - 110, 2390, lower.tail = FALSE) #2
phyper(8, 40, 5430 - 40, 1416, lower.tail = FALSE) #5
phyper(139, 219, 13258 - 219, 7702, lower.tail = FALSE) #7

arrTargetGenes_filteredAAclade25<- intersect(arrTargetGenes_clade25, filtered_clade25_gene)
phyper(5, 141, 12651 - 141, 599, lower.tail = FALSE) #7

arrTargetGenes_filteredAAclade257<- intersect(arrTargetGenes_clade257, filtered_clade257_gene)
test_overlap_filteredAA <-  intersect(filtered_clade27_gene, arrAllGenes_filteredAAclade27)
test_overlap_arrTargetGenes <-  intersect(arrTargetGenes_clade25, arrAllGenes_clade25)
phyper(30, 243, 16723 - 243, 818, lower.tail = FALSE) #7

write.table(arrTargetGenes_filteredAAclade257, file="/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/plot/plot5/alloverlap.txt", quote = F, sep=",", col.names = T, row.names = F)

arrTargetGenes_filteredAAclade2_DEG<- intersect(arrTargetGenes_filteredAAclade2, DEGs)
arrTargetGenes_filteredAAclade5_DEG<- intersect(arrTargetGenes_filteredAAclade5, DEGs)
arrTargetGenes_filteredAAclade7_DEG<- intersect(arrTargetGenes_filteredAAclade7, DEGs)


#dp_scan检验
# 读取mRNA ID列表
mrna_ids <- readLines("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/dpscan/results_mrna_ids.txt")
# 去重（虽然Python已处理，但R中再次确认）
unique_mrnas <- unique(mrna_ids)
dp_intersectall <- intersect(arrTargetGenes_clade257, unique_mrnas)

# 参数设置
n_sim <- 100000  # 模拟次数
observed_overlap <- 497  # 实际三向交集大小
sim_overlaps <- numeric(n_sim)
# 蒙特卡洛模拟
set.seed(123)
for (i in 1:n_sim) {
  sim_A <- sample(arrAllGenes_clade257, 1490)  # 随机生成集合A
  sim_B <- sample(arrAllGenes_clade257, 1294)  # 随机生成集合B
  sim_C <- sample(arrAllGenes_clade257, 2584)  # 随机生成集合C
  sim_overlaps[i] <- length(Reduce(intersect, list(sim_A, sim_B, sim_C)))
}
# 计算 p 值（单尾检验，观测值是否显著大于随机）
p_value <- sum(sim_overlaps >= observed_overlap) / n_sim
cat("蒙特卡洛模拟 p 值:", p_value, "\n")
# 可视化模拟分布
hist(sim_overlaps, breaks = 30, main = "随机模拟的三向交集大小分布",
     xlab = "重叠数量", col = "lightblue", xlim = c(0, 500))
abline(v = observed_overlap, col = "red", lwd = 2)



#卡方检验 都不显著
#clade7和DEG上调
# 数据准备
#上下调为显著的
relaxed_up <- 78       # 放松选择基因中上调的数目
relaxed_non_up <- 66    # 放松选择基因中未上调的数目
non_relaxed_up <- 634-78   # 非放松选择基因中上调的数目
non_relaxed_non_up <- 548-66 # 非放松选择基因中未上调的数目
# 构建2x2表格
contingency_table <- matrix(c(relaxed_up, relaxed_non_up, non_relaxed_up, non_relaxed_non_up),
                            nrow = 2,
                            byrow = TRUE,
                            dimnames = list(c("Relaxed", "Non_Relaxed"), c("Upregulated", "Non_Upregulated")))
# 卡方检验
chi_result <- chisq.test(contingency_table)
cat("\n contingency_table clade7:\n")
print(contingency_table)
# 输出卡方检验结果
print("卡方检验结果：")
print(chi_result)
# 如果需要单独提取p值
p_value <- chi_result$p.value
cat("P值：", p_value, "\n")


#卡方检验 都不显著
# Clade2 数据
clade2 <- matrix(c(69, 6, 2299-69, 236-6), nrow = 2, byrow = TRUE)
rownames(clade2) <- c("DEG", "NonDEG")
colnames(clade2) <- c("Relaxed", "Intensified")
# Clade5 数据
clade5 <- matrix(c(36, 5, 1323-36, 121-5), nrow = 2, byrow = TRUE)
rownames(clade5) <- c("DEG", "NonDEG")
colnames(clade5) <- c("Relaxed", "Intensified")
# Clade7 数据
clade7 <- matrix(c(144, 14, 4263-144, 320-14), nrow = 2, byrow = TRUE)
rownames(clade7) <- c("DEG", "NonDEG")
colnames(clade7) <- c("Relaxed", "Intensified")
# 分别进行卡方检验
chisq_clade2 <- chisq.test(clade2)
chisq_clade5 <- chisq.test(clade5)
chisq_clade7 <- chisq.test(clade7)
# 输出结果
cat("Clade2 卡方检验结果:\n")
print(chisq_clade2)
cat("\nClade5 卡方检验结果:\n")
print(chisq_clade5)
cat("\nClade7 卡方检验结果:\n")
print(chisq_clade7)