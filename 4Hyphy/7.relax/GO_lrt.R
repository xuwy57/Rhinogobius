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

arrTargetGenes_clade27 <- intersect(arrTargetGenes_clade2, arrTargetGenes_clade7)
arrTargetGenes_clade25 <- intersect(arrTargetGenes_clade2, arrTargetGenes_clade5)
arrTargetGenes_clade57 <- intersect(arrTargetGenes_clade5, arrTargetGenes_clade7)
arrTargetGenes_clade257 <- intersect(arrTargetGenes_clade2, arrTargetGenes_clade57)
arrAllGenes_clade27 <- intersect(arrAllGenes_clade2, arrAllGenes_clade7)
arrAllGenes_clade25 <- intersect(arrAllGenes_clade2, arrAllGenes_clade5)
arrAllGenes_clade57 <- intersect(arrAllGenes_clade5, arrAllGenes_clade7)
arrAllGenes_clade257 <- intersect(arrAllGenes_clade2, arrAllGenes_clade57)

expressed <- read.table("../../deseq2/sim_genes_lrt.csv", header=T, sep=",", fill = T, quote = "")
expressed <- expressed[,1]
DEGs <- read.table("../../deseq2/siglrt_all_output.csv", header=T, sep=",", fill = T, quote = "")
#cluster2 ↓
#cluster4 ↑
DEGs <- subset(DEGs, DEGs$cluster == '3')
DEGs <- DEGs$genes


arrTargetGenes_clade2_DEG <- intersect(arrTargetGenes_clade2, DEGs)
arrAllGenes_clade2_expressed <- intersect(arrAllGenes_clade2, expressed)
arrTargetGenes_clade5_DEG <- intersect(arrTargetGenes_clade5, DEGs)
arrAllGenes_clade5_expressed <- intersect(arrAllGenes_clade5, expressed)
arrTargetGenes_clade7_DEG <- intersect(arrTargetGenes_clade7, DEGs)
arrAllGenes_clade7_expressed <- intersect(arrAllGenes_clade7, expressed)
arrTargetGenes_clade27_DEG <- intersect(arrTargetGenes_clade27, DEGs)
arrAllGenes_clade27_expressed <- intersect(arrAllGenes_clade27, expressed)
arrTargetGenes_clade57_DEG <- intersect(arrTargetGenes_clade57, DEGs)
arrAllGenes_clade57_expressed <- intersect(arrAllGenes_clade57, expressed)
arrTargetGenes_clade25_DEG <- intersect(arrTargetGenes_clade25, DEGs)
arrAllGenes_clade25_expressed <- intersect(arrAllGenes_clade25, expressed)
arrTargetGenes_clade257_DEG <- intersect(arrTargetGenes_clade257, DEGs)
arrAllGenes_clade257_expressed <- intersect(arrAllGenes_clade257, expressed)

arrTargetGenes_clade2_expressed <- intersect(arrTargetGenes_clade2, expressed)
DEGs_Clade2 <- intersect(arrAllGenes_clade2, DEGs)
arrTargetGenes_clade5_expressed <- intersect(arrTargetGenes_clade5, expressed)
DEGs_Clade5 <- intersect(arrAllGenes_clade5, DEGs)
arrTargetGenes_clade7_expressed <- intersect(arrTargetGenes_clade7, expressed)
DEGs_Clade7 <- intersect(arrAllGenes_clade7, DEGs)



nName <- "Cluster1clade2"
filter_genes <- arrAllGenes_clade1_expressed
arrTargetGenes <- arrTargetGenes_clade1_DEG

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
datOut <- oEnrichRet@result[ oEnrichRet@result$p.adjust <=1, ]
View(fnAddGeneSymbols(datOut))
write.table(fnAddGeneSymbols(datOut), file=sOut, quote = F, sep=",", col.names = T, row.names = F)

go <- enrichGO(gene = arrTargetGenes,
               keyType="GID",
               OrgDb = org.Rformosanus.eg.db,
               ont = "ALL",
               universe = filter_genes,
               pAdjustMethod = "BH",
               pvalueCutoff = 1,
               qvalueCutoff =1)
write.table(go, paste0(nName,"_",sRelaxed,"_fdr", nFDRCutoff, "_p", nRawPCutoff, "_GO_",sCategory,".interprotGO.fdr0.05.csv"), sep = ",", col.names = T, row.names = F, quote = F)

go <- enrichGO(gene = arrTargetGenes,
               keyType="GID",
               OrgDb = org.Rformosanus.v2.eg.db,
               ont = "ALL",
               universe = filter_genes,
               pAdjustMethod = "BH",
               pvalueCutoff = 1,
               qvalueCutoff =1)
write.table(go, paste0(nName,"_",sRelaxed,"_fdr", nFDRCutoff, "_p", nRawPCutoff, "_GO_",sCategory,".interprotGOwithEnsembl.fdr0.05.csv"), sep = ",", col.names = T, row.names = F, quote = F)


length(intersect(arrAllGenes_clade57_expressed, arrTargetGenes_clade7_DEG))

## phyper test
#cluster1
#relaxed
phyper(10, 32, 5139 - 32, 27, lower.tail = FALSE) #25
phyper(34, 54, 7411 - 54, 78, lower.tail = FALSE) #27
phyper(18, 29, 5229 - 29, 51, lower.tail = FALSE) #57
#relaxed
phyper(57, 994, 7424 - 994, 348, lower.tail = FALSE) #2
phyper(29, 578, 5231 - 578, 244, lower.tail = FALSE) #5
phyper(88, 1885, 7978 - 1885, 375, lower.tail = FALSE) #7
#intensified
phyper(4, 117, 7424 - 117, 348, lower.tail = FALSE) #2
phyper(1, 56, 5231 - 56, 244, lower.tail = FALSE) #5
phyper(8, 160, 7978 - 160, 375, lower.tail = FALSE) #7

#cluster2
#relaxed
phyper(13, 20, 5139 - 20, 21, lower.tail = FALSE) #25
phyper(29, 37, 7411 - 37, 65, lower.tail = FALSE) #27
phyper(15, 22, 5229 - 22, 41, lower.tail = FALSE) #57
#relaxed
phyper(42, 994, 7424 - 994, 251, lower.tail = FALSE) #2
phyper(23, 578, 5231 - 578, 191, lower.tail = FALSE) #5
phyper(75, 1885, 7978 - 1885, 270, lower.tail = FALSE) #7
#intensified
phyper(3, 117, 7424 - 117, 251, lower.tail = FALSE) #2
phyper(1, 56, 5231 - 56, 191, lower.tail = FALSE) #5
phyper(2, 160, 7978 - 160, 270, lower.tail = FALSE) #7

#cluster3
#relaxed
phyper(13, 22, 5139 - 22, 23, lower.tail = FALSE) #25
phyper(26, 34, 7411 - 34, 71, lower.tail = FALSE) #27
phyper(19, 23, 5229 - 23, 49, lower.tail = FALSE) #57
#intensified
phyper(1, 4, 5139 - 4, 3, lower.tail = FALSE) #25
phyper(2, 4, 7411 - 4, 8, lower.tail = FALSE) #27
phyper(2, 3, 5229 - 3, 6, lower.tail = FALSE) #57
#relaxed
phyper(34, 994, 7424 - 994, 223, lower.tail = FALSE) #2
phyper(23, 578, 5231 - 578, 143, lower.tail = FALSE) #5
phyper(84, 1885, 7978 - 1885, 242, lower.tail = FALSE) #7
#intensified
phyper(4, 117, 7424 - 117, 223, lower.tail = FALSE) #2
phyper(3, 56, 5231 - 56, 143, lower.tail = FALSE) #5
phyper(8, 160, 7978 - 160, 242, lower.tail = FALSE) #7

#cluster4
#relaxed
phyper(5, 22, 5139 - 22, 17, lower.tail = FALSE) #25
phyper(22, 37, 7411 - 37, 70, lower.tail = FALSE) #27
phyper(7, 17, 5229 - 17, 45, lower.tail = FALSE) #57
#intensified
phyper(1, 7, 7411 - 7, 9, lower.tail = FALSE) #27
#relaxed
phyper(37, 994, 7424 - 994, 367, lower.tail = FALSE) #2
phyper(17, 578, 5231 - 578, 266, lower.tail = FALSE) #5
phyper(77, 1885, 7978 - 1885, 388, lower.tail = FALSE) #7
#intensified
phyper(6, 117, 7424 - 117, 367, lower.tail = FALSE) #2
phyper(3, 56, 5231 - 56, 266, lower.tail = FALSE) #5
phyper(10, 160, 7978 - 160, 388, lower.tail = FALSE) #7

#cluster5
#relaxed
phyper(4, 11, 5139 - 11, 6, lower.tail = FALSE) #25
phyper(11, 16, 7411 - 16, 21, lower.tail = FALSE) #27
phyper(4, 6, 5229 - 6, 15, lower.tail = FALSE) #57
#relaxed
phyper(16, 994, 7424 - 994, 89, lower.tail = FALSE) #2
phyper(6, 578, 5231 - 578, 69, lower.tail = FALSE) #5
phyper(24, 1885, 7978 - 1885, 94, lower.tail = FALSE) #7
#intensified
phyper(2, 117, 7424 - 117, 89, lower.tail = FALSE) #2
phyper(1, 56, 5231 - 56, 69, lower.tail = FALSE) #5
phyper(4, 160, 7978 - 160, 94, lower.tail = FALSE) #7


#cluster6
#relaxed
phyper(5, 7, 5139 - 7, 6, lower.tail = FALSE) #25
phyper(10, 12, 7411 - 12, 17, lower.tail = FALSE) #27
phyper(6, 6, 5229 - 6, 9, lower.tail = FALSE) #57
#relaxed
phyper(14, 994, 7424 - 994, 69, lower.tail = FALSE) #2
phyper(7, 578, 5231 - 578, 42, lower.tail = FALSE) #5
phyper(23, 1885, 7978 - 1885, 75, lower.tail = FALSE) #7
#intensified
phyper(0, 117, 7424 - 117, 69, lower.tail = FALSE) #2
phyper(0, 56, 5231 - 56, 42, lower.tail = FALSE) #5
phyper(1, 160, 7978 - 160, 75, lower.tail = FALSE) #7


#cluster7
#relaxed
phyper(3, 5, 5139 - 5, 5, lower.tail = FALSE) #25
phyper(6, 11, 7411 - 11, 15, lower.tail = FALSE) #27
phyper(2, 5, 5229 - 5, 10, lower.tail = FALSE) #57
#relaxed
phyper(11, 994, 7424 - 994, 62, lower.tail = FALSE) #2
phyper(5, 578, 5231 - 578, 40, lower.tail = FALSE) #5
phyper(15, 1885, 7978 - 1885, 65, lower.tail = FALSE) #7
#intensified
phyper(0, 117, 7424 - 117, 62, lower.tail = FALSE) #2
phyper(0, 56, 5231 - 56, 40, lower.tail = FALSE) #5
phyper(0, 160, 7978 - 160, 65, lower.tail = FALSE) #7
