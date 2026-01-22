.libPaths("/public3/group_crf/home/b20xuwy57/R/x86_64-pc-linux-gnu-library/")
setwd("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/deseq2")
library(org.Dr.eg.db)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)

nFDRCutoff <- 0.05
nRawPCutoff <- 1
sRefSp <- "RhinogobiusFormosanus";
sCategory <- 'MF' ; #biological process


datPGID2Ref <- read.table("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/GO/genespace/synorthos.txt", header = T, stringsAsFactors = F)

#datRNAID2ZebraFishHuman <- datPGID2Ref[, c('Human', sRefSp, 'Zebrafish', 'Human.1', paste0(sRefSp,'.1'), 'Zebrafish.1') ];
datRNAID2ZebraFishHuman <- datPGID2Ref[, c(sRefSp, paste0(sRefSp,'.1'), 'Human', 'Zebrafish', 'Zebrafish.1') ];
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

#datSymbol$Gene <- gsub("-T1", "", datSymbol$Gene)

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


library(stringr)
library(dplyr)
# 将 filter_genes 中的字符串连接成一个正则表达式
filter_genes <- expressed
pattern <- paste(filter_genes, collapse = "|")
# 使用 str_detect 筛选数据框中gene列含有 filter_genes 中的字符串
filtered_data <- datMap %>%
  filter(str_detect(gene, pattern))
datMap <- filtered_data

darDEG_fdr <- DEG@rownames
darDEG_fdr <- DEGs_up@rownames
darDEG_fdr <- DEGs_down@rownames

darDEG_fdr <- group$genes

#datMap$gene <- gsub("-T1", "", datMap$gene)

go_rich <- enricher(
    darDEG_fdr,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    #qvalueCutoff = 1,
    TERM2GENE = datMap ,
    TERM2NAME = go2term(unique(datMap$term))
  )

  datOut <- go_rich@result[ go_rich@result$p.adjust <=0.5, ]
  datOut <- go_rich@result
  View(fnAddGeneSymbols(datOut))
  
  write.table(fnAddGeneSymbols(datOut),file = "GO_DEGdown_MF.csv", quote = F, sep=",", col.names = T, row.names = F)



    
p1<-barplot(go_rich,showCategory = 10,title="EnrichmentGO_MF")+
  scale_y_discrete(labels=function(x) strwrap(x,width = 100))
p1
ggsave(p1,file= "./GO_sigLRT_genes30_MF_barplot.png")


p2<- dotplot(go_rich,showCategory = 13,title="EnrichmentGO_MF_dot")+
  scale_y_discrete(labels=function(x) strwrap(x,width = 100))
p2
ggsave(p2,file= "./GO_sigLRT_genes_MF_dotplot.png")

cnetplot(go_rich,showCategory = 3)
res1 <- enrichplot::pairwise_termsim(go_rich)
emapplot(res1)
