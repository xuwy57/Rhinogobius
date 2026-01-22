setwd("/fast3/group_crf/home/b20xuwy57/APL/compare/cafe/GO/");
library(org.Dr.eg.db)
library(org.Hs.eg.db)
library(clusterProfiler)
library(stringr)


sInDir <- "./runresults_k2_19final/";
sRefSp <- "NothobranchiusFurzeri"
sSp <- "<5>";
nFDRCutoff <- 0.2
sCategory <- 'BP' ; #biological process

datFam2Gene <- read.table("/fast3/group_crf/home/b20xuwy57/APL/compare/cafe/genespace/rundir/orthofinder/Results_Sep19/Phylogenetic_Hierarchical_Orthogroups/N0.tsv", header=T, sep="\t", stringsAsFactors = F);

datChange <- read.table(paste0(sInDir, "/Gamma_change.tab" ), sep="\t", header=T);
datChange <- datChange[, c(1, 7) ];

datBranchProb <- read.table(paste0(sInDir, "/Gamma_branch_probabilities.tab" ), sep="\t", header=T, comment.char = '/', fill = T);

datChange <- merge(datChange , datBranchProb[, c(1,7)], by.x = "FamilyID", by.y=1)

datChange$fdr <- p.adjust(datChange[, 3])

arrExpandedFam <- datChange[datChange$fdr<=nFDRCutoff & datChange[,2] > 0  ,1]
arrShrinkedFam <- datChange[datChange$fdr<=nFDRCutoff & datChange[,2] < 0  ,1]
arrControlFam <- datChange[ ,1]

fnFam2Gene <- function(arr) {
  arrGenes <- datFam2Gene[datFam2Gene$OG %in% arr, sRefSp];
  return(unique(str_trim(unlist(strsplit(x = arrGenes, split = ',')))));
}

arrExpandedGenes <- fnFam2Gene(arrExpandedFam);
arrShrinkedGenes <- fnFam2Gene(arrShrinkedFam);
arrControlGenes <- fnFam2Gene(arrControlFam);



### map & get gene names 
datPGID2Ref <- read.table("~/APL/RNA/for_annotation/synorthos.txt", header = T, stringsAsFactors = F)

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








###enrichment test:

# arrExpandedGenes 
# arrShrinkedGenes 
# arrControlGenes
                                                
oEnrichRet <- enricher(
  gene=arrExpandedGenes,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  universe = arrControlGenes ,
  minGSSize = 10,
  maxGSSize = 5000,
  qvalueCutoff = 1,
  TERM2GENE = datMap ,
  TERM2NAME = go2term(unique(datMap$term) )
)
datOut <- oEnrichRet@result[ oEnrichRet@result$p.adjust <=1, ]
View(fnAddGeneSymbols(datOut))
write.table(fnAddGeneSymbols(datOut), paste0(sSp, "_expandedfam.GO_", sCategory, ".csv"), sep = ",", col.names = T, row.names = F, quote = F);


oEnrichRet <- enricher(
  gene=arrShrinkedGenes,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  universe = arrControlGenes ,
  minGSSize = 10,
  maxGSSize = 5000,
  qvalueCutoff = 1,
  TERM2GENE = datMap ,
  TERM2NAME = go2term(unique(datMap$term) )
)
datOut <- oEnrichRet@result[ oEnrichRet@result$p.adjust <=1, ]
view(fnAddGeneSymbols(datOut))
write.table(fnAddGeneSymbols(datOut), paste0(sSp, "_shrinkedfam.GO_", sCategory, ".csv"), sep = ",", col.names = T, row.names = F, quote = F);

arrUnMappedExpanded <- arrExpandedGenes[!(arrExpandedGenes %in% unique(datMap$gene))];

arrShow <-c();
for(sUnmapped in arrUnMappedExpanded) {
  arrShow <- c(arrShow, grep(sUnmapped, datPGID2Ref$NothobranchiusFurzeri));
}

View(datPGID2Ref[unique(arrShow),])
