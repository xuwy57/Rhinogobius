setwd("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/GO")
library(stringr)
datOrtho <- read.table("./ref/Orthogroups.tsv", header = T, stringsAsFactors = F, fill=T, sep="\t")
arrSpp <- colnames(datOrtho)
arrSpp <- arrSpp[2:length(arrSpp)]
sSp <- "RhinogobiusFormosanus";
sTERM2GO <- paste0(sSp, "_term2go.tsv");

arrFromSpp <- arrSpp[arrSpp!=sSp]
bAppend<-F;
for(sFromSp in arrFromSpp) {
  sGOMap <- paste0("ref_GOout/",sFromSp, ".tsv");
  
  if (!file.exists(sGOMap)) {
    next;
  }
  cat("Reading GO from ", sFromSp, "\n");
  datGO <- read.table(sGOMap, header=F, sep="\t");
  datOrthoPairMap <- datOrtho[, c(sSp, sFromSp)]
  datOrthoPairMap <- datOrthoPairMap[complete.cases(datOrthoPairMap),];
  for(i in 1:nrow(datOrthoPairMap)) {
    
       arrFromTransID <- str_trim(unlist(strsplit(x = datOrthoPairMap[i, 2], split = ',')));
       arrSpTransID <- str_trim(unlist(strsplit(x = datOrthoPairMap[i, 1], split = ',')));
       arrGOs <- datGO[datGO$V2 %in% arrFromTransID, 'V4'];
       arrGOs <- unique(arrGOs)
       if (length(arrGOs) == 0) {
             next;
       }
       
       for(sSpTransID in arrSpTransID) {
              datToBind <- data.frame(term=arrGOs, gene=sSpTransID, stringsAsFactors = F);
              write.table(datToBind, file = sTERM2GO, append = bAppend, sep = "\t", col.names = (!bAppend), row.names = F, quote = F)
              bAppend<-T;
       }
         
  }
}

datAllGO <- read.table(sTERM2GO, sep="\t", header=T, stringsAsFactors = T);
datAllGO <- datAllGO[!duplicated(datAllGO),];
write.table(datAllGO, file = sTERM2GO, append = F, sep = "\t", col.names = T, row.names = F, quote = F);
