setwd("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/cafe5/family_test0125");
library(org.Rformosanus.clade2.eg.db)
library(org.Rformosanus.clade5.eg.db)
library(org.Rformosanus.clade7.eg.db)
library(clusterProfiler)
library(stringr)
library(readr)

clade2_significant_results <- read_csv("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/cafe5/family_test0125/clade2_vs_amphidromous_fdr0.05_results.csv")
clade5_significant_results <- read_csv("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/cafe5/family_test0125/clade5_vs_amphidromous_fdr0.05_results.csv")
clade7_significant_results <- read_csv("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/cafe5/family_test0125/clade7_vs_amphidromous_fdr0.05_results.csv")
sRefSp <- "RhinogobiusFormosanus"
datFam2Gene_clade2 <- read.table("N0_modified_clade2.tsv", header=T, sep="\t", stringsAsFactors = F);
datFam2Gene_clade5 <- read.table("N0_modified_clade5.tsv", header=T, sep="\t", stringsAsFactors = F);
datFam2Gene_clade7 <- read.table("N0_modified_clade7.tsv", header=T, sep="\t", stringsAsFactors = F);

fnFam2Gene <- function(arr, dat) {
  arrGenes <- dat[dat$OG %in% arr, sRefSp];
  return(unique(str_trim(unlist(strsplit(x = arrGenes, split = ',')))));
}

Expansion_Status <- 'Expanded'
Expansion_Status <- 'Contracted'
arrFamClade2 <- clade2_significant_results$Orthogroup[clade2_significant_results$Expansion_Status == Expansion_Status]
arrFamClade5 <- clade5_significant_results$Orthogroup[clade5_significant_results$Expansion_Status == Expansion_Status]
arrFamClade7 <- clade7_significant_results$Orthogroup[clade7_significant_results$Expansion_Status == Expansion_Status]
arrFamClade2Genes <- fnFam2Gene(arrFamClade2, datFam2Gene_clade2);
arrFamClade5Genes <- fnFam2Gene(arrFamClade5, datFam2Gene_clade5);
arrFamClade7Genes <- fnFam2Gene(arrFamClade7, datFam2Gene_clade7);
arrControlGenes_clade2 <- fnFam2Gene(datFam2Gene_clade2$OG, datFam2Gene_clade2)
arrControlGenes_clade5 <- fnFam2Gene(datFam2Gene_clade5$OG, datFam2Gene_clade5)
arrControlGenes_clade7 <- fnFam2Gene(datFam2Gene_clade7$OG, datFam2Gene_clade7)

###enrichment test:
go <- enrichGO(gene = arrFamClade2Genes,
               keyType="GID",
               OrgDb = org.Rformosanus.clade2.eg.db,
               ont = "ALL",
               universe = arrControlGenes_clade2 ,
               pAdjustMethod = "BH",
               pvalueCutoff = 1,
               qvalueCutoff =1)
write.table(go, paste0("Clade2", Expansion_Status, ".interprotGOwithEnsembl.fdr0.05.csv"), sep = ",", col.names = T, row.names = F, quote = F)


Contracted_interprotGOwithEnsembl_fdr0_1_overlapped <- merge(Clade2Contracted_interprotGOwithEnsembl_fdr0_1, Clade7Contracted_interprotGOwithEnsembl_fdr0_1, by = "Description")
Expanded_interprotGOwithEnsembl_fdr0_1_overlapped <- merge(Clade2Expanded_interprotGOwithEnsembl_fdr0_1, Clade7Expanded_interprotGOwithEnsembl_fdr0_1, by = "Description")
Contracted_interprotGOwithEnsembl_fdr0_2_overlapped <- merge(Clade2Contracted_interprotGOwithEnsembl_fdr0_2, Clade7Contracted_interprotGOwithEnsembl_fdr0_2, by = "Description")
Expanded_interprotGOwithEnsembl_fdr0_2_overlapped <- merge(Clade2Expanded_interprotGOwithEnsembl_fdr0_2, Clade7Expanded_interprotGOwithEnsembl_fdr0_2, by = "Description")
write.table(Contracted_interprotGOwithEnsembl_fdr0_1_overlapped, paste0("Contracted_interprotGOwithEnsembl_fdr0_1_overlapped.csv"), sep = ",", col.names = T, row.names = F, quote = F)
write.table(Contracted_interprotGOwithEnsembl_fdr0_2_overlapped, paste0("Contracted_interprotGOwithEnsembl_fdr0_2_overlapped.csv"), sep = ",", col.names = T, row.names = F, quote = F)
write.table(Expanded_interprotGOwithEnsembl_fdr0_1_overlapped, paste0("Expanded_interprotGOwithEnsembl_fdr0_1_overlapped.csv"), sep = ",", col.names = T, row.names = F, quote = F)
write.table(Expanded_interprotGOwithEnsembl_fdr0_2_overlapped, paste0("Expanded_interprotGOwithEnsembl_fdr0_2_overlapped.csv"), sep = ",", col.names = T, row.names = F, quote = F)
