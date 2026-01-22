
arrOutgroupSpp='Gobiopsis_macrostoma,Boleophthalmus_pectinirostris,Glossogobius_giuris,Tridentiger_bifasciatus,Tridentiger_barbatus,Tridentiger_radiatus,Odontamblyopus_rebecca,Z6_R-giurinus_Zhejianglongquanshizhulongzhen';
arrUnusedClades='HN6_R-sp_Hainanhaikou,GX1_R-sp_Guangxinanning'; #c('Aplocheilidae', 'root'); #mark as unused clades "only for hyphy relax"
arrMarkUnusedCladeChildren='F,F'; #c(T , F);

#######################################################################
arrForegroundClades='Clade7'; #mark as foreground
arrMarkForegroundChildren='T';

sTaxonRequirements="min_taxon_requirements_Rform.txt";

nMarkStyle='relax'; #either "codeml" or "relax"
sOutDIR="Relax_Clade7";

mkdir -p $sOutDIR
Rscript labelbranches.R $sOutDIR  $nMarkStyle  $sTaxonRequirements  $arrOutgroupSpp  $arrForegroundClades  $arrMarkForegroundChildren  $arrUnusedClades  $arrMarkUnusedCladeChildren > $sOutDIR/log.txt &


sOutDIR="Codeml_Clade7";
arrMarkForegroundChildren='F';
nMarkStyle='codeml';
mkdir -p $sOutDIR
Rscript labelbranches.R $sOutDIR  $nMarkStyle  $sTaxonRequirements  $arrOutgroupSpp  $arrForegroundClades  $arrMarkForegroundChildren  $arrUnusedClades  $arrMarkUnusedCladeChildren > $sOutDIR/log.txt &

########################################################################

wait
