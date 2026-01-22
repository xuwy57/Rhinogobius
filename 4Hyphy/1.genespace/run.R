#setwd("~/bahaha_assembly/synteny/genespace/example/");
library(GENESPACE)
runwd <- file.path("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/genespace/rundir")
#list.files(runwd, recursive = T, full.names = F)

gpar <- init_genespace(
  genomeIDs = c('RhinogobiusFormosanus', 'GobiopsisMacrostoma', 'GlossogobiusGiuris',
		'TridentigerBifasciatus', 'TridentigerBarbatus', 'TridentigerRadiatus',
		'BoleophthalmusPectinirostris', 
                'OdontamblyopusRebecca'),
  speciesIDs = c('Rhinogobius_formosanus', 'Gobiopsis_macrostoma', 'Glossogobius_giuris',
                'Tridentiger_bifasciatus', 'Tridentiger_barbatus', 'Tridentiger_radiatus',
		'Boleophthalmus_pectinirostris', 
                 'Odontamblyopus_rebecca'),
  versionIDs = rep('1', 8),
  ploidy = rep(1, 8),
  diamondMode = "sensitive",
  orthofinderMethod = "default",
  wd = runwd,
  nCores = 64,
  minPepLen = 30,
  gffString = "gff",
  pepString = "pep",
  path2orthofinder = "orthofinder",
  path2diamond = "diamond",
  path2mcscanx = "/public/software/MCScanX/",
  rawGenomeDir = file.path(runwd, "rawGenomes"))


parse_annotations(
  gsParam = gpar,
  gffEntryType = "gene",
  gffIdColumn = "locus",
  gffStripText = "locus=",
  headerEntryIndex = 1,
  headerSep = " ",
  headerStripText = "locus=")

gpar <- run_orthofinder(
  gsParam = gpar)

gpar <- synteny(gsParam = gpar)

#arrCol <- rainbow(23);
#regs <- data.frame(
#  genome = rep('MacropodusOpercularis', 23),
#  chr =1:23)

#datInvert <- data.frame(genome="MacropodusHongkongensis", chr=c(1, 4, 8, 10, 6, 18,19,20, 22) )

#ripdat <- plot_riparianHits(
#  gpar,onlyTheseRegions = regs, refChrCols = arrCol,minGenes2plot=50, invertTheseChrs = datInvert,
#  blackBg = F, chrFill = "orange",returnSourceData=T,
#  chrBorder = "grey", useOrder=F, labelTheseGenomes = c('MacropodusOpercularis', 'MacropodusHongkongensis') )

#dump('ripdat', file="ripdata.R");

pg <- pangenome(gpar)

pdf(file="synteny.pdf",width = 10,height = 8)
plot_riparianHits(pg)
dev.off()
