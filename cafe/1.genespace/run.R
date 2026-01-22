#setwd("~/bahaha_assembly/synteny/genespace/example/");
library(GENESPACE)
runwd <- file.path("/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/cafe5/genespace/rundir")
#list.files(runwd, recursive = T, full.names = F)

gpar <- init_genespace(
  genomeIDs = c('RhinogobiusFormosanus', 'GobiopsisMacrostoma', 'GlossogobiusGiuris',
                'OdontamblyopusRebecca', 
		'TridentigerBifasciatus', 'TridentigerBarbatus', 'TridentigerRadiatus',
		'BoleophthalmusPectinirostris',
		'AmphiprionOcellaris', 'AnabasTestudineus', 
                'AstatotilapiaCalliptera', 'AstyanaxMexicanus',
                'DanioRerio',
                'GadusMorhua', 'GasterosteusAculeatus',
                'LepisosteusOculatus', 'MaylandiaZebra',
                'OreochromisNiloticus', 
                'OryziasLatipes',
                'ScleropagesFormosus', 'SeriolaDumerili', 'StegastesPartitus', 'TakifuguRubripes', 
                'TetraodonNigroviridis', 'XiphophorusMaculatus'),
  speciesIDs = c('Rhinogobius_formosanus', 'Gobiopsis_macrostoma', 'Glossogobius_giuris',
                 'Odontamblyopus_rebecca', 
                'Tridentiger_bifasciatus', 'Tridentiger_barbatus', 'Tridentiger_radiatus',
		'Boleophthalmus_pectinirostris',
		 'Amphiprion_ocellaris', 'Anabas_testudineus',
                 'Astatotilapia_calliptera', 'Astyanax_mexicanus',
                 'Danio_rerio',
                 'Gadus_morhua', 'Gasterosteus_aculeatus',
                 'Lepisosteus_oculatus', 'Maylandia_zebra',
                 'Oreochromis_niloticus', 
                 'Oryzias_latipes',
                 'Scleropages_formosus', 'Seriola_dumerili', 'Stegastes_partitus', 'Takifugu_rubripes', 
                 'Tetraodon_nigroviridis', 'Xiphophorus_maculatus'),
  versionIDs = rep('1', 25),
  ploidy = rep(1, 25),
  diamondMode = "sensitive",
  orthofinderMethod = "default",
  wd = runwd,
  nCores = 48,
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
