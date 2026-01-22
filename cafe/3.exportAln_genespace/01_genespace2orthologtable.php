<?php
$sGeneSpace = "/data2/projects/wxu/Rhinogobius_spp/cafe5/exportAln_genespace/gffWithOgs.txt.gz";
$sGeneSymbolMap = "/data2/projects/wxu/Rhinogobius_spp/cafe5/assigngenesymbols_orthofinder/orthofinder_genesymbols.tsv";
$sRNAID2SpGeneIDMap = "/data2/projects/wxu/Rhinogobius_spp/cafe5/exportAln_genespace/preProt/allproteins.fa";
$arrTaxTranslate = array('RhinogobiusFormosanus' => 'Rhinogobius_formosanus',
	'OdontamblyopusRebecca' => 'Odontamblyopus_rebecca',
#        'PeriophthalmusMagnuspinnatus' => 'Periophthalmus_magnuspinnatus',
#        'RhinogobiusSimilis' => 'Rhinogobius_similis',
	'BoleophthalmusPectinirostris' => 'Boleophthalmus_pectinirostris',
	'GobiopsisMacrostoma' => 'Gobiopsis_macrostoma',
	'GlossogobiusGiuris' => 'Glossogobius_giuris',
	'TridentigerBifasciatus' => 'Tridentiger_bifasciatus',
	'TridentigerBarbatus' => 'Tridentiger_barbatus',
	'TridentigerRadiatus' => 'Tridentiger_radiatus',
	'AmphiprionOcellaris' => 'Amphiprion_ocellaris',
	'AnabasTestudineus' => 'Anabas_testudineus',
	'AstatotilapiaCalliptera' => 'Astatotilapia_calliptera',
	'AstyanaxMexicanus' => 'Astyanax_mexicanus',
	'DanioRerio' => 'Danio_rerio',
	'GadusMorhua' => 'Gadus_morhua',
	'GasterosteusAculeatus' => 'Gasterosteus_aculeatus',
	'LepisosteusOculatus' => 'Lepisosteus_oculatus',
	'MaylandiaZebra' => 'Maylandia_zebra',
	'OreochromisNiloticus' => 'Oreochromis_niloticus',
	'OryziasLatipes' => 'Oryzias_latipes',
	'ScleropagesFormosus' => 'Scleropages_formosus',
	'SeriolaDumerili' => 'Seriola_dumerili',
	'StegastesPartitus' => 'Stegastes_partitus',
	'TakifuguRubripes' => 'Takifugu_rubripes',
	'TetraodonNigroviridis' => 'Tetraodon_nigroviridis',
	'XiphophorusMaculatus' => 'Xiphophorus_maculatus'
);
$arrTaxTranslateBack = array_flip($arrTaxTranslate);
$sGeneSymbolIndexSp = "RhinogobiusFormosanus";

$sOut = "genespace.orthogroups.txt";

$h = popen("zcat $sGeneSpace", 'r');
$hO = fopen($sOut, 'w');

$arrRNAID2SpGeneID = fnLoadSpGeneIdMap($sRNAID2SpGeneIDMap);

$arrFam = array();
$arrSpecies = array();

while(false !== ($sLn = fgets($h)) ) {
						$sLn = trim($sLn);
						if ($sLn == "") continue;
						
					 $arrF = explode("\t", $sLn);
					 if ($arrF[0] == 'genome') continue;
					 $sSpecies = $arrF[0];
					 $sOrthoFam = $arrF[13];
					 $bIsRepresent = ($arrF[11]=='TRUE');

					 $sGeneID = $arrF[1];

					 $sSpecies = (array_key_exists($sSpecies, $arrTaxTranslate) )?  $arrTaxTranslate[$sSpecies] : $sSpecies;
					 if (!$bIsRepresent) {
					 		continue;
					 }

					 $arrSpecies[$sSpecies] = true;
					 if (!array_key_exists($sOrthoFam, $arrFam) ) {
					 		$arrFam[$sOrthoFam] = array();
					 }

					 if (!array_key_exists($sSpecies, $arrFam[$sOrthoFam])) {
					 		$arrFam[$sOrthoFam][$sSpecies] = $sGeneID;
					 } else {
					 	  echo("Warning: duplicated gene in orthogroup $sOrthoFam in $sSpecies: $sGeneID ".$arrFam[$sOrthoFam][$sSpecies]."\n");
					 }
} 

$arrSpecies = array_keys($arrSpecies);

$arrID2Symbol = fnLoadGeneSymbolMap($sGeneSymbolMap , $sGeneSymbolIndexSp);
//print_r($arrID2Symbol );
//die();

fwrite($hO, "OrthoID\tGeneSymbol\tGeneName\t");
fwrite($hO, implode("\t", $arrSpecies)."\n" );

$sGeneSymbolIndexSp = (array_key_exists($sGeneSymbolIndexSp, $arrTaxTranslate) )?  $arrTaxTranslate[$sGeneSymbolIndexSp] : $sGeneSymbolIndexSp;

foreach($arrFam as $sOrthoFam => $arrSpp) {
					$sOrthoFam = "Group_".$sOrthoFam."_1";
					 if (!array_key_exists($sGeneSymbolIndexSp, $arrSpp)) {
					 		continue; //if the focal species is missing then ignore
					 }
					 if (!array_key_exists($arrSpp[$sGeneSymbolIndexSp],$arrID2Symbol )) {
					 		fwrite($hO, "$sOrthoFam\tunknown\tunknown\t");
					 } else {
					 	 fwrite($hO, "$sOrthoFam\t".implode("\t", $arrID2Symbol[$arrSpp[$sGeneSymbolIndexSp]] ) . "\t");
					 }
					 $arrSpIDs = array();
					 foreach($arrSpecies as $sSp) {
					 				 $sSpGeneId  = $sRNAID = array_key_exists($sSp, $arrSpp)? str_replace('.','_', $arrSpp[$sSp]):'NA';	
									 if (array_key_exists($sSp, $arrRNAID2SpGeneID ) &&  array_key_exists($sRNAID, $arrRNAID2SpGeneID[$sSp] ) ) {
									 		$sSpGeneId = $sSp. "|". $arrRNAID2SpGeneID[$sSp][$sRNAID] . "|".$sRNAID;
									 }
									  
					 				 $arrSpIDs[] = $sSpGeneId;					 
					 }

					 fwrite($hO, implode("\t", $arrSpIDs)."\n" );
}

function fnLoadGeneSymbolMap($sGeneSymbolMap , $sGeneSymbolIndexSp) {
				 $h = fopen($sGeneSymbolMap, 'r');
				 $nIdx = -1;
				 $arrRet = array();
				 while(false !== ($sLn = fgets($h)) ) {
						$sLn = trim($sLn);
						if ($sLn == "") continue;
						$arrF = explode("\t", $sLn);
						if ($arrF[0] =='HOG') {
							 $arrFlip = array_flip($arrF);
							 if (!array_key_exists($sGeneSymbolIndexSp , $arrFlip)) {
							 		die("Index species $sGeneSymbolIndexSp not found in $sGeneSymbolMap\n");
							 }
							 $nIdx = $arrFlip[$sGeneSymbolIndexSp];
							 continue;
						}

						$sSymbol = $arrF[3];
						$sGeneName = $arrF[4];
						$arrIDs = explode(',', $arrF[$nIdx]);
						foreach($arrIDs as $sID) {
									 $arrRet[trim($sID)] = array($sSymbol, $sGeneName);
						}
				}

				return $arrRet;
}

function fnLoadSpGeneIdMap($sRNAID2SpGeneIDMap) {
				 $h = popen("grep '>' $sRNAID2SpGeneIDMap", 'r');
				 $arrRet = array();
				 while(false !== ($sLn = fgets($h)) ) {
						$sLn = trim($sLn);
						if ($sLn == "") continue;

						$sLn = substr($sLn, 1);
						list($sSp, $sSpGeneId, $sRNAID) = explode("|", $sLn);
						if (!array_key_exists($sSp, $arrRet)) {
							 $arrRet[$sSp] = array();
						}

						$arrRet[$sSp][$sRNAID] = $sSpGeneId;
				 }

				 return $arrRet;
}
?>
