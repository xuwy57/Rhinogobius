<?php
//only works for NCBI Refseq gff
$sGFF = $argv[1]; //"/public4/courses/ec3121/shareddata/Camellia_Sect_Chrysantha/refgenomes/Malvales_Gossypium_hirsutum/GCF_007990345.1_Gossypium_hirsutum_v2.1_genomic.gff";
$sGAF = $argv[2]; //"/public4/courses/ec3121/shareddata/Camellia_Sect_Chrysantha/refgenomes/Malvales_Gossypium_hirsutum/GCF_007990345.1_Gossypium_hirsutum_v2.1_gene_ontology.gaf.gz";
$sOut = $argv[3]; //"out.txt";
$arrTypeMap = array('P'=>'BP', 'F' => 'MF', 'C' => 'CC'); //P (biological process), F (molecular function) or C (cellular component)

$h = popen("zcat -f $sGFF", 'r');
$arrGene2RNA = array();
while(false !== ($sLn=fgets($h)) ) {
	$sLn = trim($sLn);
	$arrF = explode("\t", $sLn);
	if (count($arrF)<9) {
		continue;
	}
	if ($arrF[2] == 'mRNA') {
		$sRNAID = "";
		$sGeneID = "";
		preg_match('/ID=([^;]+)/', $arrF[8], $arrM);
		if (count($arrM)==2) {
			$sRNAID = $arrM[1];
		}
		preg_match('/GeneID:([^;,]+)/', $arrF[8], $arrM);
                if (count($arrM)==2) {
                        $sGeneID = $arrM[1];
                }

		if ($sRNAID!='' && $sGeneID!="") {
			if (!array_key_exists($sGeneID, $arrGene2RNA)) {
				$arrGene2RNA[$sGeneID] = array();
			}

			$arrGene2RNA[$sGeneID][$sRNAID] = 1; 
		}

	}
}

//print_r($arrGene2RNA);

$hGAF = popen("zcat -f $sGAF" , 'r');
$hO = fopen($sOut, 'w');
while(false !== ($sLn=fgets($hGAF)) ) {
        $sLn = trim($sLn);
	if ($sLn == '') continue;
	if ($sLn[0] == '!') continue;
        $arrF = explode("\t", $sLn);
        if (count($arrF)<13) {
                continue;
        }

	$sGene = $arrF[1];
	$sGO = $arrF[4];
	$sType = $arrTypeMap[$arrF[8]];

	if (array_key_exists($sGene, $arrGene2RNA)) {
		foreach($arrGene2RNA[$sGene] as $sRNA => $s) {
			fwrite($hO, "$sGene\t$sRNA\t$sType\t$sGO\n");
		}
	}

}

?>
