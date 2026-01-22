<?php
$sRunDir = "rundir";
/*
$sScfPrefix = "bahahascf_"; //to be removed from scaffold names, leaving only numbers
$sSp = "Bahaha_taipingensis";
$sVersion = "1.0";
$sLongestIsoformProt = "/data/projects/rcui/bahaha_assembly/release1.0/annotations/btp.longest_isoform.prot.fa";
 */
/*
$sScfPrefix = ""; //to be removed from scaffold names, leaving only numbers
$sSp = "Nothobranchius_furzeri";
$sVersion = "1";
$sLongestIsoformProt = "~/APL/compare/ref_data/Nothobranchius_furzeri/Nothobranchius_furzeri.longest_isoform.prot.fa";
*/

$sScfPrefix = ""; //to be removed from scaffold names, leaving only numbers
$sSp = "Rhinogobius_formosanus";
$sVersion = "1";
$sLongestIsoformProt = "/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/ref/ref/Rform.longest_isoform.prot.fa";

/*
$sSp_Prot = ($argv);
$sScfPrefix = ""; //to be removed from scaffold names, leaving only numbers
$sSp = $sSp_Prot[1];
$sVersion = "1";
$sLongestIsoformProt = $sSp_Prot[2];
*/
/*
echo $sSp;
echo "\n";
echo $sLongestIsoformProt;
echo "\n"
 */

$sWD = "$sRunDir/rawGenomes/$sSp/$sVersion/annotation";
exec("mkdir -p $sWD");

$hGFF = false;
if ($sScfPrefix !='') {
	$hGFF = popen("sort -k1,1n -k4,4n | gzip -c > $sWD/$sSp.gene.gff.gz", 'w');
} else {
	$hGFF = popen("sort -k1,1 -k4,4n | gzip -c > $sWD/$sSp.gene.gff.gz", 'w');
}

$hFas = popen("gzip -c > $sWD/$sSp.pep.fa.gz", 'w');

$hIn = popen("zcat -f $sLongestIsoformProt", 'r');

while(false!== ($sLn = fgets($hIn) )) {
	$sLn = trim($sLn);
	if ($sLn == '') {
		continue;
	}

	if ($sLn[0] == '>') {
		$arrF = explode(' ', substr($sLn, 1));
		$sSeqName = $arrF[0];
		$sSeqName = preg_replace('/[:|]/', '_', $sSeqName);
		$sPos = $arrF[count($arrF)-1];
		preg_match('/([^:]+):([^-]+)-([^(]+)\\((\\S)\\)/', $sPos, $arrM);
		if (count($arrM) != 5) {
			die("Failed to parse header: $sLn\n");
		}

		list($sChr, $nStart, $nEnd, $sOrient) = array_slice($arrM, 1);
		$sChr = str_replace($sScfPrefix, '', $sChr);
		fwrite($hGFF, "$sChr\tfunannotate\tgene\t$nStart\t$nEnd\t\t$sOrient\t\tlocus=$sSeqName\n");
		fwrite($hFas, ">$sSeqName\n");
	} else {
		$sLn = str_replace('*', '', $sLn);
		fwrite($hFas, $sLn."\n");
	}

	
}
?>
