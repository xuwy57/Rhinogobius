<?php
//if the neighboring records are completely identical, join them into one record to speed up the computation
$sIn = "all.normdp.gz";
$sOut = "all.normdp.neighbor.compressed.gz";

$h = popen("zcat $sIn", 'r');
$hO = popen("gzip -c > $sOut",'w');

$sPrevChr = "";
$nPrevPos = -1;
$nBlockStart = -1;
$sPrevLn = "";
do {

	$sLn = fgets($h);
	$arrF = explode("\t", trim($sLn));

	if ($arrF[0]=='CHROM') {
		$arrNewHeader = array("CHROM", "START", "END") + array_slice($arrF, 2);
		
		fwrite($hO, implode("\t",$arrNewHeader)."\n");
		continue; //write header
	}
	if ($sLn === false || (count($arrF)>2 && ($arrF[0]!=$sPrevChr || $arrF[1]!=($nPrevPos+1))) ) {
		//write prev.
		if ($nPrevPos!=-1) {
			fwrite($hO, "$sPrevChr\t$nBlockStart\t$nPrevPos\t$sPrevLn\n");
		}
		if ($sLn === false) {
			break;
		} 

		$sPrevChr = $arrF[0];
		$nPrevPos = $nBlockStart =  $arrF[1];
		$sPrevLn = implode("\t", array_slice($arrF, 2) );
	}

	$sCurrLn = implode("\t", array_slice($arrF, 2) );
	if ($sCurrLn === $sPrevLn) {
		$nPrevPos = $arrF[1];
		continue;
	} else {
		if ($nPrevPos!=-1) {
			fwrite($hO, "$sPrevChr\t$nBlockStart\t$nPrevPos\t$sPrevLn\n");
		}
                $sPrevChr = $arrF[0];
                $nPrevPos = $nBlockStart =  $arrF[1];
                $sPrevLn = implode("\t", array_slice($arrF, 2) );
	}

} while (true);

?>
