<?php
$arrFiles = glob("normdp/*.gz");
$sOut = "all.normdp.gz";

$hO = popen("gzip -c > $sOut", 'w');
$arrIn = array();

foreach($arrFiles as $i => $sF) {
	$arrIn[$i] = popen("zcat -f $sF", 'r');
}

$nSamples = count($arrIn);
while(true) {
	$sChr = "";
	$nPos = -1;
	$bEOF = false;
	$bHeader = false;
	$arrSamplesInHeader = array();
	for($i=0;$i<$nSamples;$i++) {
		$sLn = fgets($arrIn[$i]);
		if ($sLn === false) {
			$bEOF = true;
			continue;
		}
		if ($bEOF) {
			die("Error! Files are not synchronized, some file ended earlier\n".$arrFiles[$i-1]."\t".$arrFiles[$i]."\n");
		}

		list($sThisChr, $nThisPos, $nNormDP) = explode("\t", trim($sLn));
		if ($nThisPos == 'POS') {
			$bHeader = true;
			$arrSamplesInHeader[] = $nNormDP;
			continue;
		} else {
			if ($bHeader) {
				die("Error! Header expected but not found in file ".$arrFiles[$i]."\n");
			}
			if ($nPos == -1) {
				$sChr = $sThisChr;
				$nPos = $nThisPos;
			}

			if ($sChr != $sThisChr || $nPos != $nThisPos) {
				die("Error! Files are not synchronized, file coord unmatch\n".$arrFiles[$i-1]." $sChr:$nPos\t".$arrFiles[$i]." $sThisChr:$nThisPos\n");

			}
			$arrSamplesInHeader[] = $nNormDP;
		}
	}

        if ($bEOF) {
                break;
        }


	if ($bHeader) {
		fwrite($hO, "CHROM\tPOS\t".implode("\t", $arrSamplesInHeader)."\n");
		continue;
	} else {
		fwrite($hO, "$sChr\t$nPos\t".implode("\t", $arrSamplesInHeader)."\n");
                continue;
	}

}

echo("finished\n");
?>
