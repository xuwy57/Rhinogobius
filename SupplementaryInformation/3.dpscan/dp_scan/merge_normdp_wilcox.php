<?php
$sOut = "normdp.wilcox.gz";
$nTotalParts = 100;

$hO = popen("gzip -c > $sOut", 'w');

$arrIn = array();
$arrEOF = array();
for($i=0;$i<$nTotalParts;$i++) {
	$sIn = "normdpparts/normdp.wilcox.$i.of.$nTotalParts.gz";
	if (!file_exists($sIn)) {
		die("$sIn not found!\n");
	}
	$arrIn[$i] = popen("zcat -f $sIn", 'r');
	$arrEOF[$i] = 0;
}

$nLn = -1;
$bAnyEOF = false;
while(true) {
	$nLn++;
	for($i=0;$i<$nTotalParts;$i++)  {
		$sLn = fgets($arrIn[$i]);
		if ($sLn===false) {
			$arrEOF[$i] = 1;
			$bAnyEOF = true;
			continue;
		}

		$sLn = trim($sLn);
		if ($sLn == "") continue;
		if ($i!=0 && $nLn == 0) {
			continue;
		}
		fwrite($hO, $sLn."\n");
	}

	if ($bAnyEOF && array_sum($arrEOF)==$nTotalParts) {
		break;
	}
}

echo("finished\n");
?>
