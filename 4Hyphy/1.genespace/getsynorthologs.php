<?php
//ver 2, fixed bug of colliding pgid
$sIn = "rundir/results/RhinogobiusFormosus_pangenomeDB.txt.gz";
//$sIn = "debugin.txt";
$sOut = "synorthos.txt";
//$sOut = "debigout.txt";

$arrGenomes = fnGetGenomes($sIn); 
$hIn = popen("zcat -f $sIn" , 'r');
//$hIn = fopen($sIn, 'r');
$hO = fopen($sOut, 'w');

$arrOrths = array();
while(false !== ($sLn = fgets($hIn)) ) {
	$sLn = trim($sLn);
	if ($sLn == '') continue;
	$arrF = explode("\t", $sLn);
	if ($arrF[0] == 'pgChr') {
		continue;
	}

	if (count($arrF) != 14) {
	//		die("Error: $sLn\n");
		continue;
	}

	$nPGID = $arrF[2];
	$sUniqID = implode("_", array_slice($arrF, 0,4) ); //"$arrF[0]"."_".$nPGID;
	if (!array_key_exists($sUniqID , $arrOrths)) {
		$arrOrths[$sUniqID] = array('pgChr' => $arrF[0], 'pgOrd' => $arrF[1], 'pgid' => $nPGID  , 'directSynOrth' => array() , 'NSOrths' => array() );
		foreach($arrGenomes as $sSample) {
			$arrOrths[$sUniqID]['directSynOrth'][$sSample] = "NA";
			$arrOrths[$sUniqID]['NSOrths'][$sSample] = array();
		}
	}

	$sGenome = $arrF[6];
	$bIsDirectSyn = ($arrF[9]=='TRUE');
	if ($bIsDirectSyn) {
		$arrOrths[$sUniqID]['directSynOrth'][$sGenome] = $arrF[13];
	}

	$bIsNSOrtho = ($arrF[12] == 'TRUE');
	if ($bIsNSOrtho) {
		$arrOrths[$sUniqID]['NSOrths'][$sGenome][$arrF[13]] = true ;
	}
}

fwrite($hO, "pgChr\tpgOrd\tpgID\t".implode("\t",$arrGenomes)."\t".implode("\t",$arrGenomes)."\n");
foreach($arrOrths as $sUniqID => $arrOr) {
	$nPGID = $arrOr['pgid'];
	fwrite($hO, $arrOr['pgChr']."\t".$arrOr['pgOrd']."\t$nPGID\t".implode("\t",$arrOr['directSynOrth'])."\t");
	$arrNSOrths = array();
	foreach($arrOr['NSOrths'] as $arrNSOr ) {
		if (count($arrNSOr)==0) {
			$arrNSOrths[] = 'NA';
		} else {
			$arrNSOrths[] = implode(';', array_keys($arrNSOr) );
		}
	}
	fwrite($hO, implode("\t",$arrNSOrths )."\n");
}


function fnGetGenomes($sIn) {

	$h = popen("zcat -f $sIn | cut -f7 | grep -v genome | sort | uniq", 'r');
	$arrRet = array();
	while(false !== ($sLn = fgets($h)) ) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		$arrRet[] = $sLn;
	}

	pclose($h);
	return $arrRet;
}
?>
