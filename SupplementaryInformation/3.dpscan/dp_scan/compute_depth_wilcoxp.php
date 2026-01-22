<?php
$nParts = intval($argv[1]); //1;
$nThisPart = intval($argv[2]); //0;
$sIn = "all.normdp.neighbor.compressed.gz";
$sOut = "normdpparts/normdp.wilcox.$nThisPart.of.$nParts.gz";
$sSpGroups = "species_assignments.txt";
$arrGroup1 = array("1");
$arrGroup2 = array("2", "5", "7"); //species in all other assignments will be ignored

exec("mkdir -p normdpparts");
$arrSpAssgn = fnReadAssignments($sSpGroups);

//print_r($arrSpAssgn);

$arrTwoGroups = fnDivideIntoTwoGroups($arrSpAssgn , $arrGroup1, $arrGroup2);

//print_r($arrTwoGroups);
$hWilcoxEngine = null;
$arrPipes = array();
fnInitWilcoxEngine();

$arrCols2Group = array(); //map column index to group
$arrGroup2Col = array(1 => array(), 2=> array());

$h = popen("zcat -f $sIn", 'r');
$hO = popen("gzip -c > $sOut", 'w');

$nLn = -1;
while(false !==($sLn = fgets($h))) {
	$sLn = trim($sLn);
	if ($sLn == "") continue;
	$arrF = explode("\t", $sLn);
	$bHeaderWritten = false;
	if ($arrF[1]=='START') {
		//process header
		$arrHeader = array_flip($arrF);
		foreach($arrHeader as $sSp => $nCol) {
			if ($sSp == 'CHROM' || $sSp == 'START' || $sSp == 'END') {
				continue;
			}

			if (!array_key_exists($sSp, $arrTwoGroups)) {
				echo("Excluded $sSp\n");
				continue;
			}
			$nGroup = $arrTwoGroups[$sSp];
			$arrCols2Group[$nCol] = $nGroup;
			$arrGroup2Col[$nGroup][] = $nCol;
		}

		if (count($arrGroup2Col[1])<2) {
			die("Too few species in group 1\n");
		}

                if (count($arrGroup2Col[2])<2) {
                        die("Too few species in group 2\n");
                }
		$$bHeaderWritten=true;
		fwrite($hO, "CHROM\tSTART\tEND\tGroup1_mean\tGroup2_mean\tGroup1_median\tGroup2_median\tWilcoxP\n");
		continue;
	}

	$nLn++;

	if ($nLn % $nParts != $nThisPart) {
		continue;
	}
	$arrGrp1Nums = array();
	$arrGrp2Nums = array();

	foreach($arrCols2Group as $nCol => $nGroup) {
		if ($nGroup == 1) {
			$arrGrp1Nums[] = $arrF[$nCol];
		} else if ($nGroup ==2) {
			$arrGrp2Nums[] = $arrF[$nCol];
		}
	}

	$arrWilcoxRet = fnWilcox($arrGrp1Nums, $arrGrp2Nums);
	fwrite($hO, "$arrF[0]\t$arrF[1]\t$arrF[2]\t".implode("\t",$arrWilcoxRet)."\n");
}

fwrite($arrPipes[0], "END\n");
proc_close($hWilcoxEngine);


function fnInitWilcoxEngine() {
	global $hWilcoxEngine, $arrPipes;

	$descriptorspec = array(
   		0 => array("pipe", "r"),  // stdin is a pipe that the child will read from
   		1 => array("pipe", "w"),  // stdout is a pipe that the child will write to
   		2 => array("file", "/dev/null", "a") // stderr is a file to write to
	);

	$hWilcoxEngine = proc_open('Rscript wilcox_engine.R', $descriptorspec, $arrPipes);
	if (!is_resource($hWilcoxEngine)) {
		die("Failed to initiate R engine\n");
	}
}

function fnWilcox($arr1, $arr2) {
	global $hWilcoxEngine, $arrPipes;

        if (!is_resource($hWilcoxEngine)) {
                die("R engine failed\n");
        }

	fwrite($arrPipes[0], implode("\t", $arr1)."\n");
	fwrite($arrPipes[0], implode("\t", $arr2)."\n");

	$sRetVal = fgets($arrPipes[1]);
	return explode("\t", $sRetVal);
}

function fnReadAssignments($sF) {
	$h = fopen($sF, 'r');
	$arrRet = array();
	while(false !== ($sLn = fgets($h))) {
		$sLn = trim($sLn);
		$arrF = explode("\t", $sLn);
		if (array_key_exists($arrF[0], $arrRet)) {
			die("$arrF[0] already defined previously, please check $sF \n");
		}

		$arrRet[$arrF[0]] = $arrF[1];
	}

	return $arrRet;
}

function fnDivideIntoTwoGroups($arrSpAssgn , $arrGroup1, $arrGroup2) {
	$arrGroup1 = array_flip($arrGroup1);
	$arrGroup2 = array_flip($arrGroup2);

	$arrTwoGroups = array();
	foreach($arrSpAssgn as $sSp => $sAssn) {
		if (array_key_exists($sAssn, $arrGroup1)) {
			$arrTwoGroups[$sSp] = 1;
			echo("$sSp: group 1\n");
		} else if (array_key_exists($sAssn, $arrGroup2)) {
			$arrTwoGroups[$sSp] = 2;
			echo("$sSp: group 2\n");
		} else {
			echo("$sSp excluded from analysis\n");
		}
	}

	return $arrTwoGroups;
}
?>
