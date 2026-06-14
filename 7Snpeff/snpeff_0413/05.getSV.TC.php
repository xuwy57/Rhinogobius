<?php

$sSnpEffVcf = "ann.final_merged.vcf.gz";
$sCloseOutgroup = "Tridentiger_barbatus";
$sFarOutgroup = "Chaeturichtys_stigmatias";
$sOut = "SV.polarizedTC.out.txt";

$hO = fopen($sOut, 'w');

//print_r($arrConsurf);

$hVCF = popen("zcat -f $sSnpEffVcf", 'r');

$nCloseOutgroupCol = 0;
$nFarOutgroupCol = 0;
$arrHaps = array();
$arrCovMask = array();
while(false !== ($sLn = fgets($hVCF) )) {
	$sLn = trim($sLn);
	$arrF = explode("\t", $sLn);
	if ($arrF[0] == "#CHROM") {
		//read in header line
		$arrRev = array_flip($arrF);
		if (!array_key_exists($sCloseOutgroup, $arrRev)) {
			die("Outgroup sp. $sCloseOutgroup not found in vcf header!\n");
		} else {
			$nCloseOutgroupCol = $arrRev[$sCloseOutgroup];
			unset($arrRev[$sCloseOutgroup]);
		}

		if (!array_key_exists($sFarOutgroup, $arrRev)) {
			die("Outgroup sp. $sFarOutgroup not found in vcf header!\n");
		} else {
			$nFarOutgroupCol =  $arrRev[$sFarOutgroup];
			unset($arrRev[$sFarOutgroup]);
		}

		$arrHaps = array_slice($arrRev, 9);
    
    // ==================== 修改这里的表头生成 ====================
    $header = "#CHROM\tPOS\tANC_ALLELE\tDER_ALLELE\tRNA_ID\tANNOTATION\tIMPACT\tPROTEIN_VARIANT";
    foreach($arrHaps as $sSample => $nCol) {
        $header .= "\t" . $sSample;  // 输出样本ID
    }
    fwrite($hO, $header . "\n");
    // ============================================================
    
		continue;

	} 
	if ($sLn[0] == '#') {
		continue;
	}

  
  // ==================== 修改部分 ====================
  // 如果两个外群中任意一个没有数据（返回 -1），就跳过这一行
  $nCloseGeno = fnParseGenotype($arrF[$nCloseOutgroupCol]);
  $nFarGeno   = fnParseGenotype($arrF[$nFarOutgroupCol]);
  
  if ($nCloseGeno == -1 || $nFarGeno == -1) {
      continue;   // 任意一个外群没有数据 → 跳过
  }
  
  if ($nCloseGeno != $nFarGeno) {
      continue;   // 两个外群都有数据但基因型不同 → 无法确定祖先，也跳过
  }
  // ================================================


	$arrAnn = fnParseSnpEffAnn($arrF[7]);
	if (count($arrAnn) == 0) {
		continue;
	}

	
	//found overlapped consurf score and variants
	//check to make sure that both outgroups are the same genotype
	$nAncestralState = fnParseGenotype($arrF[$nCloseOutgroupCol]);
	$sAncAllele = $arrF[3];
	$sDerAllele = $arrF[4];
	if ($nAncestralState == 2) {
		$sAncAllele = $arrF[4];
		$sDerAllele = $arrF[3];
	}

	$arrPolarizedSamples = array(); 
  $bAllNA = true;                    // 新增：标记是否全部为 NA

	foreach($arrHaps as $sSample => $nCol) {
		$nGeno = fnParseGenotype($arrF[$nCol]);
		if ($nGeno==-1 ) {
			$arrPolarizedSamples[] = 'NA';
		} else if ($nGeno == 1) { //het site
			$arrPolarizedSamples[] = 1;
      $bAllNA = false;
		} else if ($nGeno == $nAncestralState) {
			$arrPolarizedSamples[] = 0;
      $bAllNA = false;
		} else {
			$arrPolarizedSamples[] = 2;
      $bAllNA = false;
		}
	}
  
  if ($bAllNA) {
        continue;
    }

	fwrite($hO, $arrF[0]."\t".$arrF[1]."\t".$sAncAllele."\t".$sDerAllele."\t".implode("\t", $arrAnn[0]) ."\t".implode("\t", $arrPolarizedSamples)."\n");

}

function fnParseGenotype($s) {
	preg_match(';([01])[|/]([01]);', $s, $arrM);
	if (count($arrM)!=3) {
		return -1;
	}

	return intval($arrM[1]) + intval($arrM[2]);
}

function fnLoadMask($arrNames) {
	$arrRet = array();
	foreach($arrNames as $sName) {
		$sF = "pairwise.simp.ref.".$sName.".txt";
		echo("Loading mask $sF ...\n");
		if (!file_exists($sF)) {
			echo("Warning: Mask $sF not found. Will not apply\n");
			continue;
		}

		$h = fopen($sF, 'r');
		while(false !== ($sLn = fgets($h))) {
			$sLn = trim($sLn);
			$arrF = explode("\t", $sLn);
			$sScf = $arrF[4];
			$nStart = $arrF[5];
			$nEnd = $arrF[6];
			if (!array_key_exists($sScf, $arrRet)) {
				$arrRet[$sScf] = array();
			}

			$arrRet[$sScf][] = array($nStart, $nEnd );
		}

		foreach(array_keys($arrRet) as $sScf) {
			usort($arrRet[$sScf], 'fnSortCoord');
		}
	}


	return $arrRet;
}

function fnSortCoord($arrCoord1, $arrCoord2) {
	$a = min($arrCoord1);
	$b = min($arrCoord2);
	if ($a == $b) {
        return 0;
    }
    return ($a < $b) ? -1 : 1;
}

function fnInRegion( $sScf, $nPos) {
	global $arrCovMask;
	if (!array_key_exists($sScf, $arrCovMask)) {
		return false;
	}

	foreach($arrCovMask[$sScf] as $arrCoords) {
		if ($nPos > $arrCoords[0] && $nPos <= $arrCoords[1]) {
			return true;
		}

	}

	return false;
}

function fnParseSnpEffAnn($s) {
	$arrF1 = explode(";", $s);
	$arrMissenseVars = array();
	foreach($arrF1 as $sField) {
		list($sKey, $sVal) = explode('=', $sField);
		if ($sKey != 'ANN') {
			continue;
		}

		$arrAnns = explode(',', $sVal);
		foreach($arrAnns as $sAnnStr) {
			$arrAnnInfo = explode('|', $sAnnStr);
			$sRNAID = $arrAnnInfo[6];
			$sProteinVar = $arrAnnInfo[10];
			$arrMissenseVars[] = array( $sRNAID , $arrAnnInfo[1], $arrAnnInfo[2], $sProteinVar );
		}
	}

	return $arrMissenseVars;
}
function fnLoadConsurf($s) {
	$h = popen("zcat -f $s", "r");
	$arrRet = array();
	while(false !== ($sLn = fgets($h) )) {
		$sLn = trim($sLn);
		$arrF = explode("\t", $sLn);
		$sCoord = $arrF[0].":".$arrF[1];
		if (!array_key_exists($sCoord, $arrRet)) {
			$arrRet[$sCoord] = array();
		} 

		$arrRet[$sCoord][$arrF[3]] = $arrF[2];
	}

	return $arrRet;
}
?>