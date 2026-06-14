<?php
ini_set("memory_limit", -1);
//Take a ref genome, and two pseudogenome haplotypes, output a gvcf file
$sRef = $argv[1]; // "/public4/courses/ec3121/shareddata/Camellia_Sect_Chrysantha/refgenomes/Camellia_tunghinensis_hapbetter_purged/ctb.softmasked.fa";
$sHap1 = $argv[2]; // "../../whole_genome_aln/ctb.vs.Camellia_chekiangoleosa/pseudogenome.ctb.Camellia_chekiangoleosa.fa";
$sHap2 = $argv[3]; // "../../whole_genome_aln/ctb.vs.Camellia_chekiangoleosa/pseudogenome.ctb.Camellia_chekiangoleosa.fa";
$sID = $argv[4]; // "ctb.Camellia_chekiangoleosa";
$sOut = $argv[5]; // "ctb.Camellia_chekiangoleosa.g.vcf.gz";
$pO = popen("bgzip -c > $sOut", 'w');

//load ref genome
list($arrRef, $arrRefLen) = fnReadFasta($sRef);
list($arrHap1, $arrHap1Len) = fnReadFasta($sHap1);
list($arrHap2, $arrHap2Len) = fnReadFasta($sHap2);

if ($arrRefLen === $arrHap1Len && $arrHap1Len === $arrHap2Len) {
  echo("Ref and hap sizes ok\n");
} else {
  die("Ref and/or provided haps have different sizes, please check!\n");
}

//write vcf header lines:
fwrite($pO, '##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=RefCall,Description="Genotyping model thinks this site is reference.">
##FILTER=<ID=LowQual,Description="Confidence in this variant being real is below calling threshold.">
##FILTER=<ID=NoCall,Description="Site has depth=0 resulting in no call.">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position (for use with symbolic alleles)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Conditional genotype quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Min Read depth in block">
##FORMAT=<ID=VAF,Number=A,Type=Float,Description="Variant allele fractions.">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer">
##PSEUDOGENOME_version=1.0
');
foreach($arrRefLen as $sContig => $nLen) {
  fwrite($pO,"##contig=<ID=$sContig,length=$nLen>\n");
} 
fwrite($pO,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sID\n");

$sHomBlockPrevContig = "";
$sHomBlockPrevCoord = -1;
$arrHomBlockOut = array();

foreach($arrRefLen as $sContig => $nLen) {
  $nLenMinus1 = $nLen -1;
  for($i=0;$i<$nLen;$i++) {
    $sRefBase = $arrRef[$sContig][$i];
    $sHap1Base = $arrHap1[$sContig][$i];
    $sHap2Base = $arrHap2[$sContig][$i];
    if ($sRefBase == 'N' || $sHap1Base == 'N' || $sHap2Base == 'N' ) {
      fnWriteHomBlock();
      continue;
    }
    
    if ($sHap1Base != $sHap2Base && $sRefBase != $sHap1Base && $sRefBase != $sHap2Base ) {
      fnWriteHomBlock();
      continue; //multiallelic, exclude
    }

    if ($sHomBlockPrevContig != "" && $sHomBlockPrevContig != $sContig) {
      fnWriteHomBlock();
    }

    if ($sHomBlockPrevContig ==  $sContig && $i == $nLenMinus1) {
      fnWriteHomBlock();
    }

    
    $sAltBase = '.';
    if ($sRefBase != $sHap1Base ) {
      $sAltBase = $sHap1Base;
    }
    
    if ($sRefBase != $sHap2Base ) {
      $sAltBase = $sHap2Base;
    }
    
    $arrOut = array($sContig, $i+1, '.', $sRefBase, "<*>", 30, '.', '.');
    $arrInfo = array();
    $arrInfo['GT'] = ($sHap1Base != $sHap2Base)? "0/1" : (($sRefBase == $sHap1Base)? "0/0":"1/1" );
    $arrInfo['GQ'] = 30;
    $arrInfo['DP'] = 60;
    $arrOut[4] = ($arrInfo['GT'] !='0/0')? "$sAltBase,<*>":"<*>";
    $arrInfo['AD'] = ($arrInfo['GT']=='0/1')? '30,30,0': (($arrInfo['GT']=='0/0')? '60,0':'0,60,0');  
    if ($arrInfo['GT']=='0/0') {
	unset($arrInfo['AD'] );
    }
    $arrInfo['PL'] = ($arrInfo['GT']=='0/1')? '99,0,99,990,990,990': (($arrInfo['GT']=='0/0')? '0,99,990':'99,99,0,990,990,990');
    $arrOut[8] = implode(":", array_keys($arrInfo));
    $arrOut[9] = implode(":", $arrInfo);

    if ($arrInfo['GT'] == "0/0") {
	//check if can be merged 
	if ($sHomBlockPrevContig == "") {
		$sHomBlockPrevContig = $sContig;
		$sHomBlockPrevCoord = $i+1;
		unset($arrInfo['DP'] );
		$arrInfo['MIN_DP'] = 60;
		$arrOut[8] = implode(":", array_keys($arrInfo));
		$arrOut[9] = implode(":", $arrInfo);
		$arrHomBlockOut = $arrOut;
		continue;
	} else {
		$sHomBlockPrevCoord = $i+1;
		continue; //merge
	}
    } else {
	      fnWriteHomBlock();
    }

    $arrOut[8] = implode(":", array_keys($arrInfo));
    $arrOut[9] = implode(":", $arrInfo);

    fwrite($pO, implode("\t", $arrOut)."\n");
  }
}

function fnWriteHomBlock() {
	global $sHomBlockPrevContig, $sHomBlockPrevCoord, $arrHomBlockOut, $pO;
	if ($sHomBlockPrevContig == "" || $sHomBlockPrevCoord ==-1) {
		return;
	}
	$arrHomBlockOut[7] = "END=$sHomBlockPrevCoord";
	fwrite($pO, implode("\t", $arrHomBlockOut)."\n");
	$sHomBlockPrevContig = "";
	$sHomBlockPrevCoord = -1;
	$arrHomBlockOut = array();

}

exec("tabix $sOut "); 

function fnParseAnn($f, $v) {
  $arrFields = explode(":", $f);
  $arrValues = explode(":", $v);
  return array_combine($arrFields, $arrValues);
}

function fnParseEND($s)      
{
	$arrS = explode('=', $s);
	if ($arrS[0] != 'END') {
		return false;
	}

	return intval($arrS[1]);
}

function fnReadFasta($s) {
echo("Loading $s ...\n");
$h = fopen($s, 'r');
$arrRet = array();
$arrLen = array();
$sSeqName = "";
$sSeq = "";
   do {
        $sLn = fgets($h);

        if ($sLn==false) {
                $arrRet[$sSeqName] = $sSeq;
                $arrLen[$sSeqName] = strlen($sSeq);
                break;
        }

        $sLn = trim($sLn);
        if (substr($sLn , 0,1) == ">") {
                if ($sSeqName!='') {
                              $arrRet[$sSeqName] = $sSeq;
                  $arrLen[$sSeqName] = strlen($sSeq);
                }
                list($sSeqName) = preg_split('/\s+/', substr($sLn , 1));
                $sSeq = "";

        } else {
                $sSeq .= strtoupper($sLn);
        }

    } while (true);

return array($arrRet , $arrLen);
}


?>
