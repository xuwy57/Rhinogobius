<?php
ini_set("memory_limit", -1);
//Take a ref genome, and two pseudogenome haplotypes, output a gvcf file
$sRef = $argv[1]; // "/public4/courses/ec3121/shareddata/Camellia_Sect_Chrysantha/refgenomes/Camellia_tunghinensis_hapbetter_purged/ctb.softmasked.fa";
$sHap1 = $argv[2]; // "../../whole_genome_aln/ctb.vs.Camellia_chekiangoleosa/pseudogenome.ctb.Camellia_chekiangoleosa.fa";
$sHap2 = $argv[3]; // "../../whole_genome_aln/ctb.vs.Camellia_chekiangoleosa/pseudogenome.ctb.Camellia_chekiangoleosa.fa";
$sID = $argv[4]; // "ctb.Camellia_chekiangoleosa";
$sOut = $argv[5]; // "ctb.Camellia_chekiangoleosa.g.vcf.gz";
$sVCF1 = $argv[6]; //public4/courses/ec3121/shareddata/Camellia_Sect_Chrysantha/whole_genome_aln_cth/cth.vs.Camellia_chekiangoleosa/out.cth.Camellia_chekiangoleosa.vcf.gz //get the original vcf back
$sVCF2 = $argv[7]; //public4/courses/ec3121/shareddata/Camellia_Sect_Chrysantha/whole_genome_aln_cth/cth.vs.Camellia_chekiangoleosa/out.cth.Camellia_chekiangoleosa.vcf.gz //get the original vcf back
$pO = popen("bgzip -c > $sOut", 'w');

//load ref genome
list($arrRef, $arrRefLen) = fnReadFasta($sRef);
list($arrHap1, $arrHap1Len) = fnReadFasta($sHap1);
list($arrHap2, $arrHap2Len) = fnReadFasta($sHap2);
$arrVarOnlyVCF1 = fnReadVCF($sVCF1); //this is a variant only VCF, this gets all the INDELs back
$arrVarOnlyVCF2 = fnReadVCF($sVCF2); //this is a variant only VCF, this gets all the INDELs back

if ($arrRefLen === $arrHap1Len && $arrHap1Len === $arrHap2Len) {
  echo("Ref and hap sizes ok\n");
} else {
  die("Ref and/or provided haps have different sizes, please check!\n");
}

//write vcf header lines:
fwrite($pO, '##fileformat=VCFv4.2
##ALT=<ID=NON_REF,Description="Represents any possible alternative allele not already represented at this location by REF and ALT">
##FILTER=<ID=LowQual,Description="Low quality">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another; will always be heterozygous and is not intended to describe called alleles">
##FORMAT=<ID=PID,Number=1,Type=String,Description="Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phasing set (typically the position of the first variant in the set)">
##FORMAT=<ID=RGQ,Number=1,Type=Integer,Description="Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)">
##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fisher\'s Exact Test to detect strand bias.">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=ExcessHet,Number=1,Type=Float,Description="Phred-scaled p-value for exact test of excess heterozygosity">
##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher\'s exact test to detect strand bias">
##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation">
##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
##INFO=<ID=RAW_MQandDP,Number=2,Type=Integer,Description="Raw data (sum of squared MQ and total depth) for improved RMS Mapping Quality calculation. Incompatible with deprecated RAW_MQ formulation.">
##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
##INFO=<ID=SOR,Number=1,Type=Float,Description="Symmetric Odds Ratio of 2x2 contingency table to detect strand bias">
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
    $nOneBasePos = $i + 1;
    //perform sanity checks
    $bSNPsite = ($sRefBase != 'N' && $sHap1Base != 'N' && $sHap2Base != 'N' && ($sRefBase!=$sHap1Base || $sRefBase!=$sHap2Base));
    if ( $bSNPsite && array_key_exists($sContig,$arrVarOnlyVCF1) && array_key_exists($nOneBasePos, $arrVarOnlyVCF1[$sContig])) {
      list($sVCFRef, $sVCFAlt) = $arrVarOnlyVCF1[$sContig][$nOneBasePos];
      if ($sVCFRef!=$sRefBase || $sVCFAlt!=$sHap1Base) {
        die("Fatal: VCF1 mistaches with pseudogenome hap1\n");
      }
    }
    
    if ( $bSNPsite && array_key_exists($sContig,$arrVarOnlyVCF2) && array_key_exists($nOneBasePos, $arrVarOnlyVCF2[$sContig])) {
      list($sVCFRef, $sVCFAlt) = $arrVarOnlyVCF2[$sContig][$nOneBasePos];
      if ($sVCFRef!=$sRefBase || $sVCFAlt!=$sHap2Base) {
        die("Fatal: VCF2 mistaches with pseudogenome hap2\n");
      }
    }
    
    
    
    //check against the original VCF to recover the INDEL variants, which were excluded in pseudogenomes
    if ($sHap1Base == 'N') {      
      if (array_key_exists($sContig,$arrVarOnlyVCF1) && array_key_exists($nOneBasePos, $arrVarOnlyVCF1[$sContig])) {
        list($sRefBase, $sHap1Base) = $arrVarOnlyVCF1[$sContig][$nOneBasePos]; //substitute the ref and alt bases (indels)
      }
    }
    
    if ($sHap2Base == 'N') {
      if (array_key_exists($sContig,$arrVarOnlyVCF2) && array_key_exists($nOneBasePos, $arrVarOnlyVCF2[$sContig])) {
        list($sRefBase, $sHap2Base) = $arrVarOnlyVCF2[$sContig][$nOneBasePos]; //substitute the ref and alt bases (indels)
      }
    }
    
    
    
    if ($sRefBase == 'N' || $sHap1Base == 'N' || $sHap2Base == 'N' ) {
      //fnWriteHomBlock();
      continue;
    }
    
    if ($sHap1Base != $sHap2Base && $sRefBase != $sHap1Base && $sRefBase != $sHap2Base ) {
      //fnWriteHomBlock();
      continue; //multiallelic, exclude
    }

    if ($sHomBlockPrevContig != "" && $sHomBlockPrevContig != $sContig) {
      //fnWriteHomBlock();
    }

    if ($sHomBlockPrevContig ==  $sContig && $i == $nLenMinus1) {
      //fnWriteHomBlock();
    }

    
    $sAltBase = '.';
    if ($sRefBase != $sHap1Base ) {
      $sAltBase = $sHap1Base;
    }
    
    if ($sRefBase != $sHap2Base ) {
      $sAltBase = $sHap2Base;
    }
    
    //NC_037590.1     767     .       A       .       .       .       .       GT:AD:DP:RGQ    0/0:0:0:0
    //NC_037590.1     4410    .       A       G       141.64  .       AC=1;AF=0.500;AN=2;BaseQRankSum=0.791;DP=16;ExcessHet=0.0000;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=50.83;MQRankSum=-3.040e+00;QD=8.85;ReadPosRankSum=0.815;SOR=0.527  GT:AD:DP:GQ:PL  0/1:10,6:16:99:149,0,259

    $arrOut = array($sContig, $i+1, '.', $sRefBase, ".", 30, '.', '.');
    $arrInfo = array();
    $arrInfo['GT'] = ($sHap1Base != $sHap2Base)? "0/1" : (($sRefBase == $sHap1Base)? "0/0":"1/1" );
    if ($arrInfo['GT']=="0/0") {
      $arrInfo['RGQ'] = 30;
    } else {
      $arrInfo['GQ'] = 30;
    }
    $arrInfo['DP'] = 60;
    $arrOut[4] = ($arrInfo['GT'] !='0/0')? "$sAltBase":".";
    $arrInfo['AD'] = ($arrInfo['GT']=='0/1')? '30,30': (($arrInfo['GT']=='0/0')? '60':'0,60');  
    if ($arrInfo['GT']!='0/0') {
	    $arrInfo['PL'] = ($arrInfo['GT']=='0/1')? '99,0,99': (($arrInfo['GT']=='0/0')? '0,99,99':'99,99,0');
    }
    $arrOut[8] = implode(":", array_keys($arrInfo));
    $arrOut[9] = implode(":", $arrInfo);

    /*if ($arrInfo['GT'] == "0/0") {
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
    } */

    $arrOut[8] = implode(":", array_keys($arrInfo));
    $arrOut[9] = implode(":", $arrInfo);

    fwrite($pO, implode("\t", $arrOut)."\n");
  }
}

fclose($pO);   // 改用 fclose，更可靠

echo "Re-compressing to ensure BGZF integrity and EOF marker...\n";
system("bgzip -@ 4 -r $sOut");   // 重新压缩一次，强制写入完整 EOF marker


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

function fnReadVCF($sVCF) {
  //Reads in the original Variant only VCF file
  echo("Loading variant-only VCFs $sVCF to recover INDELs from...\n");
  $p = popen("zcat -f $sVCF",'r');
  $arrRet = array();
  while(false !==($sLn = fgets($p))) {
    $arrF = explode("\t", trim($sLn));
    if ($arrF[0][0]=="#") continue;
    list($sChr, $nPos, $sID, $sRef, $sAlt, $nQual, $sFilter, $sInfo, $sFormat, $sSample) = $arrF;
    if (!array_key_exists($sChr, $arrRet)) {
      $arrRet[$sChr] = array();
    }
    if ($sSample!="1/1") {
      die("VCF file format unexpected, it must be produced by the minimap-pafcaller pipeline! $sVCF\nViolating line:\n$sLn\n");
    }
    $arrRet[$sChr][intval($nPos)] = array($sRef, $sAlt);
  }
  
  return $arrRet;
}
?>
