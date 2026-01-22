<?php
$sIn="/data/projects/zma/R.formosanus_data/10_funannotate/repeatmodeler/hic/Rform_hic-families.fa";
$sOut = "Rform_hic-te-families.fa";
$sSpecies = "Rhinogobius_formosanus";

$h = fopen($sIn, 'r');
$hO = fopen($sOut, 'w');

while(false !== ($sLn = fgets($h)) ) {
	if ($sLn[0] != '>') {
		fwrite($hO, $sLn);
		continue;
	}

	list($s1, $s2)  =  explode(' ' , $sLn);

	fwrite($hO, "$s1 @$sSpecies\n");
}
?>
