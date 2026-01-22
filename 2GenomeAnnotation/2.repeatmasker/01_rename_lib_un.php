<?php
$sIn="/data/projects/zma/R.formosanus_data/10_funannotate/repeatmodeler/hic/Rform-families.fa";
$sOut = "Rform-te-families.fa";
$sSpecies = "Rhinogobius_form";

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
