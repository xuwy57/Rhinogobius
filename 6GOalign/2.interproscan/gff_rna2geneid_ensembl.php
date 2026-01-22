<?php
// 输入文件
$sGOFile = $argv[1];  // GO 文件路径
$sOut = $argv[2];     // 输出文件路径

// GO类别映射
$arrTypeMap = array('biological_process' => 'BP', 'molecular_function' => 'MF', 'cellular_component' => 'CC');

// 打开 GO 文件进行解析
$hGO = fopen($sGOFile, 'r');
$hO = fopen($sOut, 'w');

// 读取并解析每一行
while (false !== ($sLn = fgets($hGO))) {
    $sLn = trim($sLn);
    if ($sLn == '') continue;
    
    $arrF = explode("\t", $sLn);
    if (count($arrF) < 8) continue; // 忽略无效行

    // 提取信息
    $sGene = $arrF[0];  // Gene stable ID
    $sRNA = "transcript_" . $arrF[2];   // Transcript stable ID
    $sGO = $arrF[4];    // GO term accession
    $sType = $arrTypeMap[$arrF[8]];  // GO category (P/F/C)

    // 输出到文件
    fwrite($hO, "$sGene\t$sRNA\t$sType\t$sGO\n");
}

// 关闭文件句柄
fclose($hGO);
fclose($hO);
?>

