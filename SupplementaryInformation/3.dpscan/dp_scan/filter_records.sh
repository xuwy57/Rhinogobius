zcat normdp.wilcox.gz | awk '{if ($8<1e-5) print $0}' > sign.1e-5.txt
cat sign.1e-5.txt | awk '{if ($4>$5) print $0}' > sign.1e-5.grp1.larger.txt
cat sign.1e-5.txt | awk '{if ($5>$4) print $0}' > sign.1e-5.grp2.larger.txt

