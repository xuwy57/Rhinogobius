grep -h CV log*.out |sort -nk4  -t ' ' |tail -n +2 |perl -p -e 's/.*K=(.+?)\):.*/$1/g' |xargs -i rm QC.{}.Q QC.{}.P log{}.out

num=`grep -h CV log*.out |sort -nk4  -t ' ' |head -n 1 |perl -p -e 's/.*K=(.+?)\):.*/$1/g' `
paste <(cut -f 1 -d ' ' Xiphi.fam  ) <(cat Xiphi.${num}.Q|tr ' ' '\t')  > QC.${num}.Q.tabq
