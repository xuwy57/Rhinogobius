for i in */pseudogenome.*.fa; do
 cat $i | grep -io 'N' | wc -l > $i.countN.txt &
done
wait
