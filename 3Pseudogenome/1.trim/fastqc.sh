for f in ~/Rhinogobius_spp/DNA/01_rawdata/*.gz
do
    echo $f
    fastqc -t 8 -o ./ $f
done
