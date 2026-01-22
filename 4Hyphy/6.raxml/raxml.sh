#ln -s ../phylipformat/full.phy ./in.phy

#/public/software/raxml/standard-RAxML/raxmlHPC-PTHREADS-AVX -s in.phy -T 96 -f a -x 23333 -N 100 -m GTRGAMMA -n allcodons -p 23333

/data/software/standard-RAxML-8.2.12/raxmlHPC-PTHREADS-AVX2 -s in.phy -T 96 -f a -x 23333 -N 100 -m GTRGAMMA -n allcodons -p 23333


