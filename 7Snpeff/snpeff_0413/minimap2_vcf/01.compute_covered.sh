module load minimap2

paftools=paftools.js

#cd ./rform_Ggi/
#$paftools call -q 60 -L10000 -l1000 pairwise.aln.Rform.Glossogobius_giuris.paf | grep -P "R\t" | cut -f2-4 > covered.regions.Rform.Glossogobius_giuris.txt  && touch cov.done

cd ./rform_Tba/
$paftools call -q 60 -L200 -l200 pairwise.aln.Rform.Tridentiger_barbatus.paf | grep -P "R\t" | cut -f2-4 > covered.regions.Rform.Tridentiger_barbatus.txt  && touch cov.done

#cd ./rform_Cstig/
#$paftools call -q 60 -L200 -l200 pairwise.aln.Rform.Chaeturichtys_stigmatias.paf | grep -P "R\t" | cut -f2-4 > covered.regions.Rform.Chaeturichtys_stigmatias.txt  && touch cov.done
