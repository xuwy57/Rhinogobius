module load bedtools

#awk -F"\t" -v OFS='\t' '
#    !/^#/ { 
#        start = $2 - 1;
#        printf("%s\t%s\t%s", $1, start, $2);
#        for (i=3; i<=NF; i++) printf("\t%s", $i);
#        print "";
#    }
#' SV.polarizedRT.out.txt > SV.polarizedRT.out.bed
grep '^Chromosome' mask.ref.Tba.bed > mask.filtered.ref.Tba.bed
bedtools intersect -a SV.polarizedRT.out.bed -b mask.filtered.ref.Tba.bed > SV.polarizedRT.filtered.bed
