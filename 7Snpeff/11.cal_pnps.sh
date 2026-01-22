module load bcftools

bcftools query -f '%CHROM\t%POS\t%INFO/ANN\t[%SAMPLE=%GT\t]\n' ann.merged.snps.indels.vcf.gz | \
awk '
{
    for(i=4; i<=NF; i++) {
        split($i,a,"=");
        sample=a[1];
        gt=a[2];
        if(gt=="0/1" || gt=="1/1") {
            if($3 ~ /missense_variant/) pN[sample]++;
            if($3 ~ /synonymous_variant/) pS[sample]++;
        }
    }
}
END {
    for(s in pN) {
        pn = pN[s]+0;
        ps = pS[s]+0;
        if(ps==0) ratio="NA"; else ratio=pn/ps;
        print s, pn, ps, ratio;
    }
}' > all_samples_pNpS.tsv

