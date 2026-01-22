sed -E "s/^((@|\+)SRR[^.]+\.[^.]+)\.(1|2)/\1/" <(zcat ../Rsimilis.paired_1.fq.gz)  > Rsimilis_rename_pass_1.fastq &
sed -E "s/^((@|\+)SRR[^.]+\.[^.]+)\.(1|2)/\1/" <(zcat ../Rsimilis.paired_2.fq.gz)  > Rsimilis_rename_pass_2.fastq &
