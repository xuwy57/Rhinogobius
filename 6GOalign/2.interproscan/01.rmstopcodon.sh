awk '{if($0 ~ /^>/) {print $0} else {gsub(/\*$/, ""); print}}' Rform.prot.fa > Rform.prot__no_stop.fasta
