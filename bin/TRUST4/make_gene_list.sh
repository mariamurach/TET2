grep ">" IMGT+C.fa | cut -f2 -d'>' | cut -f1 -d'*' | sort | uniq > bcr_tcr_gene_name.txt
