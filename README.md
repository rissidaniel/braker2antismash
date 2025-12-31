Scripts to convert Braker3 output for antismash 8

The preparation script does:
1) Create an “all CDS only” GFF3
2) Filtering the FASTA to contigs that have at least one CDS
3) Remove duplicates CDS (GeneMark), keeping AUGUSTUS
4) Convert to GenBank as input to antismash
