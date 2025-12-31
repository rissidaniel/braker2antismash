#!/bin/bash

#SBATCH -c 48                      # 1 core per job (i.e., if you need 8 cores, you would have to use "-c 8")
#SBATCH -t 90:00:00                # Runtime in D-HH:MM
#SBATCH -p long                  # Partition to submit to
#SBATCH --mem=384gb                 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --mail-user=dar21@dsmz.de
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

echo "START TIME: " `date`
echo "load conda"

eval "$(conda shell.bash hook)"

conda activate /mnt/mountdoom/mnt/work/bioinfo/home/dar21/miniconda3/envs/gff2gbk


for GFF in *.gff3; do
  S="${GFF%.gff3}"
  FA="${S}.fna"
  [ -s "$FA" ] || { echo "Missing $FA, skipping $S"; continue; }

  echo "=== $S ==="

  # 1) all CDS only
  awk -F'\t' 'BEGIN{OFS="\t"; print "##gff-version 3"} $0!~/^#/ && $3=="CDS"{print}' \
    "$GFF" > "${S}.CDS.all.gff3"

  # 2) contigs with CDS + filter FASTA
  cut -f1 "${S}.CDS.all.gff3" | sort -u > "${S}.contigs_with_CDS.txt"
  awk 'BEGIN{while((getline<ARGV[1])>0) keep[$1]=1; ARGV[1]=""}
       /^>/{id=substr($1,2); p=(id in keep)}
       {if(p) print}
  ' "${S}.contigs_with_CDS.txt" "$FA" > "${S}.CDSonly_contigs.fna"

  # 3) remove GeneMark duplicates + dedup exact loci
  awk -F'\t' '
  BEGIN{OFS="\t"}
  FNR==NR{key=$1 FS $4 FS $5 FS $7; if($2!="GeneMark.hmm3") other[key]=1; next}
  {key=$1 FS $4 FS $5 FS $7;
   if($2=="GeneMark.hmm3" && other[key]) next;
   if(seen[key]++) next;
   print}
  ' "${S}.CDS.all.gff3" "${S}.CDS.all.gff3" > "${S}.CDS.noGenemarkDups.gff3"


  # convert to GBK
   python gff3_cds_to_gbk.py "${S}.CDSonly_contigs.fna" "${S}.CDS.noGenemarkDups.gff3" "${S}.antismash_input.gbk"

done


mkdir tmp_for_antismash

mv *.CDS.all.gff3 tmp_for_antismash
mv *.contigs_with_CDS.txt tmp_for_antismash

mv *.CDS.noGenemarkDups.gff3 tmp_for_antismash
mv *.CDSonly_contigs.fna tmp_for_antismash


#antismash Arachnopeziza_araneosa.antismash_input.gbk -t fungi -c 48 --output-dir Arachnopeziza_araneosa_antismash71 --cc-mibig --cb-knownclusters --clusterhmmer --pfam2go --smcog-trees --tfbs

#for f in *.fna;
         #do SAMPLE=`basename ${f%%.fna}`;
                #echo "$SAMPLE";
                #EMBLmyGFF3   "$SAMPLE".gff3  "$SAMPLE".fasta --topology linear --molecule_type "genomic DNA" --transl_table 1 -o "$SAMPLE".embl --locus_tag AVERR --species "$SAMPLE" --project_id "$SAMPLE";
		#python gff_to_genbank.py "$SAMPLE".gff3 "$SAMPLE".fasta "genomic DNA";
		#antismash "$SAMPLE".fna  -t fungi -c 48 --genefinding-gff3 "$SAMPLE".CDSonly.gff3 --output-dir "$SAMPLE"_antismash71 --cc-mibig --cb-knownclusters --tfbs --clusterhmmer --pfam2go --smcog-trees  --verbose --logfile "${SAMPLE}_antismash71/antismash.log" ;		
#done;

echo "run completed"
echo "END TIME: " `date`
           
