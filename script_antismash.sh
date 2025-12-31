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


conda activate /mnt/mountdoom/mnt/work/bioinfo/home/dar21/miniconda3/envs/antismash8

for f in *.fna;
         do SAMPLE=`basename ${f%%.fna}`;
                echo "$SAMPLE";
                antismash "$SAMPLE".antismash_input.gbk  -t fungi -c 48 --output-dir "$SAMPLE"_antismash8 --cc-mibig --cb-knownclusters --tfbs --clusterhmmer --pfam2go --smcog-trees  --verbose --logfile "${SAMPLE}_antismash8/antismash.log" ;		
done;

echo "run completed"
echo "END TIME: " `date`
           
