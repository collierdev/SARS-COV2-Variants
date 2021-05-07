#!/bin/bash
#SBATCH --partition=bluemoon
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=128G
#SBATCH --job-name=krak
#SBATCH --output=%x_%j.out
#SBATCH --mail-user=jwcollie@uvm.edu
#SBATCH --mail-type=ALL

cd /gpfs1/home/j/w/jwcollie/metagenomics/raw_reads

source activate kraken2

for i in /gpfs1/home/j/w/jwcollie/metagenomics/raw_reads/*_1.fastq.gz

do
  SAMPLE=$(echo ${i} | sed "s/_pe_1.fastq.gz//") 


kraken2 --db /gpfs1/cl/mmg232/MMG232/References/kraken2_standard/ --threads 6 --output ${SAMPLE}.out --report ${SAMPLE}_report.txt --paired ${SAMPLE}_pe_1.fastq.gz ${SAMPLE}_pe_2.fastq.gz
 

done

conda deactivate
