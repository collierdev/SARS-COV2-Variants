source activate kraken2
cd /gpfs1/home/j/w/jwcollie/metagenomics/cleaned_reads
for i in /gpfs1/home/j/w/jwcollie/metagenomics/cleaned_reads/*.out

do
  SAMPLE=$(echo ${i} | sed "s/.out//") 

ktImportTaxonomy -q 2 -t 3 ${sample}.out -o ${sample}_krona.html
done

conda deactivate
