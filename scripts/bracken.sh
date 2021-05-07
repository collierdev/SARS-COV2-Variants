source activate kraken2

for i in cleaned_reads/*_report.txt
do
  SAMPLE=$(echo ${i} | sed "s/_report\.txt//")
  echo ${SAMPLE}_report.txt

bracken -r 150 -d /gpfs1/cl/mmg232/MMG232/References/kraken2_viralDB/ -i ${SAMPLE}_report.txt -o ${SAMPLE}_bracken-out.tsv -w ${SAMPLE}_report_bracken_species.txt 
done

conda deactivate
