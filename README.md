# Type 1 Diabetes Three Country Cohort Analysis

This purpose of this project was to analyze the viral taxonomy among a sample set from three countries via the Diabimmune Database.

## Available Scripts

In the scripts directory, you can run:

### `run download.sh`
Downloads required files.\
Ideally run from prefered download directory.

### `run fastqc.sh`
Runs fastQC quality control report on raw reads.\
Ideally run in folder with downloaded raw reads.

### `run multiqc.sh`
Runs multiQC generating a aggreagte report of samples.\

### `run kraken.sh`
Runs kraken analysis on sample reads for taxonomic classification.\

### `run bracken.sh`
Runs bracken analysis on kraken output, to determine abundance at taxonomic levels.\
Ideally run from 'cleaned_reads' folder.

### `run krona.sh`
Runs Krona to create interactive taxonomic chart from analysis.\

## Output

R Analysis on Data can be found within the 'R_Analysis' folder.

### Charts and Graphs

Ouput chart and graphs can be found within the 'Images' subfolder of 'R_Analysis'.
