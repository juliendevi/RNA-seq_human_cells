# Quality control and Kallisto alignment RNA-seq samples - SnakeMake script

## What does it do? 

This script uses FastQC v.0.11.9 (Andrews, 2010) to check the quality of the sequencing. Then it uses Kallisto 0.44.0 (Bray et al., 2016) to index the reference trasncriptome and align each sample to the reference (last GRCh38 .fasta file from NCBI has been used here). Finally all the results are summarised in a .html file using MultiQC 1.12 (Ewels et al., 2016). 

## How to run it? 

Create a data folder with a reference transcriptome (named ref_genome.fasta) and a folder samples/ with all the samples (_1.fq.gz and _2.fq.gz for per end sequencing). 
Fill the _config.yaml file. 
Run snakemake as follow to include the conda environement: 

```
snakemake --cores 28 --use_conda 
```

### References 

Andrews, S. (2010) FastQC: A Quality Control Tool for High Throughput Sequence Data.

Bray, N.L., Pimentel, H., Melsted, P. & Pachter, L. (2016) Near-optimal probabilistic RNA-seq quantification. Nature Biotechnology, 34, 525–527.

Ewels, P., Magnusson, M., Lundin, S. & Käller, M. (2016) MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32, 3047–3048.
