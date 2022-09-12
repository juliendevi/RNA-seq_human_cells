# Quality control and Kallisto alignment RNA-seq samples - SnakeMake script

## What does it do? 

This script uses FastQC v.0.11.9 (Andrews, 2010![image](https://user-images.githubusercontent.com/74775470/189668753-3fc665a5-b2eb-403a-a886-24f6fc533739.png)
) to check the quality of the sequencing. Then it uses Kallisto 0.44.0 (Bray et al., 2016)![image](https://user-images.githubusercontent.com/74775470/189668839-e5a10dc8-8a1a-457c-8123-ae51a5507c51.png)
to index the reference trasncriptome and align each sample to the reference (last GRCh38 .fasta file from NCBI has been use here). Finally all the results are summarised in a .html file using MultiQC 1.12 (Ewels et al., 2016)![image](https://user-images.githubusercontent.com/74775470/189669067-c13866fc-7ad4-451c-bd26-e9f4f47631a4.png). 

## How to run it? 
