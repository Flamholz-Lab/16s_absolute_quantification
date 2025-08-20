#!/bin/bash

# create your working environment
conda install -y -c conda-forge mamba
mamba create -y -n hb-dada2-ex-wf -c conda-forge -c bioconda -c defaults \
    cutadapt r-base r-tidyverse r-vegan r-dendextend r-viridis \
    bioconductor-phyloseq bioconductor-deseq2 bioconductor-dada2 \
    bioconductor-decipher bioconductor-decontam r-biocmanager r-matrix
conda activate hb-dada2-ex-wf

# cd into your folder containing raw reads
cd /Users/albertli/Desktop/lab/Flamholz/projects/absolute_16s_quantification/data/standard_sequencing_07222025/Albert

# unzip your raw reads
dest="../unzipped"
mkdir -p "$dest"
find . -type f -name '*.gz' -print0 |
  while IFS= read -r -d '' f; do
    gunzip -c "$f" > "$dest/$(basename "${f%.gz}")"
  done

# cutadapt to trim your reads
dest="../trimmed"
mkdir -p "$dest"
cd ../unzipped
ls *_R1_001.fastq | sed 's/_R1_001\.fastq//' > samples

for sample in $(cat samples)
do

    echo "On sample: $sample"

    cutadapt -g ^CCTACGGGNGGCNGCAG \
             -o ../trimmed/${sample}_R1_trimmed.fastq \
             ${sample}_R1_001.fastq \
             >> cutadapt_primer_trimming_stats.txt 2>&1

done

for sample in $(cat samples)
do

    echo "On sample: $sample"
    
    cutadapt -g ^GACTACNNGGGTATCTAATCC \
             -o ../trimmed/${sample}_R2_trimmed.fastq \
             ${sample}_R2_001.fastq \
             >> cutadapt_primer_trimming_stats.txt 2>&1

done

