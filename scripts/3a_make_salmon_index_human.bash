#!/bin/bash

# Load Salmon module
module load salmon/1.10.1

# Create the directory for the Salmon index
mkdir -p salmon_index
path="salmon_index"

# Download the transcriptome and move it to the designated directory
wget -P $path ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip $path/Homo_sapiens.GRCh38.cdna.all.fa.gz

# Download the GTF file and move it to the designated directory
wget -P $path ftp://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz
gunzip $path/Homo_sapiens.GRCh38.104.gtf.gz

# Create the Salmon index in the designated directory
salmon index -t $path/Homo_sapiens.GRCh38.cdna.all.fa -i $path/salmon_index_human


