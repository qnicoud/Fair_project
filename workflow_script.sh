#------------------------------------------------------------------------------------------------------------------
# Workflow script
# Author: T. Denecker, C. Toffano-Nioche
# Affiliation: I2BC
# Aim: A workflow to process RNA-seq data.
# Organism: O.tauri
# Date: Jan 2019
# Step:
#	1- Fastqc quality check
#	2- bowtie2 indexing of the samples
#	3- bowtie2 reads alignment
#-------------------------------------------------------------------------------------------------------------------

# Name of the files to analyze
dirlist=$(find Project/samples/*.fastq.gz
# Name of the file containing the reference genome
genome="./Project/genome/O.tauri_genome.fna"
# Name of the fil containing the annotation of the genome
annotation="./Project/annotations/O.tauri_annotation.gff"

for file in $dirlist
do
	# Name without path
	file_name="$(basename $file)"
	# 


echo "============================================================================================================"
echo "Quality chech with fastqc"
echo "============================================================================================================"

