##--------------------------------------------------------------------------------------------------------------------
## FAIR script
## Author: T. Denecker & C. Toffano-Nioche
## affiliation: I2BC
## Aim: A workflow to process RNAseq.
## Organism: O.tauri
## Date: Jan 2019
## Step:
##	1- Create tree structure
##	2- Sownload data from SRA
##--------------------------------------------------------------------------------------------------------------------

echo "================================================================================"
echo "Creation of tree structure"
echo "================================================================================"

mkdir Project
mkdir Project/samples
mkdir Project/annotations
mkdir Project/bowtie2
mkdir Project/fastqc
mkdir Project/genome
mkdir Project/graphics
mkdir Project/htseq
mkdir Project/reference
mkdir Project/samtools

echo "================================================================================"
echo "Download data from SRA"
echo "================================================================================"

cd Project/samples

IFS=$'\n'		# Make newlines the only separator

for j in $(tail -n +2 ../../conditions.txt)
do
	# get important informations from the line
	access=$(echo "${j}" | cut -f6)
	id=$(echo "${j}" | cut -f1)
	md5=$(echo "${j}" | cut -f7)

	echo "-----------------------------------------------------------------------"
	echo ${id}
	echo "-----------------------------------------------------------------------"

	# Download file
	wget ${access}

	# Get md5 of downloaded file
	md5_local="$(md5sum ${id}.fastq.gz | cut -d' ' -f1)"
	echo ${md5_local}
	echo ${md5}
	echo ${#md5_local}
	echo ${#md5}

	md5Test=${md5} | tr -d ' '

	# Test md5
	if [ "${md5_local}" == "${md5test}" ]
	then
		echo "Done"
	else
		echo "ERROR : ${id}.fasta.gz corrupted"
		exit 1
	fi
done

cd ../..

echo "======================================================================================"
echo " Download annotations"
echo "======================================================================================"

wget https://raw.githubusercontent.com/thomasdenecker/FAIR_Bioinfo/master/Data/O.tauri_annotation.gff -P Project/annotations

echo "======================================================================================"
echo " Download genome"
echo "======================================================================================"

wget https://raw.githubusercontent.com/thomasdenecker/FAIR_Bioinfo/master/Data/O.tauri_genome.fna -P Project/genome

