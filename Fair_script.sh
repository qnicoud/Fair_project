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

for j in $(tail -n +2 ../../PRJNA304086true.txt)
do
	# get important informations from the line
	access = $(echo "$j" | cut -f6)
	id = $(echo "$j" | cut -f1)
	md5 = $(echo "$j" | cut -f7)

	echo "-----------------------------------------------------------------------"
	echo ${id}
	echo "-----------------------------------------------------------------------"

	# Download file
	wget ${acces}

	# Get md5 of downloaded file
	md5_local = "$(md5sum $id.fastq.gz | cut -d' ' -f1)"
	echo $md_local

	# Test md5
	if [ "$md5_local" == "$md5" ]
	then
		echo "Done"
	else
		echo "Nope"
		exit
	fi
done

cd ../..

