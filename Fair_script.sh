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

if  [ -f "md5List.txt" ]
then
	touch "md5List.txt"
else
	echo "md5List.txt already exist"
fi
echo "=====================================================================================" >> "md5List.txt"
echo "Session of the " "$(date)" >> "md5List.txt"
echo "=====================================================================================" >> "md5List.txt"

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

	if [ -f Project/samples/"${id}.fastq.gz" ]
	then
		# Download file
		wget ${access}
	else
		echo "${id}.fastq.gz already exist"
	fi

	# Get md5 of downloaded file
	md5_local=$(md5sum ${id}.fastq.gz | cut -d' ' -f1)
	echo ${md5_local}
	echo ${md5}
	echo ${#md5_local}
	echo ${#md5}

	cd ../..
	echo "$id" >> "md5List.txt"
	echo "local" "$md5_local" >> "md5List.txt"
	echo "distant" "$md5" >> "md5List.txt"

	#md5Test=$( echo "${md5}" | sed 's/ //g' )
	#md5Test=$( echo "${md5}" | tr -d ' ' )
	md5Test=$(echo ${md5//[[:blank:]]/})

	echo ${#md5Test}

	# Test md5
	if [ "${md5_local}" == "${md5test}" ]
	then
		echo "Done"
		echo "md5 are similar" >> "md5List.txt"
	else
		echo "ERROR : ${id}.fasta.gz corrupted"
		echo "md5 are different" >> "md5List.txt"
	fi
	echo "---------------------------------------------------------" >> "md5List.txt"
	cd Project/samples
done

cd ../..

echo "======================================================================================"
echo " Download annotations"
echo "======================================================================================"
if [ -f Project/annotations/O.tauri_annotation.gff ]
then
	wget https://raw.githubusercontent.com/thomasdenecker/FAIR_Bioinfo/master/Data/O.tauri_annotation.gff -P Project/annotations
else
	echo "The file already exist."
fi

echo "======================================================================================"
echo " Download genome"
echo "======================================================================================"
if [ -f Project/annotations/O.tauri_genome.fna ]
then
	wget https://raw.githubusercontent.com/thomasdenecker/FAIR_Bioinfo/master/Data/O.tauri_genome.fna -P Project/genome
else
	echo "The file already exist"
fi
