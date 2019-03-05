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
dirlist=$(find Project/samples/*.fastq.gz)
echo  $dirlist
# Name of the file containing the reference genome
genome="./Project/genome/O.tauri_genome.fna"
# Name of the fil containing the annotation of the genome
annotation="./Project/annotations/O.tauri_annotation.gff"

for file in $dirlist
do
	# Name without path
	file_name="$(basename $file)"
	echo $file_name
	# Name without path and .gz
	nameFastq="${file_name%.*}"
	echo $nameFastq
	# Name without path, gz and fastq
	sample="${nameFastq%.*}"
	echo $sample

	echo "============================================================================================================"
	echo "Quality chech with fastqc of sample ${sample}"
	echo "============================================================================================================"
	fastqc Project/samples/${file_name} --outdir Project/fastqc/

	echo "============================================================================================================"
	echo "Indexing reference genome"
	echo "============================================================================================================"
	bowtie2-build ${genome} O_tauri

	echo "============================================================================================================"
	echo "Reads alignment on reference genome - sample ${sample}"
	echo "============================================================================================================"
	bowtie2 -x O_tauri -U Project/samples/${file_name} -S Project/bowtie2/bowtie-${sample}.sam 2> Project/bowtie2/bowtie-${sample}.out

	echo "============================================================================================================"
	echo "Converting files in binary, sorting and indexing aligned reads - sample ${sample}"
	echo "============================================================================================================"
	samtools view -b Project/bowtie2/bowtie-${sample}.sam > Project/samtools/bowtie-${sample}.bam
	samtools sort Project/samtools/bowtie-${sample}.bam -o Project/samtools/bowtie-${sample}.sorted.bam
	samtools index Project/samtools/bowtie-${sample}.sorted.bam

	echo "============================================================================================================"
	echo "Comptage -échantillon ${sample}"
	echo "============================================================================================================"
	htseq -count --stranded=no --type='gene' ==idattr='ID' --order=name --format=bam Project/samtools/bowtie-${sample}.sorted.bam ${annotations} > Project/htseq/count-${sample}.txt

	echo "============================================================================================================"
	echo "Cleaning useless files - sample ${sample}"
	echo "============================================================================================================"
	rm -f Project/samtools/bowtie-${sample}.sam Project/bowtie2/bowtie-${sample}.bam

done
