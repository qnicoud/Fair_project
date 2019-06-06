##--------------------------------------------------------------------------------------------------------------------
## tools_installation_script.sh
## Author: Q. Nicoud
## Affiliation: I2BC
## Aim: Install all the tools required for the workflow.
## Date: Feb 2019
## Step:
##	1- Install tree
##	2- Install conda and configure it
##	3- Install all the other required packages
##--------------------------------------------------------------------------------------------------------------------

## Install tree
sudo apt-get update
sudo apt-get install tree

## Install and configure conda
wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
rm Miniconda3-latest-Linux-x86_64.sh
	## Add in path
export PATH=~/miniconda3/bin:$PATH

	## Update
conda update conda
conda update --all

	## Add channels
conda config --add channels conda-forge
conda config --add channels bioconda

	## Install fastqc, bowtie2, htseq et samtools
conda install -y fastqc bowtie2 htseq samtools
