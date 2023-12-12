#0!/bin/bash

##################################################
####### Metabarcoding pipeline: Part 2 ###########
## Matt DeSaix
## October 31, 2023
##################################################
# This is the second part of the metabarcoding pipeline and it summarizes the barcode data.
##################################################
# To run this pipeline you NEED:
# 1) Trimmomatic and vsearch software available in conda environment (as called below) or available locally
# 2) The raw fastq files to be processed in a directory path from this script as "./01_raw/"
# 3) Read1 files must have *R1_001.fastq suffix, and R2 files identical names (except "R1" = "R2")

source ~/.bash_profile

echo "Welcome to the metabarcoding pipeline!"
# Check that this file was supplied
if [ -z "$1" ]
then
  echo "No config file supplied"
  exit 1
fi
config_file=$1

# activate Conda for relevant software
conda=$(awk -v FS="=" '$1 == "conda" {print $2}' ${config_file})
conda activate ${conda}

# file with all the names of processed/filtered fasta files
processed_file=$(awk -v FS="=" '$1 == "processed_file" {print $2}' ${config_file})

word_count=$(wc -l ${processed_file})
ind=$(echo ${word_count} | cut -f1 -d' ')
echo "Processing ${ind} individual files"

#######################################################
######### Step 1: Cluster OTUs (vsearch) ##############

indir=./02_processed/05_uniques
outdir=./03_results/01_cluster
mkdir -p ${outdir}

cluster_id=$(awk -v FS="=" '$1 == "cluster_id" {print $2}' ${config_file})

echo "Running vsearch -cluster_unoise with -id ${cluster_id}"
echo "Output to: ${outdir}/"
for uniques in `cat ${processed_file}`
do
	otus=${uniques//uniques.fasta/otus.fasta}
	echo ${otus} >> ${outdir}/otus-clusters.txt 
	vsearch -cluster_unoise ${indir}/${uniques} -id ${cluster_id} -minsize 2 -unoise_alpha 2 -sizein \
		-sizeout -fasta_width 0 -centroids ${outdir}/${otus} -relabel Otu >> ${outdir}/otus-clusters.out 2>&1
	echo ${otus} "DONE"
done


########################################################################
####### Step 2: Blast cluster/otus to database (vsearch) ###############

indir=${outdir}
outdir=./03_results/02_blast
mkdir -p ${outdir}

db=$(awk -v FS="=" '$1 == "db" {print $2}' ${config_file})
id_similarity=$(awk -v FS="=" '$1 == "id_similarity" {print $2}' ${config_file})

echo "Running vsearch -usearch_global with -d ${id_similarity} -db ${db}"
echo "Output to: ${outdir}/"
for otus in `cat ${indir}/otus-clusters.txt`
do
	blast=${otus//otus.fasta/blast.txt}
	echo ${blast} >> ${outdir}/blast-samples.txt 
	vsearch -usearch_global ${indir}/${otus} -db ${db} -id ${id_similarity} -strand both -sizein -sizeout \
		-fasta_width 0 -top_hits_only -blast6out ${outdir}/${blast} >> ${outdir}/blast-clusters.out 2>&1
	echo ${blast} "DONE"
done


########################################################################
####### Step 3: Taxonomic classification (vsearch) ###############

indir=./03_results/01_cluster
outdir=./03_results/03_taxon_class
mkdir -p ${outdir}

sintax_db=$(awk -v FS="=" '$1 == "sintax_db" {print $2}' ${config_file})


echo "Running vsearch -sintax"
echo "Output to: ${outdir}/"
for otus in `cat ${indir}/otus-clusters.txt`
do
	sintax=${otus//otus.fasta/sintax.tsv}
	echo ${sintax} >> ${outdir}/sintax-samples.txt 
	vsearch -sintax ${indir}/${otus} -db ${sintax_db} -tabbedout ${outdir}/${sintax} -strand both -sintax_cutoff 0.6 >> ${outdir}/vsearch-sintax.out 2>&1
	echo ${sintax} "DONE"
done


########################################################################
####### Step 4: Summarize ###############

outdir=./03_results/04_summary
mkdir -p ${outdir}

for blast in `ls ./03_results/02_blast/*_blast.txt`
do
	id=$(echo ${blast} | cut -f4 -d/ | sed 's/_blast.txt//')
	awk -v ID=${id} '{print ID, $1, $2, $3}' ${blast} | sed 's/;/ /g' | sed 's/size=//g' >> ./03_results/04_summary/blast_summary.txt
	
	awk -v ID=${id} '{print ID, $1, $2}' ./03_results/03_taxon_class/${id}_sintax.tsv | sed 's/,/ /g' >> ./03_results/04_summary/sintax_summary.txt
done













