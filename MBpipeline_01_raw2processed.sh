#!/bin/bash

##################################################
####### Metabarcoding pipeline: Part 1 ###########
## Matt DeSaix
## October 31, 2023
##################################################
# This is the first part of the metabarcoding pipeline and it cleans up the raw read files.
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

# file with all the names of raw read 1 fastqs
R1_file=$(awk -v FS="=" '$1 == "R1_file" {print $2}' ${config_file})
word_count=$(wc -l ${R1_file})
ind=$(echo ${word_count} | cut -f1 -d' ')
echo "Processing ${ind} individual files"
# ILLUMINACLIP file for trimmomatic
illumina_clip=$(awk -v FS="=" '$1 == "illumina_clip" {print $2}' ${config_file})

##################################################
####### Step 1: Trimming (trimmomatic) ###################

indir=./01_raw
outdir=./02_processed/01_paired
mkdir -p ${outdir}

CROP=$(awk -v FS="=" '$1 == "CROP" {print $2}' ${config_file})
SLIDINGWINDOW=$(awk -v FS="=" '$1 == "SLIDINGWINDOW" {print $2}' ${config_file})
MINLEN=$(awk -v FS="=" '$1 == "MINLEN" {print $2}' ${config_file})

echo "Running Trimmomatic with Crop:${CROP} SLIDINGWINDOW:${SLIDINGWINDOW} MINLEN:${MINLEN}"
echo "Output to: ${outdir}/"
for R1 in `cat ${R1_file}`
do
	R2=${R1//R1_001.fastq/R2_001.fastq} # MGD question: Will R1 files never have "R1" text in name apart from suffix?
	R1paired=${R1//.fastq/_paired.fastq}
	echo ${R1paired} >> ${outdir}/paired-R1.txt # save off R1paired file names
	R1unpaired=${R1//.fastq/_unpaired.fastq}
	R2paired=${R2//.fastq/_paired.fastq}
	R2unpaired=${R2//.fastq/_unpaired.fastq}
	trimmomatic PE -threads 32 -phred33 ${indir}/$R1 ${indir}/$R2 ${outdir}/$R1paired \
		${outdir}/$R1unpaired ${outdir}/$R2paired ${outdir}/$R2unpaired \
		ILLUMINACLIP:${illumina_clip}:2:30:10 CROP:250 SLIDINGWINDOW:4:15 MINLEN:75 >> ${outdir}/trimmomatic.out 2>&1
	echo ${R1paired} "DONE"
done

paste ${outdir}/paired-R1.txt <(grep "Input Read Pairs" ${outdir}/trimmomatic.out) \
	| awk '{print "Step1", $0}' > ${outdir}/trimmomatic.summary.txt


##################################################
####### Step 2: Merging reads (vsearch) ###################

indir=${outdir}
outdir=./02_processed/02_merged
mkdir -p ${outdir}

echo "Running vsearch to merge"
echo "Output to: ${outdir}/"
for R1 in `cat ${indir}/paired-R1.txt`
do
	R2=${R1//R1_001_paired.fastq/R2_001_paired.fastq}
	merged=$(echo ${R1} | sed 's/R1_001_paired.fastq/merged.fastq/' | sed 's/L001_merged/merged/')
	echo ${merged} >> ${outdir}/merged-reads.txt
	vsearch -fastq_mergepairs ${indir}/$R1 -reverse ${indir}/$R2 -fastqout ${outdir}/$merged \
		-fastq_allowmergestagger >> ${outdir}/vsearch.merge.out 2>&1
	echo ${merged} "DONE"
done

grep -A 3 "Merging reads" ${outdir}/vsearch.merge.out | grep "Pairs" > ${outdir}/pairs.column.txt
grep -A 3 "Merging reads" ${outdir}/vsearch.merge.out | grep "Merged" > ${outdir}/merged.column.txt
grep -A 3 "Merging reads" ${outdir}/vsearch.merge.out | grep "Not merged" > ${outdir}/not-merged.column.txt
paste ${outdir}/merged-reads.txt ${outdir}/pairs.column.txt ${outdir}/merged.column.txt ${outdir}/not-merged.column.txt \
	| awk '{print "Step2", $0}' > ${outdir}/vsearch.merge.summary.txt
rm ${outdir}/*.column.txt


##################################################
####### Step 3: Cut adapters (cutadapt) ###################

indir=${outdir}
outdir=./02_processed/03_cutadapt
mkdir -p ${outdir}

a=$(awk -v FS="=" '$1 == "a" {print $2}' ${config_file})
g=$(awk -v FS="=" '$1 == "g" {print $2}' ${config_file})

echo "Running cutadapt to cut adapters: ${a} ${g}"
echo "Output to: ${outdir}/"
for merged in `cat ${indir}/merged-reads.txt`
do
	stripped=${merged//merged.fastq/stripped.fastq}
	echo ${stripped} >> ${outdir}/stripped-reads.txt
	cutadapt ${a} ${g} -o ${outdir}/${stripped} --discard-untrimmed ${indir}/${merged} -u 6 \
		>> ${outdir}/cutadapt.out 2>&1
	echo ${stripped} "DONE"
done

sed -n -e 's/^.*\(Reads with adapters: \).*  /\1/p' ${outdir}/cutadapt.out > ${outdir}/cutadapt.reads_adapters.column.txt
sed -n -e 's/^.*\(Reads written \).*  /\1/p' ${outdir}/cutadapt.out > ${outdir}/cutadapt.reads_written.column.txt
sed -n -e 's/^.*\(Total written \).*  /\1/p' ${outdir}/cutadapt.out > ${outdir}/cutadapt.total_written.column.txt

paste ${outdir}/stripped-reads.txt ${outdir}/cutadapt.reads_adapters.column.txt ${outdir}/cutadapt.reads_written.column.txt \
	${outdir}/cutadapt.total_written.column.txt | awk '{print "Step3", $0}' > ${outdir}/cutadapt.summary.txt
rm ${outdir}/*.column.txt

##################################################
####### Step 4: Filter errors (vsearch) ###################

indir=${outdir}
outdir=./02_processed/04_filtered
mkdir -p ${outdir}

echo "Running vsearch to filter errors"
echo "Output to: ${outdir}/"
for stripped in `cat ${indir}/stripped-reads.txt`
do
	filtered=${stripped//stripped.fastq/filtered.fasta}
	echo ${filtered} >> ${outdir}/filtered-reads.txt
	vsearch -fastq_filter ${indir}/${stripped} -fastq_maxee 1.0 -relabel Filt -fastaout ${outdir}/${filtered} \
		>> ${outdir}/filter.out 2>&1
	echo ${filtered} "DONE"
done

paste ${outdir}/filtered-reads.txt <(grep "sequences kept" ${outdir}/filter.out) \
	| awk '{print "Step4", $0}' > ${outdir}/filter.summary.txt

##################################################
####### Step 5: Dereplicate (vsearch) ###################

indir=${outdir}
outdir=./02_processed/05_uniques
mkdir -p ${outdir}

echo "Running vsearch to dereplicate"
echo "Output to: ${outdir}/"
for filtered in `cat ${indir}/filtered-reads.txt`
do
	uniques=${filtered//filtered.fasta/uniques.fasta}
	echo ${uniques} >> ${outdir}/uniques-reads.txt
	vsearch -derep_fulllength ${indir}/${filtered} -sizeout -relabel Uniq \
		-output ${outdir}/${uniques} >> ${outdir}/uniques.out 2>&1
	echo ${uniques} "DONE"
done

paste ${outdir}/uniques-reads.txt <(grep "unique sequences" ${outdir}/uniques.out) \
	| awk '{print "Step5", $0}' > ${outdir}/uniques.summary.txt

##################################################
####### Summarize ###################

outdir=./02_processed/06_summary
mkdir -p ${outdir}

echo "Summarizing results"
echo "Output to: ${outdir}/"

cat ./02_processed/*/*.summary.txt | sort -k1,1 > ${outdir}/MBpipeline_01_raw2processed.summary.txt








