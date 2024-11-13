# Welcome to the metabarcoding pipeline!

The goal of this pipeline is to provide users a quick, easy way to run through FASTQ files and get the basic summaries they need for metabarcoding data. The pipeline is set up in a couple of main steps for the user to run:

1. Process the raw data

2. Summarize the processed data

3. Tidy the summary output

For each step, the pipeline has a configuration file (\*.config.txt) the user needs to edit and is run with the Bash script (\*.sh) in the Terminal, which the user does *not* need to edit. The details of these steps can be found in the following sections. In short, running the pipeline looks like:

```sh
# Step 0: Edit the configuration files (*config.txt) appropriately for your data

# Step 1: Process the raw data
./MBpipeline_01_raw2processed.sh MBpipeline_01_config.txt
# Step 2: Summarize the processed data
./MBpipeline_02_processed2summary.sh MBpipeline_02_config.txt
# Step 3: Tidy the summary output
# activate conda environment that has R and tidyverse
conda activate YOUR_ENVIRONMENT_NAME
# copy R summary script to folder for analysis
cp blast_summary_genbank.R ./03_results/04_summary/
# run the R summary script
cd ./03_results/04_summary/
Rscript blast_summary_genbank.R blast_summary.txt
```

The output from the `blast_summary_genbank.R` are placed in a folder `blast_summary_genbank` and are the following:

- Contingency table (`blast_species_contingency_table.txt`): Summary of number of reads for species (rows) vs. samples (columns)

- Sample depth summary (`blast_sample_depth.txt`): Summary of total reads per sample

- OTU sequence information (`blast_OTU_summary.txt`): Summary of each OTU delineated per sample and it's corresponding identity match to target, length, coverage of target, and actual sequence.

- OTU histogram plot (`OTU_histogram.png`)

*Note:* If you get a warning about files not being executable, remember to use `chmod +x "filename"` to make files executable.

### Installation

The pipeline uses the software `trimmomatic`, `cutadapt`, and `vsearch` for processing the metabarcoding data, as well as `R` for the final tidying of the data summaries. For ease of execution, the user can use the YAML file (environment.yml) provided here to install these dependencies in a [Mamba](https://mamba.readthedocs.io/en/latest/index.html)/[Conda](https://docs.conda.io/projects/conda/en/stable/) environment which will be called when the pipeline is run:

```sh
git clone https://github.com/mgdesaix/metabarcoding-pipeline.git
cd metabarcoding-pipeline
conda env create -f environment.yml
```

*Note:* If you have an Apple M1/M2 processor, you will need to install the conda environment with the following commands because many of the softwares we need don't have the ARM architecture compatibility.

```sh
CONDA_SUBDIR=osx-64 conda env create -f environment.yml
conda activate MBpipeline
conda config --env --set subdir osx-64
```

### Step 1) Process the raw data

To run the first step of the pipeline some preparation of files and directories is necessary.

First, the raw \*.fastq files **must be in a directory called "01_raw"** AND that directory needs to be in the same directory as the metabarcoding-pipeline scripts.

The configuration file (`MBpipeline_01_config.txt`) also needs to be edited. The configuration file is set up such that there is a multi-line header delineated by rows beginning with "#", followed by rows providing different parameters. Each parameter is on its own row and the name to the left of the "=" is *not to be edited* while the value to the right of the "=" is what should be edited. *NO SPACES AFTER THE '='!!* The definitions of the parameters for this step are provided below:

**conda** = The conda environment that has all necessary applications. If this was set up using the \*.yml file then it is "MBpipeline"

**R1_file** = A text file that provides the filename of all the raw read 1 files to be processed (which are located in the 01_raw directory, along with their read 2 complement). 

**illumina_clip** = The corresponding file for the `ILLUMINACLIP:` parameter from `trimmomatic`.

**CROP** = The `CROP:` parameter from `trimmomatic`.

**SLIDINGWINDOW** = The `SLIDINGWINDOW:` parameter from `trimmomatic`.

**MINLEN** = The `MINLEN:` parameter from `trimmomatic`.

**R1_primer** = The primer sequence of read 1, which is used with the `-g` option of `cutadapt` (at front of sequence) and is anchored with '^'.

**R2_primer** = The primer sequence of read 2, which is used with the `-G` option of `cutadapt` (at front of sequence) and is anchored with '^'.

**trim_bases** = The `-u` (read 1) or `-U` (read 2) parameter from `cutadapt` to cut a set number of bases at the front of the sequence prior to removing the primers. In our data we had 6 "filler" bases that always preceded the primer sequences.


### Step 2) Summarize the processed data

**conda** = [string] The conda environment that has all necessary applications. If this was set up using the \*.yml file then it is "MBpipeline"

**processed_file** = [string] The final processed reads produced from Step 1. By default, these files are saved in Step 1 as "./02_processed/05_uniques/uniques-reads.txt".

**cluster_id** = [real] `-id` parameter value from `vsearch` for clustering OTUs, default is set at 0.97.

**query_cov** = [real] `-query_cov` parameter value from `vsearch` for filtering matches by query coverage, default is set at 0.6 (i.e. matches with query coverage beneath 60% are dropped)

**db** = [string] Database file specified by `-db` parameter from `vsearch`.

**id_similarity** = [real] `-id` parameter value from `vsearch` for blast

**sintax_db** = [string] Database file from `-db` parameter from `vsearch`'s sintax specification.

**multiple_hits** = [logical] Providing the top hit from `vsearch`, or multiple hits. Default=NO, which does `-top_hits_only`. If changed to YES, then `-top_hits_only` is replaced by `-maxhits 0 -maxaccepts 0` which provides all database hits that pass the filters. This creates a `blast_summary_verbose.txt` file, which can be summarized with the corresponding R file.

### Step 3) Tidy the summary output

Run the `blast_summary_genbank.R` file from within the `./03_results/04_summary/` directory that is created by the `MBpipeline_02_processed2summary.sh` file. This R script will only work with data from queried to a Genbank style database, and will need to be tweaked to accommodate other styles of databases.

```sh
# copy R summary script to folder for analysis
cp blast_summary_genbank.R ./03_results/04_summary/
# run the R summary script
cd ./03_results/04_summary/
Rscript blast_summary_genbank.R blast_summary.txt
```

The above code creates a directory `./blast_summary/genbank` with the following files:

**blast_OTU_summary.txt** = A summary of the read depths for each OTU per sample and the associated taxonomic hit on the the OTU. Note that the OTU ids created by `vsearch` in the current workflow are not standardized across samples. I have added a column (`OTU_unique`) to this workflow that provides an OTU that is unique to the sequence and standardized across samples. The `OTU_fasta` column provides the original OTU designation that is specific to each sample and is relevant if you need to go back through the fasta files produced by `vsearch`.

**blast_sample_depth.txt** = A per-sample summary of the total number of reads (i.e., depth) analzed

**blast_species_contingency_table.txt** = A contingency table where rows are the taxonomic units identified and columns are the individual samples. Values in the table provide the total number of reads per sample associated with the given taxonomic unit.

**OTU_histogram.png** = A figure of the histogram of OTU sequence length we found to be useful for seeing the expected range of sequence length and identifying spurious outliers.

If you produced blast summary with multiple hits (`blast_summary_verbose.txt` when setting `multiple_hits=YES`), the corresponding R summarizing file, `blast_summary_genbank_verbose.R` reduces the hits to one species per sample, OTU, and match ID.

```sh
# copy R summary script to folder for analysis
cp blast_summary_genbank_verbose.R ./03_results/04_summary/
# run the R summary script
cd ./03_results/04_summary/
Rscript blast_summary_genbank_verbose.R blast_summary_verbose.txt
```




