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
Rscript blast_summary_genbank.R blast_summary.txt
```

The output from the `blast_summary_genbank.R` file are:

- Contingency table (`blast_species_contingency_table.txt`): Summary of number of reads for species (rows) vs. samples (columns)

- Sample depth summary (`blast_sample_depth.txt`): Summary of total reads per sample

- OTU sequence information (`blast_OTU_summary.txt`): Summary of each OTU delineated per sample and it's corresponding identity match to target, length, coverage of target, and actual sequence.

*Note:* If you get a warning about files not being executable, remember to use `chmod +x "filename"` to make files executable.

### Installation

The pipeline uses the software `trimmomatic`, `cutadapt`, and `vsearch` for processing the metabarcoding data, as well as `R` for the final tidying of the data summaries. For ease of execution, the user can use the YAML file (environment.yml) provided here to install these dependencies in a [Mamba](https://mamba.readthedocs.io/en/latest/index.html)/[Conda](https://docs.conda.io/projects/conda/en/stable/) environment which will be called when the pipeline is run:

```sh
git clone https://github.com/mgdesaix/metabarcoding-pipeline.git
cd metabarcoding-pipeline
conda env create -f environment.yml
```

*Note:* If you have an Apple M1/M2 processor, you may get a warning about `cutadapt` not being available through Conda - if so, follow the [Cutadapt installation instructions](https://cutadapt.readthedocs.io/en/stable/installation.html).

### Step 1) Process the raw data

To run the first step of the pipeline some preparation of files and directories is necessary.

First, the raw \*.fastq files **must be in a directory called "01_raw"** AND that directory needs to be in the same directory as the metabarcoding-pipeline scripts.

The configuration file (`MBpipeline_01_config.txt`) also needs to be edited. The configuration file is set up such that there is a multi-line header delineated by rows beginning with "#", followed by rows providing different parameters. Each parameter is on its own row and the name to the left of the "=" is *not to be edited* while the value to the right of the "=" is what should be edited. The definitions of the parameters for this step are provided below:

**conda** = The conda environment that has all necessary applications. If this was set up using the \*.yml file then it is "MBpipeline"

**R1_file** = A text file that provides the filename of all the raw read 1 files to be processed (which are located in the 01_raw directory, along with their read 2 complement). 

**illumina_clip** = The corresponding file for the `ILLUMINACLIP:` parameter from `trimmomatic`.

**CROP** = The `CROP:` parameter from `trimmomatic`.

**SLIDINGWINDOW** = The `SLIDINGWINDOW:` parameter from `trimmomatic`.

**MINLEN** = The `MINLEN:` parameter from `trimmomatic`.

**a** = The `a` parameter from `cutadapt`.

**g** = The `g` parameter from `cutadapt`.


### Step 2) Summarize the processed data

**conda** = The conda environment that has all necessary applications. If this was set up using the \*.yml file then it is "MBpipeline"

**processed_file** = The final processed reads produced from Step 1. By default, these files are saved in Step 1 as "./02_processed/05_uniques/uniques-reads.txt".

**cluster_id** = `-id` parameter value from `vsearch` for clustering OTUs, default is set at 0.97.

**db** = Database file specified by `-db` parameter from `vsearch`.

**id_similarity** = `-id` parameter value from `vsearch` for blast

**sintax_db** = Database file from `-db` parameter from `vsearch`'s sintax specification.

### Step 3) Tidy the summary output

Run the `blast_summary_genbank.R` file from within the `./03_results/04_summary/` directory that is created by the `MBpipeline_02_processed2summary.sh` file. This R script will only work with data from queried to a Genbank style database, and will need to be tweaked to accomodate other styles of databases.





