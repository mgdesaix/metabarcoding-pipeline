# Welcome to the metabarcoding pipeline!

The goal of this pipeline is to provide users a quick, easy way to run through FASTQ files and get the basic summaries they need for metabarcoding data. The pipeline is set up in a couple of main steps for the user to run:

1. Process the raw data

2. Summarize the processed data

For each step, the pipeline has a configuration file (\*.config.txt) the user needs to edit and is run with the Bash script (\*.sh) in the Terminal, which the user does *not* need to edit. The details of these steps can be found in the following sections. In short, running the pipeline looks like:

```sh
# Step 1
./MBpipeline_01_raw2processed.sh MBpipeline_01_config.txt
# Step 2
./MBpipeline_02_processed2summary.sh MBpipeline_02_config.txt
```



### Installation

The pipeline uses the software `trimmomatic` and `vsearch` for processing the metabarcoding data, as well as `R` for the final tidying of the data summaries. For ease of execution, the user can use the YAML file (environment.yml) provided here to install these dependencies in a [Mamba](https://mamba.readthedocs.io/en/latest/index.html)/[Conda](https://docs.conda.io/projects/conda/en/stable/) environment which will be called when the pipeline is run:

```sh
git clone https://github.com/mgdesaix/metabarcoding-pipeline.git
cd metabarcoding-pipeline
conda env create -f environment.yml
```

### Step 1) Process the raw data

To run the first step of the pipeline some preparation of files and directories is necessary.

First, the raw \*.fastq files **must be in a directory called "01_raw"** AND that directory needs to be in the same directory as the metabarcoding-pipeline scripts.

The configuration file (`MBpipeline_01_config.txt`) also needs to be edited. The configuration file is set up such that there is a multi-line header delineated by rows beginning with "#", followed by rows providing different parameters. Each parameter is on its own row and the name to the left of the "=" is *not to be edited* while the value to the right of the "=" is what should be edited. The parameters in this step are the following:

**conda** = 

**R1_file** = 

**illumina_clip** = 

**CROP** = 

**SLIDINGWINDOW** = 

**MINLEN** = 

**a** = 

**g** = 

### Step 2) Summarize the processed data








