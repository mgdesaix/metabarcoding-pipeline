# Welcome to the metabarcoding pipeline!

The goal of this pipeline is to provide users a quick, easy way to run through FASTQ files and get the basic summaries they need for metabarcoding data. The pipeline is set up in a couple of main steps for the user to run:

1. Process the raw data

2. Summarize the processed data

The details of these steps can be found in the following sections. In short, all the user needs to do to run the pipeline is to edit the configuration files ("*.config.txt") and then run:

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









