# Linear pangenome pipeline

This repository contains a copy of the Snakemake pipeline used for the creation of the four linear _Lactuca_ pangenomes. 

## Running the pipeline

The pipeline is written in Snakemake and contains dependencies from both conda and singularity.
It is thus necessary that both conda (mamba) and singularity are installed.

Next, for running the pipeline, go to the directory of a species and run it (I recommend running it in a screen since it can take weeks to finish):
```bash
cd Lactuca_virosa #as an example, it is the same for the other species
snakemake -nrpc81 -s ../Snakefile --use-conda --use-singularity
```
