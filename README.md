# Linear pangenome pipeline

This repository contains a copy of the Snakemake pipeline used for the creation of the four linear _Lactuca_ pangenomes. 

## Running the pipeline

The pipeline is written in Snakemake and contains dependencies from both conda and singularity.
It is thus necessary that both conda (mamba) and singularity are installed.

Next, for running the pipeline, go to the directory of a species and run it (I recommend running it in a screen since it can take weeks to finish):

```bash
cd Lactuca_virosa #as an example, it is the same for the other species
snakemake -rpc81 -s ../Snakefile --use-conda --use-singularity
```

Breakdown of the command above (and other options). Please always check the Snakemake documentation of your version of Snakemake. (This pipeline has been tested with version 6 and 7 of Snakemake.)

| Parameter              | Effect                                                                                       |
| ---------------------- | -------------------------------------------------------------------------------------------- |
| `-r`                   | Prints the reason why a rule is executed.                                                    |
| `-p`                   | Prints the command a rule executes.                                                          |
| `-c`                   | Specifies the total number of cores/threads to use.                                          |
| `-s`                   | Specifies which `Snakefile` to use.                                                          |
| `--use-conda`          | Lets Snakemake use conda/mamba.                                                              |
| `--use-singularity`    | Lets Snakemake use singularity.                                                              |
| `-n`                   | Dry-run; won't execute anything.                                                             |
| `--conda-prefix`       | Specifies where to store conda environments; useful for preventing redundant environments.   |
| `--singularity-prefix` | Specifies where to store singularity containers; useful for preventing redundant containers. |

## Example command read mapping

Since the first step of this pipeline is read mapping to the respective reference genomes, one can only perform this step like so:

```bash
cd Lactuca_sativa #as an example, it is the same for the other species
snakemake -rpc81 -s ../Snakefile --use-conda output/alignment_merged/LK122.bam.bai #using LK122 as an example, it is the same for other accessions
```

**NB**: Please note that this pipeline does not add read groups or mark duplicates.


## Copyright statement
    Copyright (C)  2024  Dirk-Jan M. van Workum

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

