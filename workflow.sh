#!/bin/bash

#SBATCH --job-name=test_snakemake
#SBATCH --partition=fast

#SBATCH --output=logtest.out

module purge
module load python/3.9

snakemake --profile test_snakemake --use-envmodules

