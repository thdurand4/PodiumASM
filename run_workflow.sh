#!/bin/sh
################################ Slurm options #################################


### Job name
#SBATCH --job-name=assembly

### Requirements
#SBATCH --partition=long



### Output
#SBATCH --output=/shared/home/tdurand/annotation_fijiensism2/assembly_worflow/log_assembly.out
#SBATCH --error=/shared/home/tdurand/annotation_fijiensism2/assembly_worflow/log_assembly.err
module purge
module load python/3.7


snakemake --profile assembly --use-envmodules
