#!/bin/bash

#SBATCH --job-name=mummer
#SBATCH --partition=fast

#SBATCH --output=logmummer.out

module load mummer4/4.0.0rc1 

nucmer -p UK001_GCA UK0001.fasta GCA_clean.fasta 

 
