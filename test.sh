#!/bin/bash

#SBATCH --job-name=test_snakemake
#SBATCH --partition=fast

#SBATCH --output=logtest.out

module load mummer4/4.0.0rc1

nucmer --delta /shared/home/sbache/annotation_fusarium/assembly/snakemake_test/output_dir/mummer/CAV3049_FLYE-STEP_CORRECTION_MEDAKA.delta /shared/home/sbache/annotation_fusarium/assembly/snakemake_test/reference/UK0001.fasta /shared/home/sbache/annotation_fusarium/assembly/snakemake_test/fasta_files/CAV3049_FLYE-STEP_CORRECTION_MEDAKA.fasta


