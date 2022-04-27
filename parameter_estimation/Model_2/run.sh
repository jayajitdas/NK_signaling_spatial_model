#!/bin/sh
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH -o slurm.%j.out # STDOUT
#SBATCH -e slurm.%j.err # STDERR

set -e

module load GCC/7.3.0-2.30 OpenMPI/3.1.1 Python/2.7.15
python run_pso.py
