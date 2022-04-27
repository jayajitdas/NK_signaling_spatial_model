#!/bin/sh
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH -o slurm.%j.out # STDOUT
#SBATCH -e slurm.%j.err # STDERR

module load GCC/7.3.0-2.30 OpenMPI/3.1.1
mpirun -np 1 ./spk_redsky < in.model3_act
