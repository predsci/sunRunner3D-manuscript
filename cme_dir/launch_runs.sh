#!/bin/bash
#SBATCH --job-name=PLUTO
#SBATCH -t 99:00:00
#SBATCH --export=ALL
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32


mpirun -np 32 ./pluto -no-x3par -restart 1 1>log 2>err



