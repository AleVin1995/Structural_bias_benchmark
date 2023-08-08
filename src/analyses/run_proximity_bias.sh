#!/usr/bin/env bash

#SBATCH --job-name=proximity_bias 
#SBATCH --mem-per-cpu=150GB
#SBATCH --partition=cpuq
#SBATCH --mail-type=NONE
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --chdir=/group/iorio/Alessandro/Structural_bias_benchmark
#SBATCH --output=/group/iorio/Alessandro/Structural_bias_benchmark/%j.out
#SBATCH --error=/group/iorio/Alessandro/Structural_bias_benchmark/%j.err

source ~/.bashrc
conda activate CN_bench_r

LIB=$1

if [[ -z "$LIB" ]]
then
	echo "No library argument provided"
	echo "Processing all libraries"

    Rscript src/analyses/proximity_bias_bw_test.r
else
    Rscript src/analyses/proximity_bias_bw_test.r $LIB
fi