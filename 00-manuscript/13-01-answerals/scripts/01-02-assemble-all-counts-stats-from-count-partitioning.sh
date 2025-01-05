#!/bin/bash
#SBATCH --job-name=assemble-counts
#SBATCH --chdir=/nemo/lab/patanir/home/users/palk/Projects/0006-phase-condensation-of-repeats-in-neurons/analysis/00-manuscript/13-01-answerals/
#SBATCH --output=/nemo/lab/patanir/home/users/palk/Projects/0006-phase-condensation-of-repeats-in-neurons/analysis/00-manuscript/13-01-answerals/slurm/assemble_counts.stdout
#SBATCH --error=/nemo/lab/patanir/home/users/palk/Projects/0006-phase-condensation-of-repeats-in-neurons/analysis/00-manuscript/13-01-answerals/slurm/assemble_counts.stderr
#SBATCH --get-user-env
#SBATCH --mail-type=ALL
#SBATCH --mail-user=palk@crick.ac.uk
#SBATCH --nodes=1
#SBATCH --partition=cpu
#SBATCH --time=72:00:00
#SBATCH --mem=128GB

source /nemo/lab/patanir/home/users/palk/anaconda3/etc/profile.d/conda.sh

conda activate r_4.3.1

Rscript scripts/01-02-assemble-all-counts-stats-from-count-partitioning.R