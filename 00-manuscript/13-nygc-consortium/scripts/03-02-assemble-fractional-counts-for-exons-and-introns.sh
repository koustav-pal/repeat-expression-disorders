#!/bin/bash
#SBATCH --job-name=nygc-assemble-counts
#SBATCH --chdir=/nemo/lab/patanir/home/users/palk/Projects/0006-phase-condensation-of-repeats-in-neurons/analysis/00-manuscript/13-nygc-consortium/
#SBATCH --output=/nemo/lab/patanir/home/users/palk/Projects/0006-phase-condensation-of-repeats-in-neurons/analysis/00-manuscript/13-nygc-consortium/slurm/assemble_counts.stdout
#SBATCH --error=/nemo/lab/patanir/home/users/palk/Projects/0006-phase-condensation-of-repeats-in-neurons/analysis/00-manuscript/13-nygc-consortium/slurm/assemble_counts.stderr
#SBATCH --get-user-env
#SBATCH --mail-type=ALL
#SBATCH --mail-user=palk@crick.ac.uk
#SBATCH --nodes=1
#SBATCH --partition=hmem
#SBATCH --time=72:00:00
#SBATCH --mem=1000GB

source /nemo/lab/patanir/home/users/palk/anaconda3/etc/profile.d/conda.sh

conda activate r_4.3.1

Rscript scripts/03-02-assemble-fractional-counts-for-exons-and-introns.R