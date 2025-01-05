#!/bin/bash
#SBATCH --job-name=publicdata_deseq_objects
#SBATCH --chdir=/nemo/lab/patanir/home/users/palk/Projects/0006-phase-condensation-of-repeats-in-neurons//analysis/00-manuscript/13-03-public-data/
#SBATCH --output=/nemo/lab/patanir/home/users/palk/Projects/0006-phase-condensation-of-repeats-in-neurons//analysis/00-manuscript/13-03-public-data/slurm/panals-deseq-object.stdout
#SBATCH --error=/nemo/lab/patanir/home/users/palk/Projects/0006-phase-condensation-of-repeats-in-neurons//analysis/00-manuscript/13-03-public-data/slurm/panals-deseq-object.stderr
#SBATCH --get-user-env
#SBATCH --mail-type=ALL
#SBATCH --mail-user=palk@crick.ac.uk
#SBATCH --nodes=1
#SBATCH --partition=cpu
#SBATCH --time=72:00:00
#SBATCH --mem=128GB

source /nemo/lab/patanir/home/users/palk/anaconda3/etc/profile.d/conda.sh

conda activate r_4.3.1

Rscript scripts/02-02-save-deseq2-differential-objects-for-all-samples-als-vs-control.R