#!/bin/bash
#SBATCH --job-name=answerals_deseq_objects
#SBATCH --chdir=/camp/lab/patanir/home/users/palk/Projects/0006-phase-condensation-of-repeats-in-neurons/analysis/00-manuscript/13-01-answerals/
#SBATCH --output=/camp/lab/patanir/home/users/palk/Projects/0006-phase-condensation-of-repeats-in-neurons/analysis/00-manuscript/13-01-answerals/slurm/save-deseq-object.stdout
#SBATCH --error=/camp/lab/patanir/home/users/palk/Projects/0006-phase-condensation-of-repeats-in-neurons/analysis/00-manuscript/13-01-answerals/slurm/save-deseq-object.stderr
#SBATCH --get-user-env
#SBATCH --mail-type=ALL
#SBATCH --mail-user=palk@crick.ac.uk
#SBATCH --nodes=1
#SBATCH --partition=cpu
#SBATCH --time=72:00:00
#SBATCH --mem=128GB

source /nemo/lab/patanir/home/users/palk/anaconda3/etc/profile.d/conda.sh

conda activate r_4.3.1

Rscript /camp/home/palk/home/Projects/0006-phase-condensation-of-repeats-in-neurons/analysis/00-manuscript/13-01-answerals/scripts/02-01-save-deseq2-differential-objects.R