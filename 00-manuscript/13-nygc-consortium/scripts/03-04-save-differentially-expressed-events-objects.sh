#!/bin/bash
#SBATCH --job-name=save_diff_expr_objects
#SBATCH --chdir=/camp/lab/patanir/home/users/palk/Projects/0006-phase-condensation-of-repeats-in-neurons/analysis/00-manuscript/13-nygc-consortium/
#SBATCH --output=/camp/lab/patanir/home/users/palk/Projects/0006-phase-condensation-of-repeats-in-neurons/analysis/00-manuscript/13-nygc-consortium/slurm/save_diff_expr_objects.stdout
#SBATCH --error=/camp/lab/patanir/home/users/palk/Projects/0006-phase-condensation-of-repeats-in-neurons/analysis/00-manuscript/13-nygc-consortium/slurm/save_diff_expr_objects.stderr
#SBATCH --get-user-env
#SBATCH --mail-type=ALL
#SBATCH --mail-user=palk@crick.ac.uk
#SBATCH --nodes=1
#SBATCH --partition=hmem
#SBATCH --time=72:00:00
#SBATCH --mem=1000GB
#SBATCH --dependency=afterok:66811996

source /nemo/lab/patanir/home/users/palk/anaconda3/etc/profile.d/conda.sh

conda activate r_4.1.2

Rscript /nemo/lab/patanir/home/users/palk/Projects/0006-phase-condensation-of-repeats-in-neurons/analysis/00-manuscript/13-nygc-consortium/scripts/03-04-save-differentially-expressed-events-objects.R