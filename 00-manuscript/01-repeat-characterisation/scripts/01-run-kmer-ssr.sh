#!/bin/bash
#SBATCH --job-name=run_kmer_ssr
#SBATCH --chdir=/camp/lab/luscomben/home/users/palk/Projects/0006-phase-condensation-of-repeats-in-neurons/analysis/00-manuscript/01-repeat-characterisation/
#SBATCH --output=/camp/home/palk/home/users/palk/Projects/0006-phase-condensation-of-repeats-in-neurons/analysis/00-manuscript/01-repeat-characterisation/slurm/kmer-ssr/run_kmer_ssr.stdout
#SBATCH --error=/camp/home/palk/home/users/palk/Projects/0006-phase-condensation-of-repeats-in-neurons/analysis/00-manuscript/01-repeat-characterisation/slurm/kmer-ssr/run_kmer_ssr.stderr
#SBATCH --get-user-env
#SBATCH --mail-type=ALL
#SBATCH --mail-user=palk@crick.ac.uk
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=8GB

cd /camp/lab/luscomben/home/users/palk/Projects/0006-phase-condensation-of-repeats-in-neurons/analysis/00-manuscript/01-repeat-characterisation/
./packages/Kmer-SSR/bin/kmer-ssr -s "CAG,GCA,GCT,GCC,GCG,CTG,GGC,CCG,CGG,GAA,CAA" -p 3 -i input_files/00-all_genic_sequences_chr1-22.fa -o input_files/02-00-nonoverlapping_genic_sequences_chr1-22_trinucleotide.fa.tsv

./packages/Kmer-SSR/bin/kmer-ssr -s "CCTG" -p 4 -i input_files/00-all_genic_sequences_chr1-22.fa -o input_files/02-00-nonoverlapping_genic_sequences_chr1-22_quad.sense.fa.tsv

./packages/Kmer-SSR/bin/kmer-ssr -s "AAGGG,ATTCT,TGGAA,TTCCA,TTTCA" -p 5 -i input_files/00-all_genic_sequences_chr1-22.fa -o input_files/02-00-nonoverlapping_genic_sequences_chr1-22_penta.sense.fa.tsv

./packages/Kmer-SSR/bin/kmer-ssr -s "GGTCCT,GGTCCA,GGTCCC,GGTCCG,GGCCCT,GGCCCA,GGCCCC,GGCCCG,GGACCT,GGACCA,GGACCC,GGACCG,GGGCCT,GGGCCA,GGGCCC,GGGCCG" -p 6 -i input_files/00-all_genic_sequences_chr1-22.fa -o input_files/02-00-nonoverlapping_genic_sequences_chr1-22_GP.sense.fa.tsv
./packages/Kmer-SSR/bin/kmer-ssr -s "GGTGCT,GGTGCA,GGTGCC,GGTGCG,GGCGCT,GGCGCA,GGCGCC,GGCGCG,GGAGCT,GGAGCA,GGAGCC,GGAGCG,GGGGCT,GGGGCA,GGGGCC,GGGGCG" -p 6 -i input_files/00-all_genic_sequences_chr1-22.fa -o input_files/02-00-nonoverlapping_genic_sequences_chr1-22_GA.sense.fa.tsv
./packages/Kmer-SSR/bin/kmer-ssr -s "GGTCGT,GGTCGC,GGTCGA,GGTCGG,GGTAGA,GGTAGG,GGCCGT,GGCCGC,GGCCGA,GGCCGG,GGCAGA,GGCAGG,GGACGT,GGACGC,GGACGA,GGACGG,GGAAGA,GGAAGG,GGGCGT,GGGCGC,GGGCGA,GGGCGG,GGGAGA,GGGAGG" -p 6 -i input_files/00-all_genic_sequences_chr1-22.fa -o input_files/02-00-nonoverlapping_genic_sequences_chr1-22_GR.sense.fa.tsv
./packages/Kmer-SSR/bin/kmer-ssr -s "CCTGCT,CCTGCA,CCTGCC,CCTGCG,CCAGCT,CCAGCA,CCAGCC,CCAGCG,CCCGCT,CCCGCA,CCCGCC,CCCGCG,CCGGCT,CCGGCA,CCGGCC,CCGGCG" -p 6 -i input_files/00-all_genic_sequences_chr1-22.fa -o input_files/02-00-nonoverlapping_genic_sequences_chr1-22_PA.sense.fa.tsv
./packages/Kmer-SSR/bin/kmer-ssr -s "CCTCGT,CCTCGC,CCTCGA,CCTCGG,CCTAGA,CCTAGG,CCACGT,CCACGC,CCACGA,CCACGG,CCAAGA,CCAAGG,CCCCGT,CCCCGC,CCCCGA,CCCCGG,CCCAGA,CCCAGG,CCGCGT,CCGCGC,CCGCGA,CCGCGG,CCGAGA,CCGAGG" -p 6 -i input_files/00-all_genic_sequences_chr1-22.fa -o input_files/02-00-nonoverlapping_genic_sequences_chr1-22_PR.sense.fa.tsv
./packages/Kmer-SSR/bin/kmer-ssr -s "GGCCTG" -p 6 -i input_files/00-all_genic_sequences_chr1-22.fa -o input_files/02-00-nonoverlapping_genic_sequences_chr1-22_GL.sense.fa.tsv
./packages/Kmer-SSR/bin/kmer-ssr -s "CCCTCT" -p 6 -i input_files/00-all_genic_sequences_chr1-22.fa -o input_files/02-00-nonoverlapping_genic_sequences_chr1-22_PS.sense.fa.tsv

./packages/Kmer-SSR/bin/kmer-ssr -s "CCCCGCCCCGCG" -p 12 -i input_files/00-all_genic_sequences_chr1-22.fa -o input_files/02-00-nonoverlapping_genic_sequences_chr1-22_dodeca.sense.fa.tsv