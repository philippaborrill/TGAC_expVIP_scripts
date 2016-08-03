#!/bin/bash
#
# SLURM script to launch fasta index
#
#SBATCH -p nbi-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 30000 # memory pool for all cores
#SBATCH -t 2-00:00 # time (D-HH:MM)
#SBATCH -o /nbi/Research-Groups/NBI/Cristobal-Uauy/expression_browser/TGAC_assembly/analysis/slurm_output/fasta_index.%N.%j.out # STDOUT
#SBATCH -e /nbi/Research-Groups/NBI/Cristobal-Uauy/expression_browser/TGAC_assembly/analysis/slurm_output/fasta_index.%N.%j.err # STDERR
#SBATCH -J fasta_index
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=philippa.borrill@jic.ac.uk # send-to address

cd  /nbi/Research-Groups/NBI/Cristobal-Uauy/CS_improved/

source samtools 1.3

samtools faidx Triticum_aestivum_CS42_TGACv1_all.fa

