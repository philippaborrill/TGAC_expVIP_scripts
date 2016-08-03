#!/bin/bash
#
# SLURM script to launch compare_promoter_similiarity_triads_second_step.pl 
#
#SBATCH -p nbi-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 30000 # memory pool for all cores
#SBATCH -t 2-00:00 # time (D-HH:MM)
#SBATCH -o /nbi/Research-Groups/NBI/Cristobal-Uauy/expression_browser/TGAC_assembly/analysis/seedling_stress_analysis/promoter_similarity/slurm_output/calc_av_sim.%N.%j.out # STDOUT
#SBATCH -e /nbi/Research-Groups/NBI/Cristobal-Uauy/expression_browser/TGAC_assembly/analysis/seedling_stress_analysis/promoter_similarity/slurm_output/calc_av_sim.%N.%j.err # STDERR
#SBATCH -J calc_av_sim
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=philippa.borrill@jic.ac.uk # send-to address

cd /nbi/Research-Groups/NBI/Cristobal-Uauy/expression_browser/TGAC_assembly/scripts_used

perl compare_promoter_similarity_triads_second_step.pl 

