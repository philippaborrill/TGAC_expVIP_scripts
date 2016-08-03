#!/bin/bash
#SBATCH -p nbi-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 100000 # memory pool for all cores
#SBATCH -t 2-00:00 # time (D-HH:MM)
#SBATCH -o /nbi/group-data/ifs/NBI/Cristobal-Uauy/PB_collaborations/TGAC_genome/scripts_used/Fu_Hao_genes.%N.%j.out # STDOUT
#SBATCH -e /nbi/group-data/ifs/NBI/Cristobal-Uauy/PB_collaborations/TGAC_genome/scripts_used/Fu_Hao_genes.%N.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=philippa.borrill@jic.ac.uk # send-to address
#SBATCH -J Fu_Hao_genes

cd /nbi/Research-Groups/NBI/Cristobal-Uauy/expression_browser/TGAC_assembly/analysis/for_Fu_Hao

#grep -f 247.gene.ids.for.synteny.txt choulet_data.txt > 247_genes_in_choulet_data.txt

grep -f 247.gene.ids.for.synteny.txt replicated_data.txt > 247_genes_in_replicated_data.txt
