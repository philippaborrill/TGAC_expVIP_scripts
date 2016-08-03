#!/bin/bash
#
# SLURM script to launch gzip
#
#SBATCH -p nbi-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 30000 # memory pool for all cores
#SBATCH -t 2-00:00 # time (D-HH:MM)
#SBATCH -o /nbi/Research-Groups/NBI/Cristobal-Uauy/expression_browser/TGAC_assembly/analysis/seedling_stress_analysis/promoter_similarity/slurm_output/gzip.%N.%j.out # STDOUT
#SBATCH -e /nbi/Research-Groups/NBI/Cristobal-Uauy/expression_browser/TGAC_assembly/analysis/seedling_stress_analysis/promoter_similarity/slurm_output/gzip.%N.%j.err # STDERR
#SBATCH -J gzip
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=philippa.borrill@jic.ac.uk # send-to address

cd /nbi/Research-Groups/NBI/Cristobal-Uauy/expression_browser/TGAC_assembly/analysis/kallisto_results_bootstrap/results

gzip -c SRR1228245/abundance.h5 > SRR1228245.h5.gz
gzip -c SRR1228246/abundance.h5 > SRR1228246.h5.gz
gzip -c SRR1228247/abundance.h5 > SRR1228247.h5.gz
gzip -c SRR1228248/abundance.h5 > SRR1228248.h5.gz
gzip -c SRR1228249/abundance.h5 > SRR1228249.h5.gz
gzip -c SRR1228250/abundance.h5 > SRR1228250.h5.gz
gzip -c SRR1228251/abundance.h5 > SRR1228251.h5.gz
gzip -c SRR1228252/abundance.h5 > SRR1228252.h5.gz
gzip -c SRR1228253/abundance.h5 > SRR1228253.h5.gz
gzip -c SRR1228254/abundance.h5 > SRR1228254.h5.gz
gzip -c SRR1228255/abundance.h5 > SRR1228255.h5.gz
gzip -c SRR1228256/abundance.h5 > SRR1228256.h5.gz
gzip -c SRR1228257/abundance.h5 > SRR1228257.h5.gz
gzip -c SRR1228258/abundance.h5 > SRR1228258.h5.gz
gzip -c SRR1228259/abundance.h5 > SRR1228259.h5.gz
gzip -c SRR1228260/abundance.h5 > SRR1228260.h5.gz
gzip -c SRR1228261/abundance.h5 > SRR1228261.h5.gz
gzip -c SRR1228262/abundance.h5 > SRR1228262.h5.gz
gzip -c SRR1228263/abundance.h5 > SRR1228263.h5.gz
gzip -c SRR1228264/abundance.h5 > SRR1228264.h5.gz
gzip -c SRR1228265/abundance.h5 > SRR1228265.h5.gz
gzip -c SRR1542404/abundance.h5 > SRR1542404.h5.gz
gzip -c SRR1542405/abundance.h5 > SRR1542405.h5.gz
gzip -c SRR1542406/abundance.h5 > SRR1542406.h5.gz
gzip -c SRR1542407/abundance.h5 > SRR1542407.h5.gz
gzip -c SRR1542408/abundance.h5 > SRR1542408.h5.gz
gzip -c SRR1542409/abundance.h5 > SRR1542409.h5.gz
gzip -c SRR1542410/abundance.h5 > SRR1542410.h5.gz
gzip -c SRR1542411/abundance.h5 > SRR1542411.h5.gz
gzip -c SRR1542412/abundance.h5 > SRR1542412.h5.gz
gzip -c SRR1542413/abundance.h5 > SRR1542413.h5.gz
gzip -c SRR1542414/abundance.h5 > SRR1542414.h5.gz
gzip -c SRR1542415/abundance.h5 > SRR1542415.h5.gz
gzip -c SRR1542416/abundance.h5 > SRR1542416.h5.gz
gzip -c SRR1542417/abundance.h5 > SRR1542417.h5.gz
gzip -c ERR392055/abundance.h5 > ERR392055.h5.gz
gzip -c ERR392056/abundance.h5 > ERR392056.h5.gz
gzip -c ERR392057/abundance.h5 > ERR392057.h5.gz
gzip -c ERR392058/abundance.h5 > ERR392058.h5.gz
gzip -c ERR392059/abundance.h5 > ERR392059.h5.gz
gzip -c ERR392060/abundance.h5 > ERR392060.h5.gz
gzip -c ERR392061/abundance.h5 > ERR392061.h5.gz
gzip -c ERR392062/abundance.h5 > ERR392062.h5.gz
gzip -c ERR392063/abundance.h5 > ERR392063.h5.gz
gzip -c ERR392064/abundance.h5 > ERR392064.h5.gz
gzip -c ERR392065/abundance.h5 > ERR392065.h5.gz
gzip -c ERR392066/abundance.h5 > ERR392066.h5.gz
gzip -c ERR392067/abundance.h5 > ERR392067.h5.gz
gzip -c ERR392068/abundance.h5 > ERR392068.h5.gz
gzip -c ERR392069/abundance.h5 > ERR392069.h5.gz
gzip -c ERR392070/abundance.h5 > ERR392070.h5.gz
gzip -c ERR392071/abundance.h5 > ERR392071.h5.gz
gzip -c ERR392072/abundance.h5 > ERR392072.h5.gz
gzip -c ERR392073/abundance.h5 > ERR392073.h5.gz
gzip -c ERR392074/abundance.h5 > ERR392074.h5.gz
gzip -c ERR392075/abundance.h5 > ERR392075.h5.gz
gzip -c ERR392076/abundance.h5 > ERR392076.h5.gz
gzip -c ERR392077/abundance.h5 > ERR392077.h5.gz
gzip -c ERR392078/abundance.h5 > ERR392078.h5.gz
gzip -c ERR392079/abundance.h5 > ERR392079.h5.gz
gzip -c ERR392080/abundance.h5 > ERR392080.h5.gz
gzip -c ERR392081/abundance.h5 > ERR392081.h5.gz
gzip -c ERR392082/abundance.h5 > ERR392082.h5.gz
gzip -c ERR392083/abundance.h5 > ERR392083.h5.gz
gzip -c ERR392084/abundance.h5 > ERR392084.h5.gz



