#!/bin/bash
#
# SLURM batch script to launch order_genes
#
#SBATCH -p nbi-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 10000 # memory pool for all cores
#SBATCH -t 2-00:00 # time (D-HH:MM)
#SBATCH -o /nbi/Research-Groups/NBI/Cristobal-Uauy/expression_browser/TGAC_assembly/analysis/slurm_output/order_genes.%N.%j.out # STDOUT
#SBATCH -e /nbi/Research-Groups/NBI/Cristobal-Uauy/expression_browser/TGAC_assembly/analysis/slurm_output/order_genes.%N.%j.err # STDERR
#SBATCH -J order_genes
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=philippa.borrill@jic.ac.uk # send-to address


cd /nbi/Research-Groups/NBI/Cristobal-Uauy/expression_browser/TGAC_assembly/analysis/ordered_genes

grep ">" /nbi/Research-Groups/NBI/Cristobal-Uauy/TGACv1_annotation_CS42_ensembl_release/Triticum_aestivum_CS42_TGACv1_scaffold.annotation.gff3.cdna.fa > transcripts_used_for_mapping.txt

sed 's/>//g' <transcripts_used_for_mapping.txt >transcripts_used_for_mapping_replace_arrow.txt

sed 's/gene=//g' <transcripts_used_for_mapping_replace_arrow.txt >transcripts_used_for_mapping_replace_arrow_replace_gene.txt


# First want to find the scaffolds which the transcripts which I used for mapping belong to. This information is contained
# within the gene IDs but they just need some re-shuffling

# gene id looks like:                     TRIAE_CS42_1AL_TGACv1_000001_AA0000010
# scffold id looks like: Triticum_aestivum_CS42_TGACv1_scaffold_000001_1AL

# need to substitute "TRIAE" for "Triticum_aestivum"
# need to remove the transcript id at the end "AA0000010"
# need to re-shuffle the other parts
# easiest method will be to split the transcript ID at the "_" and re-shuffle. Could do this in linux but will do it in excel to save time


cut -f 2 transcripts_used_for_mapping_replace_arrow_replace_gene.txt > gene_ids.txt

# open "gene_ids.txt" in excel and re-shuffle as described above
# saved results as "gene_ids_to_scaffolds.xslx" and exported the final useful column (scaffold IDs) to a text file
# called "scaffold_ids_in_gene_id_order.txt" # need to change EOL to unix in notepad++

# paste together this scaffold info and transcript info:
paste scaffold_ids_in_gene_id_order.txt transcripts_used_for_mapping_replace_arrow_replace_gene.txt > scaffolds_and_transcripts_used_for_mapping.txt

# remove genes on chromosome U which couldn't be unambiguous mapped to IWGSC chromosome arms
grep -v "CS42_U" scaffolds_and_transcripts_used_for_mapping.txt > scaffolds_and_transcripts_used_for_mapping_excl_chrU.txt

# now want to combine the scaffold/transcript info file with the cM position from Christian
# Use join to find the chromosome bin position for each transcript 
# explanation of join:
#`-1 FIELD'
   #  Join on field FIELD (a positive integer) of file 1.

#`-2 FIELD'
 #    Join on field FIELD (a positive integer) of file 2.

#`-j FIELD'
 #    Equivalent to `-1 FIELD -2 FIELD'.

#`-o FIELD-LIST'

# Otherwise, construct each output line according to the format in
# FIELD-LIST.  Each element in FIELD-LIST is either the single
# character `0' or has the form M.N where the file number, M, is `1'
# or `2' and N is a positive field number.

# this assumes files are sorted so sort them first

sort -k1,1 scaffolds_and_transcripts_used_for_mapping_excl_chrU.txt > scaffolds_and_transcripts_used_for_mapping_excl_chrU_sorted.txt

sort -k1,1 /nbi/Research-Groups/NBI/Cristobal-Uauy/PB_collaborations/TGAC_genome/cs42v1-chapman-master/CS42v1_vs_W7984.map.2016-06-08.scaff.chr.cm > CS42v1_vs_W7984.map.2016-06-08.scaff.chr.cm_sorted.txt

join -1 1 -2 1 -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3 scaffolds_and_transcripts_used_for_mapping_excl_chrU_sorted.txt CS42v1_vs_W7984.map.2016-06-08.scaff.chr.cm_sorted.txt > combined_transcripts_and_cM_position.txt

# add a header to this output file

echo -e "scaffold transcript gene biotype confidence chromosome cM" > header.txt
cat header.txt combined_transcripts_and_cM_position.txt > combined_transcripts_and_cM_position_with_header.txt

