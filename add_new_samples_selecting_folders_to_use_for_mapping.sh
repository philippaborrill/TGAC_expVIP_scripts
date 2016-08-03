cd /tgac/public/reads/triticum_aestivum

# use echo to give all sub directors under named directory, pipe the output to a file
# these are for all files with paired end reads
echo SRP056412/*/ >> /nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/expression_browser/scripts_used/add_samples_paired_dir_names.txt

# move to sub directory containing the paired_dir_names.txt
cd /nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/expression_browser/scripts_used/
# manually edit the file in notepad++ to remove the log subdirectories (which don't contain reads therefore we don't want to map)

# use tr to substitute spaces for new lines:
tr ' ' '\n' < add_samples_paired_dir_names.txt > add_samples_column_paired_dir_names_triticum_aestivum.txt

# this results in a file containing all subdirectories which contain fastq.gz, with one subdirectory per line, in the triticum_aestivum folder in /tgac/public/reads/triticum_aestivum
# e.g. one line is: ERP004505/ERR392073/

# made a column_paired_dir_names_puccinia_striiformis.txt file which contains the reads which are from puccinia striiformis experiment (only 1 study on SRA)

#### Do the same thing for the single end reads

echo ERP008767/*/ >> /nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/expression_browser/scripts_used/add_samples_single_dir_names.txt

# move to sub directory containing the single_dir_names_triticum_aestivum.txt
cd /nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/expression_browser/scripts_used/
# manually edit the file in notepad++ to remove the log subdirectories (which don't contain reads therefore we don't want to map)

# use tr to substitute spaces for new lines in resulting paired_dir_names.txt
 tr ' ' '\n' < add_samples_single_dir_names.txt > add_samples_column_single_dir_names_triticum_aestivum.txt


RESULT:

have 2 files which contain the subdirectories to use for mapping:
(copied these to /nbi/group-data/ifs/NBI/Research-Groups/Cristobal-Uauy/expression_browser/alignments)

1) add_samples_column_single_dir_names_triticum_aestivum.txt
2) add_samples_column_paired_dir_names_triticum_aestivum.txt
