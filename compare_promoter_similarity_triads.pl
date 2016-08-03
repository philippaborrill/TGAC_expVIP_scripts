#!/usr/bin/perl -w

# Philippa.borrill@jic.ac.uk
#
# Aim of script is to calculate percentage similarity between promoters

#### paths and references:
my $path = '/nbi/Research-Groups/NBI/Cristobal-Uauy/';
my $output_dir = "$path/expression_browser/TGAC_assembly/analysis/seedling_stress_analysis/promoter_similarity";
my $bed_dir = "bed_dir";
my $fasta_dir = "fasta_dir";
my $clustalo_dir = "clustalo_dir";

my $TGAC_fasta = "/nbi/Research-Groups/NBI/Cristobal-Uauy/CS_improved/Triticum_aestivum_CS42_TGACv1_all.fa";
# make the fasta index to get the lengths using samtools (samtools faidx)
my $TGAC_fasta_lengths = "/nbi/Research-Groups/NBI/Cristobal-Uauy/CS_improved/Triticum_aestivum_CS42_TGACv1_all.fa.fai";
my $TGAC_gff = "/nbi/Research-Groups/NBI/Cristobal-Uauy/TGACv1_annotation_CS42_ensembl_release/Triticum_aestivum_CS42_TGACv1_scaffold.annotation.gff3";


#my $list_of_triads = "triads_over_90_percent.txt";
#my $list_of_triads = "test_triads_over_90_percent.txt";
#my $list_of_triads = "triads_over_90_percent_1st_2500.txt";
#my $list_of_triads = "triads_over_90_percent_2nd_2500.txt";
#my $list_of_triads = "triads_over_90_percent_3rd_2500.txt";
my $list_of_triads = "triads_over_90_percent_4th_2500.txt";

###############

chdir("$output_dir") or die "couldn't move to output directory";


open (INPUT_FILE, "$list_of_triads") || die "couldn't open the input file $list_of_triads!";
		    while (my $line = <INPUT_FILE>) {
			chomp $line;
 		my @triad = split (/\t/, $line);
#print "\n my triad array is @triad";	

my $group = $triad[0];
my $Agene = $triad[1];
my $Bgene = $triad[2];
my $Dgene = $triad[3];		

#my $feature = "CDS";
my $feature = "mRNA";
my $length = "1000";

# set up header for SLURM control script
my $SLURM_header = <<"SLURM";
#!/bin/bash
#
# SLURM batch script to launch parallel kallisto tasks
#
#SBATCH -p nbi-short # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 1000 # memory pool for all cores
#SBATCH -t 0-00:15 # time (D-HH:MM)
#SBATCH -o $output_dir/slurm_output/promoter_sim.JOBNAME.%N.%j.out # STDOUT
#SBATCH -e $output_dir/slurm_output/promoter_sim.JOBNAME.%N.%j.err # STDERR
#SBATCH -J JOBNAME_promoter_sim
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=philippa.borrill\@jic.ac.uk # send-to address
SLURM



 my $tmp_file = "$output_dir/tmp/compare_promoter_sim.$group.$feature.$length.bp";

  open (SLURM, ">$tmp_file") or die "Couldn't open temp file\n";
  
  $SLURM_header =~ s/JOBNAME/$group/g;
  print SLURM "$SLURM_header\n\n";
  print SLURM "\ncd $output_dir/$bed_dir\n";


# 1) get co-ordinates of A, B and D gene in new file

  print SLURM "\n grep $Agene $TGAC_gff | grep $feature >> $group.$feature.$length.gff";
  print SLURM "\n grep $Bgene $TGAC_gff | grep $feature >> $group.$feature.$length.gff";
  print SLURM "\n grep $Dgene $TGAC_gff | grep $feature >> $group.$feature.$length.gff";

# 2) get 500 bp upstream of A, B and D mRNA (so upstream of transcription start site, excludes 5' UTR) in new file

  print SLURM "\n source bedtools-2.24.0";
  print SLURM "\n bedtools flank -s -i $group.$feature.$length.gff -g $TGAC_fasta_lengths -l $length -r 0 > $group.$length.bp_promoter_$feature.bed";

# 3) extract fasta sequence for A, B and D gene promoter

# need to have the -s so if the original gene was in reverse complement it will reverse complement the sequence
  print SLURM "\n bedtools getfasta -s -fi $TGAC_fasta -bed $group.$length.bp_promoter_$feature.bed -fo $output_dir/$fasta_dir/$group.$length.bp_promoter_$feature.fa";

# 4) run clustalo to compare percentage ID between A, B and D promoter

  print SLURM "\n cd $output_dir/$clustalo_dir";
  print SLURM "\n source clustalo-1.2.0";
  print SLURM "\n clustalo -i $output_dir/$fasta_dir/$group.$length.bp_promoter_$feature.fa --distmat-out=$group.$feature.$length.bp.pim --full --percent-id -o $group.$feature.$length.bp.clustalo.txt --outfmt=clustal";

  close SLURM;
  system("sbatch $tmp_file");
 # unlink $tmp_file;
}


	    close(INPUT_FILE); 


# after loop runs
# 5) make table of perc ID for each group (will need to parse the clustalo output)





