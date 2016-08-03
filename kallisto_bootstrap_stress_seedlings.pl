#!/usr/bin/perl -w

# Philippa.borrill@jic.ac.uk
#
# Aim of script is to run kallisto on RNA-seq for multiple samples to a common reference to calculate expression levels

#### paths and references:
my $path = '/nbi/Research-Groups/NBI/Cristobal-Uauy/';
my $tgac_path_triticum = '/tgac/public/reads/triticum_aestivum/'; 
my $expr_browser_path = "$path/expression_browser/";

my $ref = "/nbi/Research-Groups/NBI/Cristobal-Uauy/TGACv1_annotation_CS42_ensembl_release/Triticum_aestivum_CS42_TGACv1_scaffold.annotation.gff3.cdna.fa";
my $index = "/nbi/Research-Groups/NBI/Cristobal-Uauy/TGACv1_annotation_CS42_ensembl_release/Triticum_aestivum_CS42_TGACv1_scaffold.annotation.gff3.cdna";

# NB make index by kallisto index -i Triticum_aestivum.IWGSC2.26.cdna.all Triticum_aestivum.IWGSC2.26.cdna.all.fa

my $output_dir = "$path/expression_browser/TGAC_assembly/analysis/seedling_stress_analysis/";

### lists of samples (text file containing directory/subdirectory with .fastq to map e.g. each line should look like: ERP004505/ERR392073/ in these subdirectories are the fastq.gz - text file must be in $output_dir):
my $triticum_TGAC_paired_list = 'column_triticum_TGAC_paired_dir_names.txt';
my $triticum_expr_browser_paired_list = 'column_expression_browser_Uauy_paired_dir_names.txt';

#############################

### first do $triticum_TGAC_paired_list
#open the input file and go through the lines one by one so go to each directories where the fastq.gz should be located
chdir("$output_dir") or die "couldn't move to output directory";

#mkdir "fastqc" or  die "couldn't make directory fastqc in $output_dir\n";
open (INPUT_FILE, "$triticum_TGAC_paired_list") || die "couldn't open the input file $triticum_TGAC_paired_list!";
		    while (my $dir = <INPUT_FILE>) {
			chomp $dir;
			
print "\nmy dir: $dir\n";


chdir("$tgac_path_triticum/$dir") or die "couldn't move to tgac specific read directory";

my $SLURM_header = <<"SLURM";
#!/bin/bash
#
# SLURM batch script to launch parallel kallisto tasks
#
#SBATCH -p nbi-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 30000 # memory pool for all cores
#SBATCH -t 2-00:00 # time (D-HH:MM)
#SBATCH -o $output_dir/slurm_output/kallisto.JOBNAME.%N.%j.out # STDOUT
#SBATCH -e $output_dir/slurm_output/kallisto.JOBNAME.%N.%j.err # STDERR
#SBATCH -J JOBNAME_kallisto
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=philippa.borrill\@jic.ac.uk # send-to address
SLURM


##Change lines with # sign to use for single end reads

#get names of .fastq.gz files
opendir(DIR, ".");

@files = grep(/([\w\+\-]+)_1.fastq.gz/,readdir(DIR));
closedir(DIR);

foreach my $file(@files) {
  chomp($file);
print "file is $file\n";

#use this line if paired end reads:
  ($prefix) = $file =~ /([\w\+\-]+)_1.fastq.gz/;
print "prefix is $prefix\n";

#use this line if single end reads:
#  ($prefix) = $file =~ /([\w\+\-]+).fastq.gz/;

   next if ($prefix =~ /^#/);

 my $tmp_file = "$output_dir/tmp/kallisto_paired_triticum_slurm.$prefix";

  open (SLURM, ">$tmp_file") or die "Couldn't open temp file\n";
  my $SLURM_header = $SLURM_header;
  $SLURM_header =~ s/JOBNAME/$prefix/g;
  print SLURM "$SLURM_header\n\n";
  print SLURM "\ncd $tgac_path_triticum/$dir\n";

  print SLURM "source kallisto-0.42.3\n";
  print SLURM "source samtools-0.1.19\n";

#use this line if paired end reads:
	print SLURM "kallisto quant -i $index -o $output_dir/$prefix -b 30 -t 8 $prefix"."_1.fastq.gz $prefix"."_2.fastq.gz \n";

  close SLURM;
  system("sbatch $tmp_file");
 # unlink $tmp_file;
}



## need to close loop which goes through all of the directories in the list
	}
	    close(INPUT_FILE); 



########################


### second do $triticum_expr_browser_paired_list
#open the input file and go through the lines one by one so go to each directories where the fastq.gz should be located
chdir("$output_dir") or die "couldn't move to output directory";

#mkdir "fastqc" or  die "couldn't make directory fastqc in $output_dir\n";
open (INPUT_FILE, "$triticum_expr_browser_paired_list") || die "couldn't open the input file $triticum_expr_browser_paired_list!";
		    while (my $dir = <INPUT_FILE>) {
			chomp $dir;
			
print "\nmy dir: $dir\n";


chdir("$expr_browser_path/$dir") or die "couldn't move to Uauy specific read directory";

my $SLURM_header = <<"SLURM";
#!/bin/bash
#
# SLURM batch script to launch parallel kallisto tasks
#
#SBATCH -p nbi-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 30000 # memory pool for all cores
#SBATCH -t 2-00:00 # time (D-HH:MM)
#SBATCH -o $output_dir/slurm_output/kallisto.JOBNAME.%N.%j.out # STDOUT
#SBATCH -e $output_dir/slurm_output/kallisto.JOBNAME.%N.%j.err # STDERR
#SBATCH -J JOBNAME_kallisto
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=philippa.borrill\@jic.ac.uk # send-to address
SLURM

##Change lines with # sign to use for single end reads

#get names of .fastq.gz files
opendir(DIR, ".");

@files = grep(/([\w\+\-]+)_1.fastq.gz/,readdir(DIR));
closedir(DIR);

foreach my $file(@files) {
  chomp($file);
print "file is $file\n";

#use this line if paired end reads:
  ($prefix) = $file =~ /([\w\+\-]+)_1.fastq.gz/;
print "prefix is $prefix\n";

#use this line if single end reads:
#  ($prefix) = $file =~ /([\w\+\-]+).fastq.gz/;

   next if ($prefix =~ /^#/);

 my $tmp_file = "$output_dir/tmp/kallisto_paired_expr_browser_SLURM.$prefix";

  open (SLURM, ">$tmp_file") or die "Couldn't open temp file\n";
  my $SLURM_header = $SLURM_header;
  $SLURM_header =~ s/JOBNAME/$prefix/g;
  print SLURM "$SLURM_header\n\n";
  print SLURM "\ncd $expr_browser_path/$dir\n";

  print SLURM "source kallisto-0.42.3\n";
  print SLURM "source samtools-0.1.19\n";

#use this line if paired end reads:
	print SLURM "kallisto quant -i $index -o $output_dir/$prefix -b 30 -t 8 $prefix"."_1.fastq.gz $prefix"."_2.fastq.gz \n";

  close SLURM;
  system("sbatch $tmp_file");
 # unlink $tmp_file;
}

## need to close loop which goes through all of the directories in the list
	}
	    close(INPUT_FILE); 




