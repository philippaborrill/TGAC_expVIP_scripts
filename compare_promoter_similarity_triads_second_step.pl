#!/usr/bin/perl -w

# Philippa.borrill@jic.ac.uk
#
# Aim of script is to calculate average percentage similarity between promoters

my $path = '/nbi/Research-Groups/NBI/Cristobal-Uauy/';
my $output_dir = "$path/expression_browser/TGAC_assembly/analysis/seedling_stress_analysis/promoter_similarity";
#my $clustalo_dir = "clustalo_dir";
my $clustalo_dir = "clustalo_dir";
my $output_file = "list_of_triads_with_perc_id_1000_mRNA.csv";

chdir ("$output_dir") or die "couldn't move to output directory";

push @files, `ls -1 $output_dir/$clustalo_dir/*pim`;

foreach my $file(@files) {
  chop($file);
 ($prefix) = $file =~ /([\w\+\-\:\.]+).pim/;

  next if ($prefix =~ /^#/);

	open (FH, "< $file") or die "Can't open $file for read: $!";
	my @lines = <FH>;
	close FH or die "Cannot close $file: $!"; 

#	print $lines[0];
#	print $lines[1];
#	print $lines[2];
#	print $lines[3];

my @A_genomeID = split (/\s+/, $lines[1]);
my @B_genomeID = split (/\s+/, $lines[2]);
my @D_genomeID = split (/\s+/, $lines[3]);

#print $A_genomeID[1];

open (OUTPUT, ">>$output_file") or die "Couldn't open output file\n";
print OUTPUT "\n$prefix, $A_genomeID[2], $A_genomeID[3], $B_genomeID[3]";
#print "$prefix, $A_genomeID[2], $A_genomeID[3], $B_genomeID[3]\n ";


}

close (OUTPUT);

  
