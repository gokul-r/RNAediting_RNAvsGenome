#use stranded RNA-seq data to determine whether each variant comes from the "positive" or "negative" strand (READ1.bam/READ2.bam are from paired-end sequencing)

use warnings;
use strict;

if (@ARGV != 4) {
	die "need to provide 3 input:Variant list, READ1.bam, READ2.bam and output file name\n";
}
my ($inputfile, $bamfile1, $bamfile2, $outputfile) = ($ARGV[0], $ARGV[1], $ARGV[2], $ARGV[3]);

my $minbasequal = 25;


open (my $INPUT , "<", $inputfile) or die "error opening inputfile: $!\n";
open (my $OUTPUT, ">", $outputfile);

while (<$INPUT>) {
	chomp;
	my $line = $_;
	my @fields = split;
	my $TEMPNAME = join '', $outputfile,'_tmp';
	my ($chr, $position) = ($fields[0], $fields[1]);
	my $bamposition = join '', $chr,':',$position,'-',$position;
	system("/srv/gs1/projects/li/shared-data/software/samtools-0.1.16/samtools view $bamfile2 $bamposition > $TEMPNAME"); #get reads that overlap each candidate site, from READ2.bam (same sequence as original RNA molecule)
	my $editnuc = $fields[4];
	my ($posi, $nega) = (0,0);

	open(my $TMPFILE, "<", $TEMPNAME);
	while (<$TMPFILE>) {
		chomp;
		my @bamfields = split;
		my ($alignment, $readstart, $cigar, $sequence, $qualities) = ($bamfields[1], $bamfields[3], $bamfields[5], $bamfields[9], $bamfields[10]);
		my @sequencebases = split(//,$sequence);
		my @qualscores = split(//,$qualities);
		my ($currentpos, $readpos) = ($readstart,1);
		my $base_readpos;
		my @cigarnums = split(/[MIDNSHP]/, $cigar);
		my @cigarletters = split(/[0-9]+/,$cigar);
		shift @cigarletters;

		for (my $i = 0; $i < @cigarnums; $i++) {
			if ($cigarletters[$i] =~ m/[ISH]/) {
				$readpos = $readpos + $cigarnums[$i];
			}
			elsif ($cigarletters[$i] =~ m/[DN]/) {
				$currentpos = $currentpos + $cigarnums[$i];
			}
			elsif ($cigarletters[$i] =~ m/M/) {
				for (my $j = 0; $j < $cigarnums[$i]; $j++) {
					$base_readpos = $readpos if ($currentpos == $position);
					$currentpos++;
					$readpos++;	
				}	
			}
		}
		next unless $base_readpos;
		my $revstrand = 0;
		$revstrand = 1 if ($alignment & 16);
		if ($revstrand == 0) {
			if (ord($qualscores[$base_readpos-1]) >= $minbasequal+33) {
				$posi++ if ($sequencebases[$base_readpos-1] eq $editnuc);	
			}
		} else {
			if (ord($qualscores[$base_readpos-1]) >= $minbasequal+33) {
				$nega++ if ($sequencebases[$base_readpos-1] eq $editnuc);	
			}
		}	
	}
	system("rm $TEMPNAME");

	system("/srv/gs1/projects/li/shared-data/software/samtools-0.1.16/samtools view $bamfile1 $bamposition > $TEMPNAME"); #get reads that overlap each candidate site, from READ1.bam (antisense sequence as original RNA molecule)
	open(my $TMPFILE2, "<", $TEMPNAME);
	while (<$TMPFILE2>) {
		chomp;
		my @bamfields = split;
		my ($alignment, $readstart, $cigar, $sequence, $qualities) = ($bamfields[1], $bamfields[3], $bamfields[5], $bamfields[9], $bamfields[10]);
		my @sequencebases = split(//,$sequence);
		my @qualscores = split(//,$qualities);
		my ($currentpos, $readpos) = ($readstart,1);
		my $base_readpos;
		my @cigarnums = split(/[MIDNSHP]/, $cigar);
		my @cigarletters = split(/[0-9]+/,$cigar);
		shift @cigarletters;

		for (my $i = 0; $i < @cigarnums; $i++) {
			if ($cigarletters[$i] =~ m/[ISH]/) {
				$readpos = $readpos + $cigarnums[$i];
			}
			elsif ($cigarletters[$i] =~ m/[DN]/) {
				$currentpos = $currentpos + $cigarnums[$i];
			}
			elsif ($cigarletters[$i] =~ m/M/) {
				for (my $j = 0; $j < $cigarnums[$i]; $j++) {
					$base_readpos = $readpos if ($currentpos == $position);
					$currentpos++;
					$readpos++;	
				}	
			}
		}
		next unless $base_readpos;
		my $revstrand = 0;
		$revstrand = 1 if ($alignment & 16);
		if ($revstrand == 0) {
			if (ord($qualscores[$base_readpos-1]) >= $minbasequal+33) {
				$nega++ if ($sequencebases[$base_readpos-1] eq $editnuc);	
			}
		} else {
			if (ord($qualscores[$base_readpos-1]) >= $minbasequal+33) {
				$posi++ if ($sequencebases[$base_readpos-1] eq $editnuc);	
			}
		}	
	}
	system("rm $TEMPNAME");

	if ($posi > 1 || $nega > 1) {
		if ($posi >= $nega) {
			print $OUTPUT "$line\tpos\n";
		} else {
			print $OUTPUT "$line\tneg\n";
		}
	}
}
close $INPUT;	
close $OUTPUT;
