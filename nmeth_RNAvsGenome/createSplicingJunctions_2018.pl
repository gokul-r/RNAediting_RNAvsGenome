#!/usr/bin/perl

use lib qw(BioPerl-1.6.1);
$| = 1;
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use File::Spec;
use File::Basename;
use Data::Dumper;
use File::Basename;
use Bio::DB::Fasta;

################################################################################
################################################################################
# $Revision: $
# Authors: Robert Piskol ( piskol@stanford.edu )
# Last modification $Author: piskol $
#
# get splicing junctions from gene annotations
# remove redundant ones
# save fasta file of these junctions

######## VARIABLES TO CHANGE!! #######

my $INDIR = "annotations/hg19/"; # directory containing gene annotation and genome fasta files

my @bedFiles = qw(hg19_UCSCGene_2011_07_13.txt hg19_refGene_2011_09_27.txt); # Gene annotation files (from UCSC browser?)

my $WINDOWSIZE = 70; # Splice junction window size

my $GENOMEFILE = "hg19all_93chrom.fa"; # genome fasta file

###############################

my $JUNCTIONS;
foreach(@bedFiles){
	print STDERR $_."\n";
	run($INDIR.$_);
}

foreach my $key(sort(keys %$JUNCTIONS)){
	print ">$key\n";
	print $JUNCTIONS->{$key}."\n";
}



sub run{
	my $annotfile = shift;
	my $junctions;
	
	############
	# read genome file
	my $seqs;
	print STDERR "reading genome file...\n";
	my $seqDb =  Bio::DB::Fasta->new($GENOMEFILE);
	#print_dumper($seqDb);
	
	############
	# read annotation file
	my $annot = readAnnotFile($annotfile);
	#print_dumper($annot);
	#exit;
	
	############
	# go through all spliceforms
	foreach my $key(sort(keys %$annot)){
		
		#print $key."\n";
		my @exonStarts = @{$annot->{$key}->{exonStarts}};
		my @exonEnds = @{$annot->{$key}->{exonEnds}};
		my ($chromosome) = $key =~ /(.*):.*/;
		
		next if($#exonStarts == 0);
		#print $chromosome."\n";
		#print_dumper($annot->{$key});
		for(my $i = 1; $i <= $#exonStarts; $i++){
			
			#print "#######\n#\n";
			my ($lSeq,$lSeqPos) = getLseq(\@exonStarts,\@exonEnds, $i-1,$seqDb,$chromosome);
			my ($rSeq,$rSeqPos) = getRseq(\@exonStarts,\@exonEnds, $i,$seqDb,$chromosome);

			my $junctionId =  $chromosome."-".$lSeqPos."-".$rSeqPos;

			print STDERR "$junctionId\r";
			$JUNCTIONS->{$junctionId} = $lSeq.$rSeq;
		}
		
	}
	
	#print_dumper($junctions);
	print "\n";
}

sub getLseq{
	my ($exonStarts, $exonEnds, $i,$seqDb,$chromosome) = @_;
	
	my $remaining = $WINDOWSIZE;
	my $lSeq = "";
	my @lSeqPosArr = ();
	
	while($remaining>0){
		if($i>=0){
			
			my $winStart = $exonEnds->[$i]-$remaining < $exonStarts->[$i] ? $exonStarts->[$i] : $exonEnds->[$i]-$remaining;
			$lSeq = uc($seqDb->get_Seq_by_id($chromosome)->subseq($winStart+1=>$exonEnds->[$i])).$lSeq;
			unshift @lSeqPosArr, ($winStart+1)."-".$exonEnds->[$i];
			$remaining -= ($exonEnds->[$i] - $winStart);
			$i -= 1;
		}
		else{
			return ($lSeq,join("-",@lSeqPosArr));
		}
	}
	return ($lSeq,join("-",@lSeqPosArr));	
}

sub getRseq{
	my ($exonStarts, $exonEnds, $i,$seqDb,$chromosome) = @_;
	
	my $remaining = $WINDOWSIZE;
	my $rSeq = "";
	my @rSeqPosArr = ();
	
	while($remaining>0){
		if($i<=scalar(@{$exonStarts})-1){
			my $winEnd = $exonStarts->[$i]+$remaining > $exonEnds->[$i] ? $exonEnds->[$i] : $exonStarts->[$i]+$remaining;
			$rSeq = $rSeq.uc($seqDb->get_Seq_by_id($chromosome)->subseq($exonStarts->[$i]+1=>$winEnd));
			push @rSeqPosArr, ($exonStarts->[$i]+1)."-".$winEnd;
			$remaining -= ($winEnd - $exonStarts->[$i]);
			$i += 1;
		}
		else{
			return ($rSeq,join("-",@rSeqPosArr));			
		}
	}
	return ($rSeq,join("-",@rSeqPosArr));		
}


sub readAnnotFile{
	my $annotfile = shift;
	my $annot;
	print STDERR $annotfile."\n";
	open(IN, "<".$annotfile);
	my @filecnt = <IN>;
	shift @filecnt;
	#print @filecnt; 
	foreach(@filecnt){
		my $line = $_;
		#print $line."\n";
		my @cnt = split("\t", $line);
		if($cnt[2] !~ /chr.*_.*/){
			$annot->{$cnt[2].":".$cnt[4]."-".$cnt[5].$cnt[1]}->{exonStarts} = [split(",",$cnt[9])];
			$annot->{$cnt[2].":".$cnt[4]."-".$cnt[5].$cnt[1]}->{exonEnds} = [split(",",$cnt[10])];
			$annot->{$cnt[2].":".$cnt[4]."-".$cnt[5].$cnt[1]}->{gene} = $cnt[1];
			#$annot->{$cnt[2].":".$cnt[4]."-".$cnt[5]}->{exonStarts} = [split(",",$cnt[9])];
			#$annot->{$cnt[2].":".$cnt[4]."-".$cnt[5]}->{exonEnds} = [split(",",$cnt[10])];
			#$annot->{$cnt[2].":".$cnt[4]."-".$cnt[5]}->{gene} = $cnt[1];
		}
	}
	
	#print_dumper($annot);
	return $annot;
}

sub print_dumper
{
	my $struct = shift;
	
	my $dumper = Data::Dumper -> new([$struct]);
	my $dump  = $dumper->Dump();
	print $dump;
}

