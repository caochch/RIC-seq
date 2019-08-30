#!/usr/bin/perl
die "perl $0 rRNA.loops want_hs\n" if(@ARGV != 1);
my $in_rRNA_loops=shift;

my $read_id;
my %matrix;
open(IN,$in_rRNA_loops) || die;
while(my $line=<IN>){
	chomp $line;
	my @sub=split/\s+/,$line;
	$read_id++;	
	my $loci_one=$sub[1];
	my $loci_two=$sub[4];
	print $read_id,"\t",16,"\t$sub[0]\t",$loci_one,"\t0\t";
	print 16,"\t$sub[3]\t",$loci_two,"\t1\t50\t50\n";
}

