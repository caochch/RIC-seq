#!/usr/bin/perl

die "perl $0 interaction_from_pair_mapped_reads_1.sort.sam interaction_from_pair_mapped_reads_2.sort.sam interaction_from_gapped_reads.sam pcp_rep1_read1Chimeric.out.processed.sam pcp_rep1_read2Chimeric.out.processed.sam\n" if(@ARGV != 5);

my $paired_read_1_sam=shift;
my $paired_read_2_sam=shift;
my $gapped_sam=shift;
my $chimeric_1=shift;
my $chimeric_2=shift;

my $pair_1=read_count($paired_read_1_sam);
my $pair_2=read_count($paired_read_2_sam);
my $gap_pair=read_count($gapped_sam);
my $C1=read_count($chimeric_1);
my $C2=read_count($chimeric_2);

if($pair_1 != $pair_2){
	die "pair reads not equal\n";
}
else{
	print "pair\t",$pair_1,"\n";
}

print "gapped\t";
print $gap_pair/2,"\n";

print "C1\t";
print $C1/2,"\n";

print "C2\t";
print $C2/2,"\n";




sub read_count{
	my $sam=shift;
	my $reads_count;
	open(SM,$sam) || die;
	while(my $line=<SM>){
		chomp $line;
		if($line=~/^@/){
			next;
		}
		else{
			$reads_count++;
		}
	}
	return $reads_count;
}
