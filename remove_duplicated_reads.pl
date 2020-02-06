#!/usr/bin/perl
die "perl $0 <read1.fq> <read2.fq> <read1.rmDup.fq> <read2.rmDup.fq>\n" if(@ARGV != 4);
my $read1_fq=shift;
my $read2_fq=shift;
my $read1_rmDup_fq=shift;
my $read2_rmDup_fq=shift;


my %unique;
open(RA,$read1_fq) || die;
open(RB,$read2_fq) || die;
open(OA,">$read1_rmDup_fq") || die;
open(OB,">$read2_rmDup_fq") || die;

while(my $id_a=<RA>){
	my $id_b=<RB>;
	my $seq_a=<RA>;
	my $seq_b=<RB>;
	my $symbol_a=<RA>;
	my $symbol_b=<RB>;
	my $qual_a=<RA>;
	my $qual_b=<RB>;
	
	my $whole_reads=$seq_a.$seq_b;
	if($unique{$whole_reads}){
		next;
	}
	else{
		print OA $id_a,$seq_a,$symbol_a,$qual_a;
		print OB $id_b,$seq_b,$symbol_b,$qual_b;
		$unique{$whole_reads}=1;
	}
}
close OA;
close OB;
