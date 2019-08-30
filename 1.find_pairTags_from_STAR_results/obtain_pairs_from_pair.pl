#!/usr/bin/perl
die "perl $0 <read_1Aligned.out.sam> <read_2Aligned.out.sam>"  if(@ARGV != 2);
my $align_sam_1=shift;
my $align_sam_2=shift;

#pair id
my $all_pairs_id;

#header
my $sam_head=header($align_sam_1);

#pair from paired reads
my %common_id;
all_id($align_sam_1);
all_id($align_sam_2);
out($align_sam_1,"interaction_from_pair_mapped_reads_1.sam");
out($align_sam_2,"interaction_from_pair_mapped_reads_2.sam");

sub out{
	my $sam=shift;
	my $out=shift;
	open(OUTS,">$out") || die;
	print OUTS $sam_head;
	open(SM,$sam) || die;
	while(my $line=<SM>){
		if($line=~/^@/){
			next;
		}
		else{
			my @sub=split/\s+/,$line;
			if($common_id{$sub[0]}==2){	#happen in both read1 and read2
				print OUTS $line;
			}
		}
	}
	close OUTS;
	close SM;
}


sub header{
	my $sam=shift;
	my $out;
	open(SM,$sam) || die;
	while(my $line=<SM>){
		if($line=~/^@/){
			$out.=$line;
		}
	}
	close SM;
	return $out;
}

sub all_id{
	my $sam=shift;
	open(SM,$sam) || die;
	while(my $line=<SM>){
		chomp $line;
		if($line=~/^@/){
			next;
		}
		else{
			my @sub=split/\s+/,$line;
			$common_id{$sub[0]}++;
		}
	}
	close SM;
}

