#!/usr/bin/perl
die "perl $0 <Chimeric.sam> > normal.sam" if(@ARGV != 1);
my $sam=shift;

my %reads_id;
open(IN,$sam) || die;
while(my $line=<IN>){
	if($line=~/^@/){
	}
	else{
		my @sub=split/\s+/,$line;
		$reads_id{$sub[0]}++;
	}
}

close IN;

my $head_or_tail=1;
open(IN,$sam) || die;
while(my $line=<IN>){
	if($line=~/^@/){
		print $line;
	}
	else{
		my @sub=split/\s+/,$line;
		if($reads_id{$sub[0]} != 2){
		}
		else{
			if($head_or_tail){
				print "Head_",$line;
				$head_or_tail=0;
			}
			else{
				print "Tail_",$line;
				$head_or_tail=1;
			}
		}
	}
}
