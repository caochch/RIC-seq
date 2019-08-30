#!/usr/bin/perl
die "perl $0 <read_1Aligned.out.sam> <read_2Aligned.out.sam> <1Chimeric.out.sam> <2Chimeric.out.sam> <out.sam>"  if(@ARGV != 5);
my $align_sam_1=shift;
my $align_sam_2=shift;
my $chimeric_sam_1=shift;
my $chimeric_sam_2=shift;
my $out_sam=shift;

open(OUTS,">$out_sam") || die;
open(NUM,">num_of_interactions_from_part.list") || die;

#pair id
my $all_pairs_id;

#header
my $sam_head=header($align_sam_1);

#gapped reads
print OUTS $sam_head;
my $tmp_count;

$tmp_count=pair_part($align_sam_1);
print NUM "Part_from_Align_Read1:\t",$tmp_count,"\n";
$tmp_count=pair_part($align_sam_2);
print NUM "Part_from_Align_Read2:\t",$tmp_count,"\n";
$tmp_count=pair_part($chimeric_sam_1);
print NUM "Part_from_Chimeric_Read1:\t",$tmp_count,"\n";
$tmp_count=pair_part($chimeric_sam_2);
print NUM "Part_from_Chimeric_Read2:\t",$tmp_count,"\n";


sub pair_part{
	my $sam=shift;
	my $count=0;
	open(SM,$sam) || die;
	while(my $line=<SM>){
		chomp $line;
		if($line=~/^@/){
			next;
		}
		else{
			my @sub=split/\s+/,$line;
			if($sub[5]=~/[^0-9IMNS]/){#too complicated
				print $sub[5],"\n";
				next;
			}
			if($sub[5]=~/(\d+)N/){	#gapped reads
				my @content;
				my @lens;
				while($sub[5]=~/(\d+)(\w)/g){
					push (@lens,$1);
					push (@content,$2);
				}
				my $most_left=$sub[3];
				my $sum_mid_N=0;
				my $already_match=0;
				foreach my $i (0..$#content){
					my $first_match;
					my $second_match;
					my $first_loci;
					my $second_loci;
					my $mid_N;
					my $mid_info;
					if ($content[$i] eq "M"){
						$first_match=$lens[$i];
						$mid_info.=$lens[$i].$content[$i];
						foreach my $j ($i+1..$#content){
							if($content[$j] eq "M"){
								$second_match=$lens[$j];
								$mid_info.=$lens[$j];
								$mid_info.="M";
								last;
							}
							elsif($content[$j] eq "N"){
								$mid_info.=$lens[$j];
								$mid_info.="N";
								$mid_N=$lens[$j];
							}
							elsif($content[$j] eq "I"){
								$already_match+=$lens[$j];
								$mid_info.=$lens[$j];
								$mid_info.="I";
							}
						}
					}
					else{
						if($content[$i] eq "S" ){
							$already_match+=$lens[$j];
						}
						
					}
					if($second_match){
						my $loci_a=$most_left+$sum_mid_N;
						if($mid_info=~/I/){
						}
						else{
							foreach (0..$#sub){
	                                 	       		if($_ == 5){
	                                 	        	       	print OUTS $first_match,"M\t"
	                                 	       		}
		                                 	       	elsif($_ == 9 or $_ == 10){
        	                         	        	       	print OUTS substr($sub[$_],$already_match,$first_match),"\t";
       		                          	       		}
								elsif($_ == 3){
									print OUTS $loci_a,"\t";
								}
                                 		       		else{
                                 	        		       	print OUTS $sub[$_],"\t";
                                 	       			}
							}
							print OUTS "\n";
						}
						$already_match+=$first_match;
						$sum_mid_N+=$mid_N;
						$sum_mid_N+=$first_match;
						my $loci_b=$most_left+$sum_mid_N;
						if($mid_info=~/I/){
						}
						else{
	                                		foreach (0..$#sub){
	                                        		if($_ == 5){
	                                        	       		print OUTS $second_match,"M\t"
	                                        		}
	                                        		elsif($_ == 9 or $_ == 10){
	                                                		print OUTS substr($sub[$_],$already_match,$second_match),"\t";
	                                        		}
	                                        		elsif($_ == 3){
	                                                		print OUTS $loci_b,"\t";
	                                        		}
	                                        		else{
	                                                		print OUTS $sub[$_],"\t";
	                                        		}
	                                		}
	                                		print OUTS "\n";
							$count++;
						}
					}
				}
			}
		}
	}
	close SM;
	return $count;
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

