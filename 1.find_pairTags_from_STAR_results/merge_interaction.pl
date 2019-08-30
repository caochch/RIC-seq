#!/usr/bin/perl
die "perl $0 num_of_interactions_from_part.list interaction_from_pair_mapped_reads_1.sort.sam interaction_from_pair_mapped_reads_2.sort.sam interaction_from_gapped_reads.sam pcp_rep1_read1Chimeric.out.processed.sam pcp_rep1_read2Chimeric.out.processed.sam interaction_merge.sam\n" if(@ARGV != 7);

my $num_of_interactions_from_part_list=shift;
my $paired_read_1_sam=shift;
my $paired_read_2_sam=shift;
my $gapped_sam=shift;
my $chimeric_1=shift;
my $chimeric_2=shift;
my $out_sam=shift;

open(OUTS,">$out_sam") || die;

my $part_from_AlignPair_readA_start=1;
my $part_from_AlignPair_readA_end;
my $part_from_AlignPair_readB_start;
my $part_from_AlignPair_readB_end;
my $part_from_ChimericPair_readA_start;
my $part_from_ChimericPair_readA_end;
my $part_from_ChimericPair_readB_start;
my $part_from_ChimericPair_readB_end;
open(NIP,$num_of_interactions_from_part_list) || die;
while(my $line=<NIP>){
	chomp $line;
	if($line=~/Part_from_Align_Read1:\s+(\d+)/){
		$part_from_AlignPair_readA_end=$1;
		$part_from_AlignPair_readB_start=$part_from_AlignPair_readA_end+1;
	}
	elsif($line=~/Part_from_Align_Read2:\s+(\d+)/){
		$part_from_AlignPair_readB_end=$part_from_AlignPair_readB_start-1+$1;
		$part_from_ChimericPair_readA_start=$part_from_AlignPair_readB_end+1;
	}
	elsif($line=~/Part_from_Chimeric_Read1:\s+(\d+)/){
		$part_from_ChimericPair_readA_end=$part_from_ChimericPair_readA_start-1+$1;
		$part_from_ChimericPair_readB_start=$part_from_ChimericPair_readA_end+1;
	}
	elsif($line=~/Part_from_Chimeric_Read2:\s+(\d+)/){
		$part_from_ChimericPair_readB_end=$part_from_ChimericPair_readB_start-1+$1;
	}
	else{
		die;
	}
}

#Check Begin#
print $part_from_AlignPair_readA_start,"\t",$part_from_AlignPair_readA_end,"\n";
print $part_from_AlignPair_readB_start,"\t",$part_from_AlignPair_readB_end,"\n";
print $part_from_ChimericPair_readA_start,"\t",$part_from_ChimericPair_readA_end,"\n";
print $part_from_ChimericPair_readB_start,"\t",$part_from_ChimericPair_readB_end,"\n";
#Check Done#

my %AlignPair_readA_head_segment;
my %AlignPair_readB_head_segment;
my %Chimeric_readA_junction_segment;
my %Chimeric_readB_junction_segment;
#first store gap information for alignpair reads
my $gapped_pair_id;
open(GS,$gapped_sam) || die;
while(my $frag_a=<GS>){
        if($frag_a=~/^@/){

                next;
        }
        else{
		$gapped_pair_id++;
                my $frag_b=<GS>;
                my $id_a=(split/\s+/,$frag_a)[0];
                my $id_b=(split/\s+/,$frag_b)[0];
                if($id_a ne $id_b){
                        die "wrong format\n";
                }
                else{
			my @sub_a=split/\s+/,$frag_a;
			my @sub_b=split/\s+/,$frag_b;
			my $strand_a=$sub_a[1];
			my $strand_b=$sub_b[1];
			my $loci_a=$sub_a[3];
			my $loci_b=$sub_b[3];
			$strand_a=$strand_a > 255 ? $strand_a-256 : $strand_a;
			$strand_b=$strand_b > 255 ? $strand_b-256 : $strand_b;
			if($strand_a ne $strand_b){
				print $frag_a,$frag_b;
				die "part segments but different strands\n";
			}
			if($gapped_pair_id >= $part_from_AlignPair_readA_start and $gapped_pair_id <= $part_from_AlignPair_readA_end){
				if($strand_a == 0){	#plus strand; save the left most
					if(exists $AlignPair_readA_head_segment{$sub_a[0]}){
						if($loci_a < $AlignPair_readA_head_segment{$sub_a[0]}{1}){
							$AlignPair_readA_head_segment{$sub_a[0]}{0}=$frag_a;
							$AlignPair_readA_head_segment{$sub_a[0]}{1}=$loci_a;
						}
						else{
						}
					}
					else{
						$AlignPair_readA_head_segment{$sub_a[0]}{0}=$frag_a;
						$AlignPair_readA_head_segment{$sub_a[0]}{1}=$loci_a;
					}
				}
				else{			#minus strand; save the right most
					if(exists $AlignPair_readA_head_segment{$sub_a[0]}){
						if($loci_b > $AlignPair_readA_head_segment{$sub_a[0]}{1}){
							$AlignPair_readA_head_segment{$sub_a[0]}{0}=$frag_b;
							$AlignPair_readA_head_segment{$sub_a[0]}{1}=$loci_b;
						}
						else{
						}
					}
					else{
						$AlignPair_readA_head_segment{$sub_a[0]}{0}=$frag_b;
						$AlignPair_readA_head_segment{$sub_a[0]}{1}=$loci_b;
					}
				}
					
			}
			elsif($gapped_pair_id >= $part_from_AlignPair_readB_start and $gapped_pair_id <= $part_from_AlignPair_readB_end){
                                if($strand_a == 0){     #plus strand; save the left most
                                        if(exists $AlignPair_readB_head_segment{$sub_a[0]}){
                                                if($loci_a < $AlignPair_readB_head_segment{$sub_a[0]}{1}){
                                                        $AlignPair_readB_head_segment{$sub_a[0]}{0}=$frag_a;
                                                        $AlignPair_readB_head_segment{$sub_a[0]}{1}=$loci_a;
                                                }
                                                else{
                                                }
                                        }
                                        else{
                                                $AlignPair_readB_head_segment{$sub_a[0]}{0}=$frag_a;
                                                $AlignPair_readB_head_segment{$sub_a[0]}{1}=$loci_a;
                                        }
                                }
                                else{                   #minus strand; save the right most
                                        if(exists $AlignPair_readB_head_segment{$sub_a[0]}){
                                                if($loci_b > $AlignPair_readB_head_segment{$sub_a[0]}{1}){
                                                        $AlignPair_readB_head_segment{$sub_a[0]}{0}=$frag_b;
                                                        $AlignPair_readB_head_segment{$sub_a[0]}{1}=$loci_b;
                                                }
                                                else{
                                                }
                                        }
                                        else{
                                                $AlignPair_readB_head_segment{$sub_a[0]}{0}=$frag_b;
                                                $AlignPair_readB_head_segment{$sub_a[0]}{1}=$loci_b;
                                        }
                                }
			}
			elsif($gapped_pair_id >= $part_from_ChimericPair_readA_start and $gapped_pair_id <= $part_from_ChimericPair_readA_end){
				if(($strand_a == 0 and $sub_a[0] =~ /Head/) or ($strand_a == 16 and $sub_a[0] =~ /Tail/)){	#save the right most
					if(exists $Chimeric_readA_junction_segment{$sub_a[0]}){
						if($loci_b > $Chimeric_readA_junction_segment{$sub_a[0]}{1}){
							$Chimeric_readA_junction_segment{$sub_a[0]}{0}=$frag_b;
							$Chimeric_readA_junction_segment{$sub_a[0]}{1}=$loci_b;
						}
					}
					else{
						$Chimeric_readA_junction_segment{$sub_a[0]}{0}=$frag_b;
						$Chimeric_readA_junction_segment{$sub_a[0]}{1}=$loci_b;
					}
				}
				elsif(($strand_a == 16 and $sub_a[0] =~ /Head/) or ($strand_a == 0 and $sub_a[0] =~ /Tail/)){	#save the left most
					if(exists $Chimeric_readA_junction_segment{$sub_a[0]}){
						if($loci_a < $Chimeric_readA_junction_segment{$sub_a[0]}{1}){
							$Chimeric_readA_junction_segment{$sub_a[0]}{0}=$frag_a;
							$Chimeric_readA_junction_segment{$sub_a[0]}{1}=$loci_a;
						}
					}
					else{
						$Chimeric_readA_junction_segment{$sub_a[0]}{0}=$frag_a;	
						$Chimeric_readA_junction_segment{$sub_a[0]}{1}=$loci_a;
					}
				}
			}
			elsif($gapped_pair_id >= $part_from_ChimericPair_readB_start and $gapped_pair_id <= $part_from_ChimericPair_readB_end){
                                if(($strand_a == 0 and $sub_a[0] =~ /Head/) or ($strand_a == 16 and $sub_a[0] =~ /Tail/)){      #save the right most
                                        if(exists $Chimeric_readB_junction_segment{$sub_a[0]}){
                                                if($loci_b > $Chimeric_readB_junction_segment{$sub_a[0]}{1}){
                                                        $Chimeric_readB_junction_segment{$sub_a[0]}{0}=$frag_b;
                                                        $Chimeric_readB_junction_segment{$sub_a[0]}{1}=$loci_b;
                                                }
                                        }
                                        else{
                                                $Chimeric_readB_junction_segment{$sub_a[0]}{0}=$frag_b;
                                                $Chimeric_readB_junction_segment{$sub_a[0]}{1}=$loci_b;
                                        }
                                }
                                elsif(($strand_a == 16 and $sub_a[0] =~ /Head/) or ($strand_a == 0 and $sub_a[0] =~ /Tail/)){   #save the left most
                                        if(exists $Chimeric_readB_junction_segment{$sub_a[0]}){
                                                if($loci_a < $Chimeric_readB_junction_segment{$sub_a[0]}{1}){
                                                        $Chimeric_readB_junction_segment{$sub_a[0]}{0}=$frag_a;
                                                        $Chimeric_readB_junction_segment{$sub_a[0]}{1}=$loci_a;
                                                }
                                        }
                                        else{
                                                $Chimeric_readB_junction_segment{$sub_a[0]}{0}=$frag_a;
                                                $Chimeric_readB_junction_segment{$sub_a[0]}{1}=$loci_a;
                                        }
                                }

			}
			else{
				print $frag_a,$frag_b;
				die "source for this pair could not be found\n";
			}
                }
        }
}
close GS;


my $interaction_id=0;
open(PSA,$paired_read_1_sam) || die;
open(PSB,$paired_read_2_sam) || die;
while(my $frag_a=<PSA>){
	my $frag_b=<PSB>;
	if($frag_a=~/^@/){
		print OUTS $frag_a;
	}
	else{
		my $id_a=(split/\s+/,$frag_a)[0];
		my $id_b=(split/\s+/,$frag_b)[0];
		my $strand_a=(split/\s+/,$frag_a)[1];
		my $strand_b=(split/\s+/,$frag_b)[1];
		my $strand_a_symbol=read_strand_to_symbol(1,$strand_a);
		my $strand_b_symbol=read_strand_to_symbol(2,$strand_b);
		if($id_a ne $id_b){
			die "wrong format\n";
		}
		else{
			$interaction_id++;
			if(exists $AlignPair_readA_head_segment{$id_a} and exists $AlignPair_readB_head_segment{$id_b}){	#both
				print OUTS "AlignPairSegment_",$interaction_id,"_",$strand_a_symbol,"_",$AlignPair_readA_head_segment{$id_a}{0};
				print OUTS "AlignPairSegment_",$interaction_id,"_",$strand_b_symbol,"_",$AlignPair_readB_head_segment{$id_b}{0};
			}
			elsif(exists $AlignPair_readA_head_segment{$id_a} and !exists $AlignPair_readB_head_segment{$id_b}){
				print OUTS "AlignPairSegment_",$interaction_id,"_",$strand_a_symbol,"_",$AlignPair_readA_head_segment{$id_a}{0};
				print OUTS "AlignPairSegment_",$interaction_id,"_",$strand_b_symbol,"_",$frag_b;
			}
			elsif(!exists $AlignPair_readA_head_segment{$id_a} and exists $AlignPair_readB_head_segment{$id_b}){
				print OUTS "AlignPairSegment_",$interaction_id,"_",$strand_a_symbol,"_",$frag_a;
				print OUTS "AlignPairSegment_",$interaction_id,"_",$strand_b_symbol,"_",$AlignPair_readB_head_segment{$id_b}{0};
			}
			else{
				print OUTS "AlignPairWhole_",$interaction_id,"_",$strand_a_symbol,"_",$frag_a;
				print OUTS "AlignPairWhole_",$interaction_id,"_",$strand_b_symbol,"_",$frag_b;
			}
		}
	}
}


my $gapped_pair_id;
open(GS,$gapped_sam) || die;
while(my $frag_a=<GS>){
	if($frag_a=~/^@/){
		next;
	}
	else{
		my $frag_b=<GS>;
		my $id_a=(split/\s+/,$frag_a)[0];
		my $id_b=(split/\s+/,$frag_b)[0];
		if($id_a ne $id_b){
			die "wrong format\n";
		}
		else{
			$gapped_pair_id++;
			$interaction_id++;
			chomp $frag_a;
			chomp $frag_b;
			$frag_a=~s/\s+$//;
			$frag_b=~s/\s+$//;
	                my $strand_a=(split/\s+/,$frag_a)[1];
	                my $strand_a_symbol;

			if($gapped_pair_id >= $part_from_AlignPair_readA_start and $gapped_pair_id <= $part_from_AlignPair_readA_end){
				$strand_a_symbol=read_strand_to_symbol(1,$strand_a);
			}
			elsif($gapped_pair_id >= $part_from_AlignPair_readB_start and $gapped_pair_id <= $part_from_AlignPair_readB_end){
				$strand_a_symbol=read_strand_to_symbol(2,$strand_a);
			}
			elsif($gapped_pair_id >= $part_from_ChimericPair_readA_start and $gapped_pair_id <= $part_from_ChimericPair_readA_end){
				$strand_a_symbol=read_strand_to_symbol(1,$strand_a);
			}
			elsif($gapped_pair_id >= $part_from_ChimericPair_readB_start and $gapped_pair_id <= $part_from_ChimericPair_readB_end){
				$strand_a_symbol=read_strand_to_symbol(2,$strand_a);
			}

			print OUTS "Part_",$interaction_id,"_",$strand_a_symbol,"_",$frag_a,"\n";
			print OUTS "Part_",$interaction_id,"_",$strand_a_symbol,"_",$frag_b,"\n";
		}
	}
}

open(CSA,$chimeric_1) || die;
while(my $frag_a=<CSA>){
	if($frag_a=~/^@/){
		next;
	}
	else{
		my $frag_b=<CSA>;
		my $id_a=(split/\s+/,$frag_a)[0];
		my $id_b=(split/\s+/,$frag_b)[0];
                my $strand_a=(split/\s+/,$frag_a)[1];
                my $strand_b=(split/\s+/,$frag_b)[1];
                my $strand_a_symbol=read_strand_to_symbol(1,$strand_a);
                my $strand_b_symbol=read_strand_to_symbol(1,$strand_b);
		$interaction_id++;
		if(exists $Chimeric_readA_junction_segment{$id_a} and exists $Chimeric_readA_junction_segment{$id_b}){
			print OUTS "ChimericSegment_",$interaction_id,"_",$strand_a_symbol,"_",$Chimeric_readA_junction_segment{$id_a}{0};
			print OUTS "ChimericSegment_",$interaction_id,"_",$strand_b_symbol,"_",$Chimeric_readA_junction_segment{$id_b}{0};
		}
		elsif(exists $Chimeric_readA_junction_segment{$id_a} and !exists $Chimeric_readA_junction_segment{$id_b}){
			print OUTS "ChimericSegment_",$interaction_id,"_",$strand_a_symbol,"_",$Chimeric_readA_junction_segment{$id_a}{0};
			print OUTS "ChimericSegment_",$interaction_id,"_",$strand_b_symbol,"_",$frag_b;
		}
		elsif(!exists $Chimeric_readA_junction_segment{$id_a} and exists $Chimeric_readA_junction_segment{$id_b}){
			print OUTS "ChimericSegment_",$interaction_id,"_",$strand_a_symbol,"_",$frag_a;
			print OUTS "ChimericSegment_",$interaction_id,"_",$strand_b_symbol,"_",$Chimeric_readA_junction_segment{$id_b}{0};
		}
		else{
			print OUTS "ChimericWhole_",$interaction_id,"_",$strand_a_symbol,"_",$frag_a;
			print OUTS "ChimericWhole_",$interaction_id,"_",$strand_b_symbol,"_",$frag_b;
		}
	}
}

	
open(CSB,$chimeric_2) || die;
while(my $frag_a=<CSB>){
        if($frag_a=~/^@/){
                next;
        }
        else{
                my $frag_b=<CSB>;
                my $id_a=(split/\s+/,$frag_a)[0];
                my $id_b=(split/\s+/,$frag_b)[0];
                my $strand_a=(split/\s+/,$frag_a)[1];
                my $strand_b=(split/\s+/,$frag_b)[1];
                my $strand_a_symbol=read_strand_to_symbol(2,$strand_a);
                my $strand_b_symbol=read_strand_to_symbol(2,$strand_b);
                $interaction_id++;
                if(exists $Chimeric_readB_junction_segment{$id_a} and exists $Chimeric_readB_junction_segment{$id_b}){
                        print OUTS "ChimericSegment_",$interaction_id,"_",$strand_a_symbol,"_",$Chimeric_readB_junction_segment{$id_a}{0};
                        print OUTS "ChimericSegment_",$interaction_id,"_",$strand_b_symbol,"_",$Chimeric_readB_junction_segment{$id_b}{0};
                }
                elsif(exists $Chimeric_readB_junction_segment{$id_a} and !exists $Chimeric_readB_junction_segment{$id_b}){
                        print OUTS "ChimericSegment_",$interaction_id,"_",$strand_a_symbol,"_",$Chimeric_readB_junction_segment{$id_a}{0};
                        print OUTS "ChimericSegment_",$interaction_id,"_",$strand_b_symbol,"_",$frag_b;
                }
                elsif(!exists $Chimeric_readB_junction_segment{$id_a} and exists $Chimeric_readB_junction_segment{$id_b}){
                        print OUTS "ChimericSegment_",$interaction_id,"_",$strand_a_symbol,"_",$frag_a;
                        print OUTS "ChimericSegment_",$interaction_id,"_",$strand_b_symbol,"_",$Chimeric_readB_junction_segment{$id_b}{0};
                }
                else{
                        print OUTS "ChimericWhole_",$interaction_id,"_",$strand_a_symbol,"_",$frag_a;
                        print OUTS "ChimericWhole_",$interaction_id,"_",$strand_b_symbol,"_",$frag_b;
                }

        }
}


sub read_strand_to_symbol{
	my $read_A_or_B=shift;
	my $read_strand=shift;
	if($read_A_or_B == 1){
		if($read_strand == 0 or $read_strand == 256){
			return "Minus";
		}
		elsif($read_strand == 16 or $read_strand == 272){
			return "Plus";
		}
		else{
			return "Ambiguous";
		}
	}
	elsif($read_A_or_B == 2){
		if($read_strand == 0 or $read_strand == 256){
			return "Plus";
		}
		elsif($read_strand == 16 or $read_strand == 272){
			return "Minus";
		}
		else{
			return "Ambiguous";
		}
	}
	else{
		die "only can be 1/2";
	}
}
