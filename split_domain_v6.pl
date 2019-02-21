#!/usr/bin/perl
die "perl $0 gene_dis.list outprefix\n" if(@ARGV != 2);
my $gene_dis_list=shift;
my $out_prefix=shift;

my $minimum_domain_percent=0.05;

my %trans_info;
my %link;
open(GDL,$gene_dis_list) || die;
$/=">";
<GDL>;
while(my $block=<GDL>){
        chomp $block;
        my @lines=split/\n/,$block;
        my $first=shift @lines;
        my ($trans,$gene,$gene_type,$trans_len,$contact_num)=split/\s+/,$first;
	$trans_info{$trans}=$first;
		
        foreach (@lines){
                my ($start,$end)=split/\s+/,$_;
		$link{$trans}{$start."\t".$end}=1;
	}
}

open(DRAW,">$out_prefix.domains_inGeneLoci_forDrawFigures.list") || die;
open(DOMAINS,">$out_prefix.domains_inGeneAtGenomeLoci.bed") || die;

my $potential_cutoff=0.6;
my $full_covered_cutoff=0.4;

foreach my $trans_key (keys %trans_info){
	my ($trans,$gene,$gene_type,$trans_len,$contact_num)=split/\s+/,$trans_info{$trans_key};
	if($contact_num < 500){
		next;
	}
	#if($gene ne "PDE3A" and $gene ne "IMMP2L" and $gene ne "PVT1" and $gene ne "FTX"){
	#if($gene ne "GPC5"){
	#if($gene ne "PDE10A"){
	#if($gene ne "CUX1"){
	#if($gene ne "TSHZ2"){
	#if($gene ne "TRIO"){
	#if($gene ne "HS6ST3"){
	#	next;
	#}

	#if($gene ne "PDE3A" and $gene ne "IMMP2L" and $gene ne "PVT1" and $gene ne "FTX" and $gene ne "GPC5" and $gene ne "PDE10A" and $gene ne "CUX1" and $gene ne "TSHZ2" and $gene ne "TRIO" and $gene ne "HS6ST3"){
	#	next;
	#}


	my @raw_links=keys %{$link{$trans}};

	#initial
	my $max_bin=250;
	my $potential_min_len=int($minimum_domain_percent*$max_bin);
	my $interaction_ruler=1/$max_bin;
	my @boundarys=(0,$max_bin);
	my $total_pairs=$max_bin*($max_bin+1)/2;
	my $last_ratio=0;
	my $best_ratio=0;
	my $inter_pairs_for_best=0;

	#initial links
	my $total_linked_pairs;
	my %unique;
	foreach (@raw_links){
		my ($s,$e)=split/\s+/,$_;
		my $start_bin=int(($s/$trans_len)/$interaction_ruler);
		my $end_bin=int(($e/$trans_len)/$interaction_ruler);
		if($unique{$start_bin}{$end_bin}){
			next;
		}
		push (@links,$_);
		$unique{$start_bin}{$end_bin}=1;
		$total_linked_pairs++;
	}
	
	#testBegin
	#print $#links+1,"\taaalinks\n";
	my $test=1;
	if($test==1){
		open(MAT,">check.matrix") || die;
		print MAT "Bin\t";	
		foreach my $c (0..$max_bin-1){
			print MAT $c,"\t";
		}
		print MAT "\n";
		foreach my $r (0..$max_bin-1){
			print MAT $r,"\t";
			foreach my $c (0..$max_bin-1){
				print MAT $unique{$r}{$c}+0,"\t";
			}
			print MAT "\n";
		}
	}
	#testOver
	
	my %already_boundary=(0=>1,$max_bin=>1);

	my $round_turn;
	my $improvement_for_lastRound;

	my %nopotential_domains;
	$nopotential_domains{250}{"Right"}=1;
	$nopotential_domains{0}{"Left"}=1;
	my %full_domains;

	print ">$trans_key\n";

	while(1){	#iteration
	#foreach (1..7){
		$round_turn++;
		print "Round\t$round_turn\n";

		my %new_Boundary_ratio;
		my %new_Boundary_kept_pairs;
		my %new_Boundary_inter_pairs;
		
		#------------------------------Check,Check,Check------------------------------------#
		#1.Check potential first
		foreach my $i (0..$#boundarys-1){
                        my $no_potential=1;
			my $max_covered_ratio;
			 my $special_potential_min_len;
			if($boundarys[$i] == 0){
				$special_potential_min_len=int($potential_min_len/2);
			}
			else{
				$special_potential_min_len=$potential_min_len;
			}
                        if($boundarys[$i+1]-$boundarys[$i] <= $special_potential_min_len){
                                $no_potential=1;
                        }
                        else{
                                foreach my $j ($boundarys[$i]..$boundarys[$i+1]-$special_potential_min_len){
                                        my $covered_ratio;
                                        foreach my $k ($j..$j+$special_potential_min_len-1){
                                                foreach my $l ($k..$j+$special_potential_min_len-1){
                                                        if($unique{$k}{$l}){
                                                                $covered_ratio++;
                                                        }
                                                }
                                        }
                                        $tmp_total_potential_window=$special_potential_min_len*($special_potential_min_len+1)/2;
                                        $covered_ratio=$covered_ratio/$tmp_total_potential_window;

                                        if($covered_ratio >= $max_covered_ratio){
						$max_covered_ratio=$covered_ratio;
                                        }
                                }
				if($max_covered_ratio >= $potential_cutoff){
					$no_potential=0;
				}
				
                        }
                        if($no_potential){
				$nopotential_domains{$boundarys[$i+1]}{"Left"}=1; #
				$nopotential_domains{$boundarys[$i]}{"Right"}=1;  #
                                next;
                        }
                        else{
                        }
		}
		#2.Check full or not, Second	
		foreach my $i (0..$#boundarys-1){
			my $no_need_to_further_divide;
			foreach my $k ($boundarys[$i]..$boundarys[$i+1]-1){
				foreach my $l ($boundarys[$i]..$boundarys[$i+1]-1){
					if($unique{$k}{$l}){
						$no_need_to_further_divide++;
					}
				}
			}
			my $tmp_total_window=($boundarys[$i+1]-$boundarys[$i])*($boundarys[$i+1]-$boundarys[$i])/2;
			$no_need_to_further_divide=$no_need_to_further_divide/$tmp_total_window;
			if($no_need_to_further_divide >= $full_covered_cutoff){
				$full_domains{$boundarys[$i]}{"Right"}=1; #The right of boundary $i
				next;
			}
		}
		#------------------------------Check Over------------------------------------#
		######################Test#########################
		foreach my $i (0..$#boundarys-1){
			print $boundarys[$i],"--to--",$boundarys[$i+1],"\t";
			if($full_domains{$boundarys[$i]}{"Right"}){
				print "isFull\t";
			}
			else{
				print "CanBeDevide\t";
			}
			if($nopotential_domains{$boundarys[$i]}{"Right"} and $nopotential_domains{$boundarys[$i+1]}{"Left"}){
				print "NoPotential\n";
			}
			else{
				print "HavePotential\n";
			}
		}
		###################Test Over#######################

		#3.find now candidate boundary
		foreach my $i (0..$#boundarys-1){
			if($full_domains{$boundarys[$i]}{"Right"}){
				next;
			}
			if($nopotential_domains{$boundarys[$i]}{"Right"} and $nopotential_domains{$boundarys[$i+1]}{"Left"}){
				next;
			}
			my $new_Boundary_StartBin=$boundarys[$i]+1;
			my $new_Boundary_Endbin=$boundarys[$i+1]-1;
			foreach my $new_Boundary_bin ($new_Boundary_StartBin..$new_Boundary_Endbin){
				#------------------------------Check,Check,Check------------------------------------#
				#1.Chenk both length and potential
				my $left_too_short=1;
				my $No_potential_for_left_new=1;
				my $right_too_short=1;
				my $No_potential_for_right_new=1;
				my $special_potential_min_len;
				if($boundarys[$i] == 0){
					$special_potential_min_len=int($potential_min_len/2);
					$special_potential_min_len=$special_potential_min_len > 5 ? $special_potential_min_len : 5;
				}
				else{
					$special_potential_min_len=$potential_min_len;
				}
				if($new_Boundary_bin-$boundarys[$i] <= $special_potential_min_len){
				}
				else{
					$left_too_short=0;
		                       	foreach my $j ($boundarys[$i]..$new_Boundary_bin-$special_potential_min_len){
	                               	        my $covered_ratio;
	                               	        foreach my $k ($j..$j+$special_potential_min_len-1){
	                               	                foreach my $l ($k..$j+$special_potential_min_len-1){
	                               	                        if($unique{$k}{$l}){
	                               	                                $covered_ratio++;
	                               	                        }
	                               	                }
	                               	        }
	                               	        $tmp_total_potential_window=$special_potential_min_len*($special_potential_min_len+1)/2;
	                               	        $covered_ratio=$covered_ratio/$tmp_total_potential_window;
	                                        if($covered_ratio >= $potential_cutoff){
							$No_potential_for_left_new=0;
							last;
						}
	                   		}
	                   	}	
				if($boundarys[$i+1]-1-$new_Boundary_bin <= $potential_min_len){
				}
				else{
					$right_too_short=0;
                                       	foreach my $j ($new_Boundary_bin..$boundarys[$i+1]-$potential_min_len){
                                               	my $covered_ratio;
                                               	foreach my $k ($j..$j+$potential_min_len-1){
                                                       	foreach my $l ($k..$j+$potential_min_len-1){
                                                      		if($unique{$k}{$l}){
                                                                        $covered_ratio++;
                                                                }
                                                        }
                                                }
                                                $tmp_total_potential_window=$potential_min_len*($potential_min_len+1)/2;
                                                $covered_ratio=$covered_ratio/$tmp_total_potential_window;
                                                if($covered_ratio >= $potential_cutoff){
							$No_potential_for_right_new=0;
                                                        last;
                                             	}
	                            	}
				}
				if($left_too_short or $right_too_short){	#both should be long
					next;
				}
				if($nopotential_domains{$boundarys[$i]}{"Left"} and $nopotential_domains{$boundarys[$i+1]}{"Right"}){		#both side are noPotential;
					if($No_potential_for_left_new and $No_potential_for_right_new){
						next;
					}
				}
				elsif($nopotential_domains{$boundarys[$i]}{"Left"} and !$nopotential_domains{$boundarys[$i+1]}{"Right"}){	#left side are noPotential
					if($No_potential_for_right_new){
						next;
					}
				}
				elsif(!$nopotential_domains{$boundarys[$i]}{"Left"} and $nopotential_domains{$boundarys[$i+1]}{"Right"}){	#right side are noPotential
					if($No_potential_for_left_new){
						next;
					}
				}
				elsif(!$nopotential_domains{$boundarys[$i]}{"Left"} and !$nopotential_domains{$boundarys[$i+1]}{"Right"}){	#both are potential
					if($No_potential_for_left_new or $No_potential_for_right_new){
						next;
					}
				}
				
				
				#------------------------------Check Over------------------------------------#
				##################################Test########################################
				#print $new_Boundary_bin,"\t",$left_too_short,"\t",$right_too_short,"\t",$No_potential_for_left_new,"\t",$No_potential_for_right_new,"\n";
				###############################Test Over######################################

				#2.filter and select new boundary candidate
				my @tmp_boundarys=@boundarys;
				push (@tmp_boundarys,$new_Boundary_bin);
				@tmp_boundarys=sort {$a<=>$b} @tmp_boundarys;

				#whole_pairs
				my $intra_pairs;
				foreach my $j (0..$#tmp_boundarys-1){
					my $tmp_domain_bin_num=$tmp_boundarys[$j+1]-$tmp_boundarys[$j];
					$intra_pairs+=($tmp_domain_bin_num+1)*$tmp_domain_bin_num/2;
				}
				my $inter_pairs=$total_pairs-$intra_pairs;

				#calcluate intra/inter ratio
				my $total_intra;
				my $total_inter;
				foreach my $start_bin (keys %unique){
					my $start_bin_Domain_ID=whichDomain($start_bin,@tmp_boundarys);
					if($start_bin_Domain_ID =~ /Bad/){
						print $start_bin,"\tstart_bin\n";
						die "incorrect bin and domains\n";
					}
					foreach my $end_bin (keys %{$unique{$start_bin}}){
						my $end_bin_Domain_ID=whichDomain($end_bin,@tmp_boundarys);
						if($start_bin_Domain_ID == $end_bin_Domain_ID){
							$total_intra++;
						}
						else{
							$total_inter++;
						}
					}
				}
				my $intra_density=$total_intra/$intra_pairs;
				my $inter_density=$total_inter/$inter_pairs;
				if(!$inter_density){
					$intraInterRatio=10000;	#set as max
				}
				else{
					$intraInterRatio=$intra_density/$inter_density;
				}
				$new_Boundary_ratio{$new_Boundary_bin}=$intraInterRatio;
				$new_Boundary_kept_pairs{$new_Boundary_bin}=$total_intra/$total_linked_pairs;
				$new_Boundary_inter_pairs{$new_Boundary_bin}=$total_inter;
			}
		}

		my @candidate_Boundary=sort {$new_Boundary_ratio{$b} <=> $new_Boundary_ratio{$a}} keys %new_Boundary_ratio;
		my $best_next_boundary=$candidate_Boundary[0];

		print $best_next_boundary,"\tBest\t";
		print $new_Boundary_ratio{$best_next_boundary},"\tBest\t";
		print $new_Boundary_kept_pairs{$best_next_boundary},"\tkept\t";
		print $new_Boundary_inter_pairs{$best_next_boundary},"\ttotalInter\n";
		
		if(!$best_next_boundary){	#no candidate anymore
			last;
		}
		if($new_Boundary_ratio{$best_next_boundary} <= 2){	#ratio too low
			last;
		}
		if($inter_pairs_for_best > 200 and $new_Boundary_ratio{$best_next_boundary} < $best_ratio-0.25 ){
			last;
		}
		if($new_Boundary_kept_pairs{$best_next_boundary} < 0.6){
			last;	
		}
		else{
			push (@boundarys,$best_next_boundary);
			@boundarys=sort {$a<=>$b} @boundarys;
			if($new_Boundary_inter_pairs{$best_next_boundary} > 200 and $new_Boundary_ratio{$best_next_boundary} > $best_ratio){
				$best_ratio=$new_Boundary_ratio{$best_next_boundary};
				$inter_pairs_for_best=$new_Boundary_inter_pairs{$best_next_boundary};
			}
		}
		
		print "Best:\t$best_ratio\t$inter_pairs_for_best\n";

	}

	if($#boundarys < 2){	#at least 2 domains: NEAT1 
		next;
	}

	#4.merge no potential domains
	my %check_final_domains_potential;
        foreach my $i (0..$#boundarys-1){
         	my $no_potential=1;
          	my $max_covered_ratio;
		my $special_potential_min_len;
		if($boundarys[$i] == 0){
			$special_potential_min_len=int($potential_min_len/2);
		}
		else{
			$special_potential_min_len=$potential_min_len;
		}
         	if($boundarys[$i+1]-$boundarys[$i] <= $special_potential_min_len){
       			$no_potential=1;
       		}
       		else{
        		foreach my $j ($boundarys[$i]..$boundarys[$i+1]-$special_potential_min_len){
        	 		my $covered_ratio;
        			foreach my $k ($j..$j+$special_potential_min_len-1){
         				foreach my $l ($k..$j+$special_potential_min_len-1){
         					if($unique{$k}{$l}){
                  					$covered_ratio++;
                				}
         				}
        			}
        			$tmp_total_potential_window=$special_potential_min_len*($special_potential_min_len+1)/2;
             			$covered_ratio=$covered_ratio/$tmp_total_potential_window;
             			if($covered_ratio >= $max_covered_ratio){
            				$max_covered_ratio=$covered_ratio;
 	         		}
 	      		}
			if($max_covered_ratio >= $potential_cutoff){
      				$no_potential=0;
      			}
     		}
  		if($no_potential){
      			$check_final_domains_potential{$boundarys[$i+1]}{"Left"}=1; #
         		$check_final_domains_potential{$boundarys[$i]}{"Right"}=1;  #
		}
	}

	my @final_merged_boundarys;
	foreach my $i (0..$#boundarys){
		if($check_final_domains_potential{$boundarys[$i]}{"Left"} and $check_final_domains_potential{$boundarys[$i]}{"Right"}){
			#if($boundarys[$i]-$boundarys[$i-1] <= $potential_min_len*1.2 and $boundarys[$i+1]-$boundarys[$i] <= $potential_min_len*1.2){
			#}
			#else{
			#	push (@final_merged_boundarys,$boundarys[$i]);
			#}
		}
		else{
			push (@final_merged_boundarys,$boundarys[$i]);
		}
	}

	@boundarys=@final_merged_boundarys;	#replace

	if($#boundarys < 2){    #at least 2 domains: NEAT1
		next;
	}
	print DRAW ">",$trans_info{$trans_key},"\n";
	foreach (0..$#boundarys-1){
		print DRAW $boundarys[$_],"\t",0-$boundarys[$_],"\n";
		print DRAW $boundarys[$_+1],"\t",0-$boundarys[$_],"\n";
		print DRAW $boundarys[$_+1],"\t",0-$boundarys[$_+1],"\n";
	}

	my @trans_genome_loci_info=split/#/,$trans;
	my $numOfdomans=$#boundarys;
	foreach (0..$#boundarys-1){
		my $start_bin=$boundarys[$_];
		my $end_bin=$boundarys[$_+1];
		if($trans_genome_loci_info[4] eq "+"){
			print DOMAINS $trans_genome_loci_info[0],"\t";
			my $start_loci=int($trans_genome_loci_info[1]+$start_bin*$interaction_ruler*$trans_len);
			my $end_loci=int($trans_genome_loci_info[1]+$end_bin*$interaction_ruler*$trans_len);
			print DOMAINS $start_loci,"\t",$end_loci,"\t",$gene,"_",$contact_num,"_Domain_$numOfdomans\t","$_\t+\n";
		}
		elsif($trans_genome_loci_info[4] eq "-"){
			print DOMAINS $trans_genome_loci_info[0],"\t";
			my $start_loci=int($trans_genome_loci_info[2]-$end_bin*$interaction_ruler*$trans_len);
			my $end_loci=int($trans_genome_loci_info[2]-$start_bin*$interaction_ruler*$trans_len);
			print DOMAINS $start_loci,"\t",$end_loci,"\t",$gene,"_",$contact_num,"_Domain_$numOfdomans\t","$_\t-\n";
		}
		else{
			print $trans,"\n";
			die "wrong strand info\n";
		}
	}
}

sub whichDomain{
	my $bin_id=shift;
	my @in_boundary=@_;
	foreach (0..$#in_boundary-1){
		#print $bin_id,"\t",$in_boundary[$_],"\t",$in_boundary[$_+1],"\n";
		if($in_boundary[$_] <= $bin_id and $in_boundary[$_+1] > $bin_id){
			return $_;
		}
	}
	return "Bad";
}
