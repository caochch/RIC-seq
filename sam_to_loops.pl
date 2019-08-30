#!/usr/bin/perl
die "perl $0 pcp_rep1_interaction.sam want_chr want_start want_end> xxx.bedpe\n" if(@ARGV != 4);
my $contact_sam=shift;
my $want_chr=shift;
my $want_start=shift;
my $want_end=shift;

my %unique;
my %gap_len;
open(CS,$contact_sam) || die;
while(my $frag_a=<CS>){
        if($frag_a=~/^@/){
                next;
        }
        my $frag_b=<CS>;
	my @sub_a=split/\s+/,$frag_a;
	my @sub_b=split/\s+/,$frag_b;
	my @id_info_a=split/_/,$sub_a[0];
	my @id_info_b=split/_/,$sub_b[0];

     	if($id_info_a[0]."_".$id_info_a[1] ne $id_info_b[0]."_".$id_info_b[1]){
        	die "wrong format\n";
      	}
     	else{#same read name
        	my $chr_a=$sub_a[2];
        	my $loci_a=$sub_a[3];
      		my $cigar_a=$sub_a[5];
    		$cigar_a=~/(\d+)M/;
  		my $match_a=$1;
  		my $end_a=$loci_a+$match_a-1;
		$sub_a[1] = $sub_a[1] > 255 ? $sub_a[1]-256 : $sub_a[1];

  		my $chr_b=$sub_b[2];
  		my $loci_b=$sub_b[3];
 		my $cigar_b=$sub_b[5];
 		$cigar_b=~/(\d+)M/;
  		my $match_b=$1;
 		my $end_b=$loci_b+$match_b-1;
		$sub_b[1] = $sub_b[1] > 255 ? $sub_b[1]-256 : $sub_b[1];

		if($id_info_a[2] ne $id_info_b[2] or $id_info_a[2] ne "Plus"){	#same reads; same RNA; same strand
			next;
		}
		
		if($chr_a ne $want_chr or $chr_b ne $want_chr){
			next;
		}

                if($sub_a[0] =~ /^AlignPair/){
			next;
			my @two=($loci_a+$match_a,$loci_b);
			if($two[0] < $two[1]){
				#3'readsAt5'
                        	#print $want_chr,"\t",$two[0],"\t",$two[0]+1,"\t";
				#print $want_chr,"\t",$two[1],"\t",$two[1]+1,"\n";
			}
			else{
				#normal but > 600
				print $want_chr,"\t",$two[1],"\t",$two[1]+1,"\t";
				print $want_chr,"\t",$two[0],"\t",$two[0]+1,"\n";
			}
                }
                elsif($sub_a[0] =~ /^Part/){
			#next;
                        if($end_a >= $loci_b){	#impossible situation
                                print $frag_a,$frag_b;
                                die "this Part-pair is incorrect\n";
                        }
			my $tmp_gap_len=$loci_b-$end_a-1;
			$gap_len{$tmp_gap_len}++;
			#if($tmp_gap_len < 5){
			#	next;
			#}
			print $want_chr,"\t",$end_a,"\t",$end_a+1,"\t";
			print $want_chr,"\t",$loci_b-1,"\t",$loci_b,"\n";
		}
                elsif($sub_a[0] =~ /^Chimeric/){
			#next;
			if($sub_a[1] eq "0"){
				if($end_b > $loci_a){	#maybe some Bugs in chimeric read sets; only a little, so filtered.
					#print $frag_a,$frag_b;
					next;
				}
				print $want_chr,"\t",$loci_b,"\t",$loci_b+1,"\t";
				print $want_chr,"\t",$end_a,"\t",$end_a+1,"\n";
			}
			elsif($sub_a[1] eq "16"){
				if($end_a > $loci_b){	#maybe some Bugs in chimeric read sets; only a little, so filtered.
					#print $frag_a,$frag_b;
					next;
				}
				print $want_chr,"\t",$loci_a,"\t",$loci_a+1,"\t";
				print $want_chr,"\t",$end_b,"\t",$end_b+1,"\n";
			}
			else{
				print $frag_a,$frag_b;
				die "wrong strand\n";
			}
		
		}
		else{
			print $frag_a,$frag_b;
			die "wrong reads class\n";
		}
	}
}

#foreach (sort {$a<=>$b} keys %gap_len){
#	print STDERR $_,"\t",$gap_len{$_},"\n";
#}
